#!/usr/bin/env Rscript
# Shiba fixed-effect beta-binomial backend. Invoked once for design and once for
# all count batches, never once per event. Components are fitted separately.
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: beta_binomial_glm.R design|fit config.json")
cfg <- fromJSON(args[2], simplifyVector = FALSE)
out <- cfg$output
write_tsv <- function(x, file, append = FALSE) {
    write.table(x, file, sep = "\t", quote = TRUE, row.names = FALSE,
                col.names = !append, append = append, na = "NA")
}
read_tsv <- function(file) read.delim(file, check.names = FALSE, stringsAsFactors = FALSE,
                                    colClasses = "character", na.strings = "", quote = '"')
# Limit formulas to fixed-effect algebra over metadata columns. In particular,
# formula evaluation must not execute arbitrary R calls from a configuration.
validate_formula <- function(expr, variables) {
    if (is.symbol(expr)) {
        if (!as.character(expr) %in% variables) stop("Unknown formula variable: ", expr)
    } else if (is.call(expr)) {
        if (!as.character(expr[[1]]) %in% c("~", "+", "-", "*", ":", "/", "^", "("))
            stop("Only fixed-effect formula operators are supported")
        for (arg in as.list(expr)[-1]) validate_formula(arg, variables)
    } else if (!is.numeric(expr)) stop("Invalid formula expression")
}

prepare_design <- function() {
    if (!requireNamespace("glmmTMB", quietly = TRUE))
        stop("Install R packages glmmTMB and jsonlite to use beta-binomial regression")
    meta <- read_tsv(cfg$metadata)
    if (anyDuplicated(meta$sample) || anyNA(meta$sample)) stop("Duplicate or missing sample IDs")
    text <- cfg$formula
    if (!startsWith(trimws(text), "~")) text <- paste("~", text)
    f <- as.formula(text, env = baseenv())
    if (length(f) != 2) stop("Supply a right-hand-side formula, without a response")
    validate_formula(f, setdiff(names(meta), "sample"))
    vars <- all.vars(f)
    refs <- cfg$reference_levels
    categorical <- unlist(cfg$categorical, use.names = FALSE)
    continuous <- unlist(cfg$continuous, use.names = FALSE)
    if (length(intersect(categorical, continuous))) stop("A variable cannot be both categorical and continuous")
    if (length(setdiff(c(names(refs), categorical, continuous), vars)))
        stop("Variable type/reference declarations must name formula variables")
    for (v in vars) {
        if (anyNA(meta[[v]]) || any(!nzchar(meta[[v]]))) stop("Missing metadata for ", v)
        if (v %in% categorical || v %in% names(refs)) {
            meta[[v]] <- factor(meta[[v]], levels = sort(unique(meta[[v]])))
        } else if (v %in% continuous) {
            meta[[v]] <- suppressWarnings(as.numeric(meta[[v]]))
        } else {
            meta[[v]] <- type.convert(meta[[v]], as.is = TRUE)
            if (!is.numeric(meta[[v]])) meta[[v]] <- factor(meta[[v]], levels = sort(unique(meta[[v]])))
        }
        if (is.numeric(meta[[v]]) && any(!is.finite(meta[[v]]))) stop("Non-finite numeric metadata: ", v)
        if (v %in% names(refs)) {
            if (!refs[[v]] %in% levels(meta[[v]])) stop("Unknown reference level for ", v)
            meta[[v]] <- relevel(meta[[v]], refs[[v]])
        }
    }
    frame <- model.frame(f, meta, na.action = na.fail)
    X <- model.matrix(f, frame)
    if (!ncol(X) || qr(X)$rank != ncol(X)) stop("Design matrix is rank deficient or empty")
    if (nrow(X) <= ncol(X) + 1L) stop("Insufficient residual degrees of freedom for mean and dispersion")
    contrasts <- list()
    for (coef in unlist(cfg$coefficients, use.names = FALSE)) {
        if (!coef %in% colnames(X)) stop("Unknown coefficient: ", coef, ". Available: ", paste(colnames(X), collapse = ", "))
        cvec <- setNames(rep(0, ncol(X)), colnames(X)); cvec[coef] <- 1
        contrasts[[coef]] <- cvec
    }
    for (name in names(cfg$contrasts)) {
        if (name %in% names(contrasts)) stop("Duplicate contrast name: ", name)
        weights <- cfg$contrasts[[name]]
        if (is.null(names(weights)) || !length(weights) || length(setdiff(names(weights), colnames(X))))
            stop("Contrast must map known coefficient names to weights: ", name)
        cvec <- setNames(rep(0, ncol(X)), colnames(X))
        cvec[names(weights)] <- as.numeric(unlist(weights))
        if (any(!is.finite(cvec)) || !any(cvec != 0)) stop("Invalid contrast: ", name)
        contrasts[[name]] <- cvec
    }
    if (!length(contrasts) && !isTRUE(cfg$design_only)) stop("Specify --coef or --contrast-file")
    prediction <- cfg$prediction
    # Automatic profiles for a selected one- or two-factor binary effect. The
    # other covariates retain the same empirical distribution in every profile.
    if (is.null(prediction) && length(contrasts) == 1L) {
        columns <- which(contrasts[[1]] != 0)
        terms_used <- unique(attr(X, "assign")[columns])
        labels <- attr(terms(f), "term.labels")
        selected <- unique(unlist(strsplit(labels[terms_used[terms_used > 0]], ":", fixed = TRUE)))
        if (length(selected) %in% 1:2 && all(selected %in% vars) &&
            all(vapply(meta[selected], function(x) is.factor(x) && nlevels(x) == 2L, logical(1)))) {
            grid <- expand.grid(lapply(meta[selected], levels), stringsAsFactors = FALSE)
            profiles <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, , drop = FALSE]))
            names(profiles) <- apply(grid, 1, function(x) paste(paste(selected, x, sep = "="), collapse = ";"))
            effects <- list()
            if (length(selected) == 1L) {
                effects[[paste0("delta:", selected)]] <- setNames(as.list(c(-1, 1)), names(profiles))
            } else {
                # expand.grid varies the first factor fastest; differences are
                # second-factor changes at each first-factor level.
                effects[[paste0("delta:", selected[2], "|", selected[1], "=", levels(meta[[selected[1]]])[1])]] <-
                    setNames(as.list(c(-1, 1)), names(profiles)[c(1, 3)])
                effects[[paste0("delta:", selected[2], "|", selected[1], "=", levels(meta[[selected[1]]])[2])]] <-
                    setNames(as.list(c(-1, 1)), names(profiles)[c(2, 4)])
                effects[[paste0("delta_delta:", paste(selected, collapse = ":"))]] <-
                    setNames(as.list(c(1, -1, -1, 1)), names(profiles))
            }
            prediction <- list(profiles = profiles, effects = effects)
        }
    }
    matrices <- list()
    if (!is.null(prediction)) {
        if (is.null(names(prediction$profiles)) || !length(prediction$profiles)) stop("Prediction profiles must be named")
        for (name in names(prediction$profiles)) {
            settings <- prediction$profiles[[name]]
            if (length(setdiff(names(settings), vars))) stop("Unknown prediction variable in ", name)
            newdata <- meta
            for (v in names(settings)) {
                if (length(settings[[v]]) != 1) stop("Profile values must be scalars")
                if (is.factor(meta[[v]])) {
                    if (!settings[[v]] %in% levels(meta[[v]])) stop("Unknown prediction level for ", v)
                    newdata[[v]] <- factor(rep(settings[[v]], nrow(meta)), levels = levels(meta[[v]]))
                } else {
                    value <- suppressWarnings(as.numeric(settings[[v]]))
                    if (!is.finite(value)) stop("Invalid prediction value for ", v)
                    newdata[[v]] <- value
                }
            }
            matrices[[name]] <- model.matrix(f, newdata, contrasts.arg = attr(X, "contrasts"))[, colnames(X), drop = FALSE]
        }
        for (name in names(prediction$effects)) {
            w <- prediction$effects[[name]]
            if (is.null(names(w)) || !length(w) || length(setdiff(names(w), names(matrices))) ||
                any(!is.finite(as.numeric(unlist(w))))) stop("Invalid PSI effect: ", name)
            if (abs(sum(as.numeric(unlist(w)))) > 1e-10) stop("PSI effect weights must sum to zero: ", name)
        }
    }
    write_tsv(data.frame(sample = meta$sample, X, check.names = FALSE), file.path(out, "design_matrix.tsv"))
    ct <- data.frame(coefficient = colnames(X), check.names = FALSE)
    for (name in names(contrasts)) ct[[name]] <- contrasts[[name]]
    write_tsv(ct, file.path(out, "contrasts.tsv"))
    write_json(list(formula = text, coefficients = colnames(X),
                    levels = lapply(meta[vars][vapply(meta[vars], is.factor, logical(1))], levels),
                    prediction = prediction, R = R.version.string,
                    glmmTMB = as.character(packageVersion("glmmTMB")),
                    hypothesis = "component logit-scale contrast equals zero",
                    prediction_population = "all design samples, equal weights"),
               file.path(out, "design.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
    saveRDS(list(X = X, samples = meta$sample, contrasts = contrasts,
                 predictions = matrices, effects = prediction$effects), cfg$design_rds)
}

fit_matrix <- function(X, k, e) {
    dat <- data.frame(k = k, e = e)
    predictors <- paste0("x", seq_len(ncol(X)))
    if (ncol(X)) for (i in seq_len(ncol(X))) dat[[predictors[i]]] <- X[, i]
    f <- if (ncol(X)) reformulate(predictors, response = "cbind(k,e)", intercept = FALSE) else cbind(k,e) ~ 0
    warnings <- character()
    attempt <- function(control) tryCatch(withCallingHandlers(
        glmmTMB::glmmTMB(f, data = dat, family = glmmTMB::betabinomial(link = "logit"),
                        dispformula = ~1, REML = FALSE, control = control),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
        error = function(e) { warnings <<- c(warnings, conditionMessage(e)); NULL })
    good <- function(m) !is.null(m) && isTRUE(m$sdr$pdHess) && m$fit$convergence == 0 &&
        is.finite(as.numeric(logLik(m))) && all(is.finite(glmmTMB::fixef(m)$cond))
    fit <- attempt(glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1500), parallel = 1L))
    if (!good(fit)) fit <- attempt(glmmTMB::glmmTMBControl(optimizer = optim,
                      optArgs = list(method = "BFGS"), optCtrl = list(maxit = 1500), parallel = 1L))
    list(fit = fit, ok = good(fit), warnings = paste(unique(warnings), collapse = " | "))
}

fit_event <- function(dat, design) {
    X <- design$X
    ids <- unique(dat$component_id)
    components <- lapply(ids, function(id) {
        d <- dat[dat$component_id == id, ]
        if (anyDuplicated(d$sample) || !setequal(d$sample, design$samples)) stop("Invalid component sample IDs")
        d[match(design$samples, d$sample), ]
    })
    mask <- Reduce(`&`, lapply(components, function(d) d$inclusion + d$exclusion >= cfg$minimum_reads & d$inclusion + d$exclusion > 0))
    x <- X[mask, , drop = FALSE]
    status <- if (nrow(x) <= ncol(x) + 1L) "insufficient_samples" else if (qr(x)$rank < ncol(x)) "rank_deficient" else "ok"
    stats <- list(); predictions <- list()
    for (i in seq_along(components)) {
        d <- components[[i]]; k <- d$inclusion[mask]; e <- d$exclusion[mask]
        state <- status; full <- NULL
        if (state == "ok" && (all(k == 0) || all(e == 0))) state <- "boundary_response"
        if (state == "ok") {
            full <- fit_matrix(x, k, e)
            if (!full$ok) state <- "full_fit_failed"
        }
        if (state == "ok") {
            beta <- glmmTMB::fixef(full$fit)$cond
            covariance <- as.matrix(vcov(full$fit)$cond)
            phi <- as.numeric(sigma(full$fit)); rho <- 1 / (1 + phi)
            profile <- vapply(design$predictions, function(m) mean(plogis(m %*% beta)), numeric(1))
            values <- c(profile, vapply(design$effects, function(w) sum(profile[names(w)] * as.numeric(unlist(w))), numeric(1)))
            if (length(values)) predictions[[length(predictions) + 1L]] <- data.frame(
                event_type = d$event_type[1], event_id = d$event_id[1], component_id = ids[i],
                kind = c(rep("profile", length(profile)), rep("effect", length(design$effects))),
                name = c(names(profile), names(design$effects)), value = unname(values))
        }
        for (name in names(design$contrasts)) {
            cvec <- design$contrasts[[name]]
            row <- data.frame(event_type = d$event_type[1], event_id = d$event_id[1], component_id = ids[i],
                contrast = name, estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                p_component = NA_real_, lr = NA_real_, rho = NA_real_, phi = NA_real_,
                loglik_full = NA_real_, loglik_null = NA_real_, n_samples = sum(mask),
                full_convergence = if (is.null(full$fit)) NA_integer_ else full$fit$fit$convergence,
                full_pd_hessian = if (is.null(full$fit)) NA else isTRUE(full$fit$sdr$pdHess),
                null_convergence = NA_integer_, null_pd_hessian = NA,
                status = state, diagnostic = if (is.null(full)) "" else full$warnings)
            if (state == "ok") {
                row$estimate <- sum(cvec * beta)
                row$se <- sqrt(as.numeric(t(cvec) %*% covariance %*% cvec))
                row$ci_low <- row$estimate - qnorm(.975) * row$se
                row$ci_high <- row$estimate + qnorm(.975) * row$se
                row$phi <- phi; row$rho <- rho; row$loglik_full <- as.numeric(logLik(full$fit))
                # A basis for c'beta=0 supports arbitrary scalar contrasts, not
                # merely dropping an interaction column or term.
                Q <- qr.Q(qr(matrix(cvec, ncol = 1)), complete = TRUE)
                basis <- Q[, -1, drop = FALSE]
                null <- fit_matrix(x %*% basis, k, e)
                if (!is.null(null$fit)) {
                    row$null_convergence <- null$fit$fit$convergence
                    row$null_pd_hessian <- isTRUE(null$fit$sdr$pdHess)
                }
                row$diagnostic <- paste(Filter(nzchar, c(row$diagnostic, null$warnings)), collapse = " | ")
                if (!null$ok) row$status <- "null_fit_failed" else {
                    row$loglik_null <- as.numeric(logLik(null$fit))
                    lr <- 2 * (row$loglik_full - row$loglik_null)
                    if (!is.finite(lr) || lr < -1e-5 || !is.finite(row$se) || row$se <= 0) {
                        row$status <- "invalid_likelihood_or_covariance"
                    } else {
                        row$lr <- max(0, lr)
                        row$p_component <- pchisq(row$lr, df = 1, lower.tail = FALSE)
                    }
                }
            }
            stats[[length(stats) + 1L]] <- row
        }
    }
    list(stats = do.call(rbind, stats), predictions = if (length(predictions)) do.call(rbind, predictions) else NULL)
}

fit_batches <- function() {
    suppressPackageStartupMessages(library(glmmTMB))
    design <- readRDS(cfg$design_rds)
    first_stats <- TRUE; first_pred <- TRUE
    for (file in unlist(cfg$batches, use.names = FALSE)) {
        dat <- read_tsv(file)
        for (v in c("inclusion", "exclusion")) {
            dat[[v]] <- as.numeric(dat[[v]])
            if (any(!is.finite(dat[[v]]) | dat[[v]] < 0 | dat[[v]] != floor(dat[[v]]))) stop("Invalid integer counts")
        }
        groups <- split(dat, interaction(dat$event_type, dat$event_id, drop = TRUE, lex.order = TRUE))
        message("Fitting ", length(groups), " events from ", basename(file))
        result <- if (cfg$processes > 1L && .Platform$OS.type != "windows")
            parallel::mclapply(groups, fit_event, design = design, mc.cores = cfg$processes, mc.preschedule = TRUE) else
            lapply(groups, fit_event, design = design)
        for (r in result) {
            if (inherits(r, "try-error")) stop("Event worker failed: ", r)
            write_tsv(r$stats, file.path(out, "component_statistics.tsv"), append = !first_stats); first_stats <- FALSE
            if (!is.null(r$predictions)) {
                write_tsv(r$predictions, file.path(out, "component_predictions.tsv"), append = !first_pred); first_pred <- FALSE
            }
        }
    }
}
if (args[1] == "design") prepare_design() else if (args[1] == "fit") fit_batches() else stop("Unknown backend stage")

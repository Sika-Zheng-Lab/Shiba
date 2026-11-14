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
    started <- proc.time()[3]; attempts <- 0L
    dat <- data.frame(k = k, e = e)
    predictors <- paste0("x", seq_len(ncol(X)))
    if (ncol(X)) for (i in seq_len(ncol(X))) dat[[predictors[i]]] <- X[, i]
    f <- if (ncol(X)) reformulate(predictors, response = "cbind(k,e)", intercept = FALSE) else cbind(k,e) ~ 0
    warnings <- character()
    attempt <- function(control) {
        attempts <<- attempts + 1L
        tryCatch(withCallingHandlers(
        glmmTMB::glmmTMB(f, data = dat, family = glmmTMB::betabinomial(link = "logit"),
                        dispformula = ~1, REML = FALSE, control = control),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
        error = function(e) { warnings <<- c(warnings, conditionMessage(e)); NULL })
    }
    good <- function(m) !is.null(m) && isTRUE(m$sdr$pdHess) && m$fit$convergence == 0 &&
        is.finite(as.numeric(logLik(m))) && all(is.finite(glmmTMB::fixef(m)$cond))
    fit <- attempt(glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1500), parallel = 1L))
    if (!good(fit)) fit <- attempt(glmmTMB::glmmTMBControl(optimizer = optim,
                      optArgs = list(method = "BFGS"), optCtrl = list(maxit = 1500), parallel = 1L))
    list(fit = fit, ok = good(fit), warnings = paste(unique(warnings), collapse = " | "),
         attempts = attempts, seconds = unname(proc.time()[3] - started))
}

fit_component <- function(task, design) {
    started <- proc.time()[3]
    components <- list(task$data); ids <- task$data$component_id[1]
    mask <- task$mask; x <- design$X[mask, , drop = FALSE]; status <- task$status
    full_attempts <- null_attempts <- 0L
    full_seconds <- null_seconds <- 0
    stats <- list(); predictions <- list()
    for (i in seq_along(components)) {
        d <- components[[i]]; k <- d$inclusion[mask]; e <- d$exclusion[mask]
        state <- status; full <- NULL
        if (state == "ok" && (all(k == 0) || all(e == 0))) state <- "boundary_response"
        if (state == "ok") {
            full <- fit_matrix(x, k, e)
            full_attempts <- full$attempts; full_seconds <- full$seconds
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
                basis <- design$null_bases[[name]]
                null <- fit_matrix(x %*% basis, k, e)
                null_attempts <- null_attempts + null$attempts
                null_seconds <- null_seconds + null$seconds
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
    list(stats = do.call(rbind, stats), predictions = if (length(predictions)) do.call(rbind, predictions) else NULL,
         timing = data.frame(event_type = task$data$event_type[1], event_id = task$data$event_id[1],
             component_id = ids[1], cache_hit = FALSE, worker_pid = Sys.getpid(),
             elapsed_seconds = unname(proc.time()[3] - started), full_seconds = full_seconds,
             null_seconds = null_seconds, full_attempts = full_attempts, null_attempts = null_attempts,
             retries = max(0L, full_attempts - 1L) + max(0L, null_attempts - length(design$contrasts))))
}

# A run has one immutable design, contrast set, prediction population and optimizer
# configuration. A cache key therefore needs the exact ordered sample mask and the
# included response counts. Hex serialization avoids hash collisions/rounding.
component_key <- function(task) {
    paste(as.character(serialize(list(mask = task$mask,
        k = task$data$inclusion[task$mask], e = task$data$exclusion[task$mask]),
        NULL, version = 2)), collapse = "")
}

prepare_tasks <- function(dat, design) {
    groups <- split(dat, interaction(dat$event_type, dat$event_id, drop = TRUE, lex.order = TRUE))
    tasks <- list()
    for (event in groups) {
        components <- lapply(unique(event$component_id), function(id) {
            d <- event[event$component_id == id, ]
            if (anyDuplicated(d$sample) || !setequal(d$sample, design$samples)) stop("Invalid component sample IDs")
            d[match(design$samples, d$sample), ]
        })
        mask <- Reduce(`&`, lapply(components, function(d) d$inclusion + d$exclusion >= cfg$minimum_reads &
                                   d$inclusion + d$exclusion > 0))
        x <- design$X[mask, , drop = FALSE]
        status <- if (nrow(x) <= ncol(x) + 1L) "insufficient_samples" else if (qr(x)$rank < ncol(x)) "rank_deficient" else "ok"
        for (d in components) tasks[[length(tasks) + 1L]] <- list(data = d, mask = mask, status = status)
    }
    tasks
}

fit_worker <- function(task) fit_component(task, worker_design)

fit_batches <- function() {
    started <- proc.time()[3]
    suppressPackageStartupMessages(library(glmmTMB))
    design <- readRDS(cfg$design_rds)
    design$null_bases <- lapply(design$contrasts, function(cvec)
        qr.Q(qr(matrix(cvec, ncol = 1)), complete = TRUE)[, -1, drop = FALSE])
    worker_design <- design
    cluster <- NULL
    if (cfg$processes > 1L) {
        # Persistent workers avoid forking/loading packages for every input batch.
        cluster <- parallel::makeCluster(cfg$processes,
            type = if (.Platform$OS.type == "windows") "PSOCK" else "FORK")
        on.exit(parallel::stopCluster(cluster), add = TRUE)
        parallel::clusterExport(cluster, c("fit_matrix", "fit_component", "worker_design"), envir = environment())
        parallel::clusterEvalQ(cluster, suppressPackageStartupMessages(library(glmmTMB)))
    }
    cache_limit <- if (is.null(cfg$component_cache_size)) 10000L else as.integer(cfg$component_cache_size)
    cache_values <- vector("list", cache_limit)
    cache_order <- character(cache_limit); cache_position <- 1L
    total <- hits <- attempts <- retries <- 0L
    first_stats <- first_pred <- first_timing <- TRUE
    first_batch <- TRUE
    for (file in unlist(cfg$batches, use.names = FALSE)) {
        batch_start <- proc.time()[3]
        dat <- read_tsv(file)
        for (v in c("inclusion", "exclusion")) {
            dat[[v]] <- as.numeric(dat[[v]])
            if (any(!is.finite(dat[[v]]) | dat[[v]] < 0 | dat[[v]] != floor(dat[[v]]))) stop("Invalid integer counts")
        }
        tasks <- prepare_tasks(dat, design)
        message("Fitting ", length(tasks), " components from ", basename(file))
        keys <- if (cache_limit > 0L) vapply(tasks, component_key, character(1)) else as.character(seq_along(tasks))
        unique_indices <- which(!duplicated(keys))
        cached <- match(keys, cache_order)
        missing <- unique_indices[is.na(cached[unique_indices])]
        # Components are independently fitted only after the event-wide mask is fixed.
        chunk <- max(1L, min(8L, ceiling(length(missing) / (4L * cfg$processes))))
        fitted <- if (!length(missing)) list() else if (is.null(cluster))
            lapply(tasks[missing], fit_component, design = design) else
            parallel::parLapplyLB(cluster, tasks[missing], fit_worker, chunk.size = chunk)
        # Character-vector matching supports arbitrary sample sizes and avoids
        # interning large keys as R environment symbols (which cannot be freed).
        unique_keys <- keys[unique_indices]
        values <- vector("list", length(unique_indices))
        fitted_indices <- match(unique_indices, missing)
        for (j in seq_along(unique_indices)) {
            index <- unique_indices[j]
            values[[j]] <- if (!is.na(cached[index])) cache_values[[cached[index]]] else fitted[[fitted_indices[j]]]
        }
        value_indices <- match(keys, unique_keys)
        reused_flags <- rep(TRUE, length(tasks)); reused_flags[missing] <- FALSE
        statistics <- predictions <- timings <- vector("list", length(tasks))
        for (j in seq_along(tasks)) {
            value <- values[[value_indices[j]]]
            reused <- reused_flags[j]
            ids <- tasks[[j]]$data[1, c("event_type", "event_id", "component_id"), drop = FALSE]
            for (column in names(ids)) {
                value$stats[[column]] <- ids[[column]]
                if (!is.null(value$predictions)) value$predictions[[column]] <- ids[[column]]
                value$timing[[column]] <- ids[[column]]
            }
            if (reused) {
                value$timing$cache_hit <- TRUE
                value$timing$worker_pid <- NA_integer_
                value$timing[, c("elapsed_seconds", "full_seconds", "null_seconds", "full_attempts", "null_attempts", "retries")] <- 0
            }
            statistics[[j]] <- value$stats; predictions[j] <- list(value$predictions); timings[[j]] <- value$timing
        }
        # FIFO bounded cache; entries are compact results, never glmmTMB fit objects.
        if (cache_limit > 0L) for (j in seq_along(missing)) {
            key <- keys[missing[j]]
            cache_values[[cache_position]] <- fitted[[j]]
            cache_order[cache_position] <- key
            cache_position <- cache_position %% cache_limit + 1L
        }
        timing <- do.call(rbind, timings)
        write_tsv(do.call(rbind, statistics), file.path(out, "component_statistics.tsv"), append = !first_stats)
        first_stats <- FALSE
        prediction <- do.call(rbind, predictions)
        if (!is.null(prediction)) {
            write_tsv(prediction, file.path(out, "component_predictions.tsv"), append = !first_pred)
            first_pred <- FALSE
        }
        write_tsv(timing, file.path(out, "component_timings.tsv"), append = !first_timing); first_timing <- FALSE
        batch <- data.frame(batch = basename(file), components = length(tasks), fitted_components = length(missing),
            cache_hits = sum(timing$cache_hit), full_attempts = sum(timing$full_attempts),
            null_attempts = sum(timing$null_attempts), retries = sum(timing$retries),
            elapsed_seconds = unname(proc.time()[3] - batch_start))
        write_tsv(batch, file.path(out, "batch_timings.tsv"), append = !first_batch); first_batch <- FALSE
        message("Completed ", basename(file), ": ", length(tasks), " components, ", batch$cache_hits,
                " cache hits, ", batch$retries, " retries, ", round(batch$elapsed_seconds, 3), " seconds")
        total <- total + length(tasks); hits <- hits + batch$cache_hits
        attempts <- attempts + batch$full_attempts + batch$null_attempts; retries <- retries + batch$retries
    }
    write_json(list(processes = cfg$processes, scheduling = "dynamic components, persistent workers",
        component_cache_size = cache_limit, components = total, cache_hits = hits,
        optimizer_attempts = attempts, retries = retries, elapsed_seconds = unname(proc.time()[3] - started)),
        file.path(out, "backend_performance.json"), auto_unbox = TRUE, pretty = TRUE)
}
if (args[1] == "design") prepare_design() else if (args[1] == "fit") fit_batches() else stop("Unknown backend stage")

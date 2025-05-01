library(DESeq2)
library(data.table)

args <- commandArgs(trailingOnly = T)
experiment_table <- args[1]
counts_table <- args[2]
REFGROUP <- args[3]
ALTGROUP <- args[4]
output <- args[5]

all_colnames <- names(fread(experiment_table, sep = "\t", nrows = 0))
target_cols <- c("sample", "bam", "group")
existing_cols <- intersect(target_cols, all_colnames)
experiment <- fread(
    experiment_table,
    sep = "\t",
    select = existing_cols
)

counts <- read.table(
    counts_table,
    sep = "\t",
    head = T,
    row.names = "gene_name",
    check.names = FALSE,
    stringsAsFactors = FALSE
)
counts_t <- t(counts)

# Merged
merged <- merge(counts_t, experiment, by.x = 0, by.y = "sample")
merged <- merged[merged$group == REFGROUP | merged$group == ALTGROUP, ]
group <- data.frame(con = factor(merged$group))
counts <- merged[, -which (colnames(merged) %in% c("bam", "group"))]
rownames(counts) <- counts$Row.names
counts <- counts[, colnames(counts) != "Row.names"]
counts <- t(counts)
counts <- as.matrix(counts)

# At least six counts
counts <- counts[apply(counts, 1, sum)>6,]

# DEseq
dds <- DESeqDataSetFromMatrix(countData = counts, colData = group, design = ~ con)
dds$con <- relevel(dds$con, ref = REFGROUP)
dds <- DESeq(dds)
res <- results(dds)
res$gene_name <- row.names(res)
res <- res[, c("gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
res <- res[order(res$padj), ]

# save
write.table(
    res,
    file = output,
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
)

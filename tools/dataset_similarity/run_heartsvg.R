args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript run_heartsvg.R <input_tsv> <output_tsv> <p_adj_threshold> <pval_threshold>")
}

input_tsv <- args[[1]]
output_tsv <- args[[2]]
p_adj_threshold <- as.numeric(args[[3]])
pval_threshold <- as.numeric(args[[4]])

suppressPackageStartupMessages({
  library(HEARTSVG)
})

data <- read.delim(input_tsv, check.names = FALSE)
result <- HEARTSVG::heartsvg(data, scale = TRUE)

if (!("pval" %in% colnames(result)) || !("p_adj" %in% colnames(result))) {
  stop("HEARTSVG result must contain columns named 'pval' and 'p_adj'")
}

svg_result <- result[(result$p_adj <= p_adj_threshold) & (result$pval <= pval_threshold), , drop = FALSE]
out <- data.frame(num_svg = nrow(svg_result))
write.table(out, file = output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

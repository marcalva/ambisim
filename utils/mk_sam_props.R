#!/usr/bin/env Rscript

# Create sample prob file from a list of samples.
# Probabilities are uniform across samples
#
# Argument 1: File containing sample IDs. No header or row names.
#
# Argument 2: Output file
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2)
    stop("arguments are required")

sam_fn <- args[1]
out_fn <- args[2]

samples <- readLines(sam_fn)
sam_df <- data.frame(prob = rep(1/length(samples), length(samples)),
                     row.names = samples)

write.table(sam_df, out_fn, row.names = TRUE, col.names = FALSE, quote = FALSE,
            sep = '\t')

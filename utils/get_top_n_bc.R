#!/usr/bin/env Rscript

library(Matrix)

# Argument 1: path to directory with raw counts output, containing files
#  features.tsv.gz, matrix.mtx.gz, barcodes.tsv.gz.
# 
# Argument 2: number of barcodes to output 
#
# Argument 3: path to output file

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2)
    stop("arguments are required")

dir <- args[1]
n_bcs <- as.numeric(args[2])
fn_out <- args[3]

fn = file.path(dir, "features.tsv.gz")
feats = read.table(fn, sep = "\t")

fn = file.path(dir, "matrix.mtx.gz")
mtx = readMM(fn)

fn = file.path(dir, "barcodes.tsv.gz")
bcs = readLines(fn)
bcs = sapply(bcs, function(x) substring(x, 1, 16))
bcs = unname(bcs)

k_rna = feats[,3] == "Gene Expression"
mtx_rna = mtx[k_rna,]
rownames(mtx_rna) = feats[k_rna, 1]
colnames(mtx_rna) = bcs

# run clustering on top 3,000 barcodes
o = order(colSums(mtx_rna), decreasing=TRUE)
top_n = colnames(mtx_rna)[o][seq_len(n_bcs)]

writeLines(top_n, fn_out)


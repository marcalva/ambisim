#!/usr/bin/env Rscript

library(Matrix)
library(Seurat)

# Get cell types, their frequencies, their gene and peak probabilities
# from a real data set
#
# Argument 1: path to directory with raw counts output, containing files
#  features.tsv.gz, matrix.mtx.gz, barcodes.tsv.gz.
# 
# Argument 2: path to barcode file.
#  Consists of barcodes to run clustering on.
#
# Argument 3: output directory

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3)
    stop("arguments are required")

dir <- args[1]
bc_fn <- args[2]
dir_out <- args[3]
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

bcs_flt <- readLines(bc_fn)

fn <- file.path(dir, "features.tsv.gz")
feats <- read.table(fn, sep = "\t")

fn <- file.path(dir, "matrix.mtx.gz")
mtx <- readMM(fn)

fn <- file.path(dir, "barcodes.tsv.gz")
bcs <- readLines(fn)
bcs <- sapply(bcs, function(x) substring(x, 1, 16))
bcs <- unname(bcs)

k_atac <- feats[, 3] == "Peaks"
mtx_atac <- mtx[k_atac, ]
mtx_atac <- as(mtx_atac, "CsparseMatrix")
rownames(mtx_atac) <- feats[k_atac, 1]
colnames(mtx_atac) <- bcs

k_rna <- feats[, 3] == "Gene Expression"
mtx_rna <- mtx[k_rna, ]
mtx_rna <- as(mtx_rna, "CsparseMatrix")
rownames(mtx_rna) <- feats[k_rna, 1]
colnames(mtx_rna) <- bcs

# run clustering on top n barcodes
mtx_rna_n <- mtx_rna[, bcs_flt]
mtx_atac_n <- mtx_atac[, bcs_flt]

set.seed(1)
seur <- CreateSeuratObject(mtx_rna_n, min.cells = 1, verbose = FALSE)
seur <- NormalizeData(seur, normalization.method = "LogNormalize",
                      scale.factor = 1e3, verbose = FALSE)
seur <- FindVariableFeatures(seur, selection.method = "vst",
                             nfeatures = 2000, verbose = FALSE)
seur <- ScaleData(seur, verbose = FALSE)
seur <- RunPCA(seur, verbose = FALSE)
n_dims <- 10
seur <- FindNeighbors(seur, dims = 1:n_dims, verbose = FALSE)
seur <- FindClusters(seur, resolution = c(0.1, 0.2, 0.5), verbose = FALSE)
seur <- RunUMAP(seur, dims = 1:n_dims, verbose = FALSE)

# pcs_dist <- dist(scale(pcs[, 1:n_dims]))
# mhc <- hclust(pcs_dist, method = "ward.D2")
# seur$mhc <- cutree(mhc, k = 8)

# cell type probs
ct_col <- "RNA_snn_res.0.2"
cts <- unique(seur@meta.data[, ct_col])
psum <- list()
gsum <- list()
for (ct in cts){
    ctc <- as.character(ct)
    kct <- which(seur@meta.data[, ct_col] == ct)
    kct <- rownames(seur@meta.data)[kct]
    gsum[[ctc]] <- rowSums(mtx_rna_n[, kct])
    psum[[ctc]] <- rowSums(mtx_atac_n[, kct])
}

gsum <- do.call(cbind, gsum)
gsum <- gsum[, sort(colnames(gsum))]

psum <- do.call(cbind, psum)
psum <- psum[, sort(colnames(psum))]

gsum <- sweep(gsum, 2, colSums(gsum), "/")
psum <- sweep(psum, 2, colSums(psum), "/")

# ambient
amb_bcs <- setdiff(colnames(mtx_rna), bcs_flt)
gsuma <- proportions(rowSums(mtx_rna[, amb_bcs]))
psuma <- proportions(rowSums(mtx_atac[, amb_bcs]))

gsum <- cbind(gsum, gsuma)
psum <- cbind(psum, psuma)

# cell type frequencies
ct_freq <- as.data.frame(proportions(table(seur@meta.data[, ct_col])))
ct_freq <- data.frame(prob = ct_freq[, 2],
                      row.names = ct_freq[, 1])

# splice prob
n_ct <- nlevels(seur@meta.data[, ct_col])
spl_prob <- rep(0.4, n_ct)
spl_prob <- c(spl_prob, 0.6)
spl_prob <- t(as.matrix(spl_prob))

# write output
fn_out <- file.path(dir_out, "ct_freq.txt")
write.table(ct_freq, fn_out, row.names = TRUE,
            col.names = FALSE, quote = FALSE, sep = "\t")

fn_out <- file.path(dir_out, "cell_types.txt")
writeLines(rownames(ct_freq), fn_out)

fn_out <- file.path(dir_out, "gene_probs.txt")
write.table(gsum, fn_out, row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")

fn_out <- file.path(dir_out, "gene_ids.txt")
writeLines(rownames(gsum), fn_out)

fn_out <- file.path(dir_out, "peak_probs.txt")
write.table(psum, fn_out, row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")

fn_out <- file.path(dir_out, "spl_probs.txt")
write.table(spl_prob, fn_out, row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")

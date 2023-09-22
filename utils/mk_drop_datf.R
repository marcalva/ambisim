#!/usr/bin/env Rscript

# generate random BC file
#
# Argument 1: File containing RNA barcodes
#
# Argument 2: File containing ATAC barcodes, corresponding to RNA.
#
# Argument 3: File with sample IDs in the first column and proportions in
#  the second column. No header.
#
# Argument 4: File with cell type indices in the first column and
#  proportions in the second column. No header.
#
# Argument 5: random seed number. 0 for NULL to obtain a different output for
#  each run.
#
# Argument 6: Output file name.

# ==============================================================================
# parameters
# ==============================================================================

n_nuc <- 10e3 # num. of droplets with nuclei
dbl_rate <- 0.1 # proportion of nuclei that are doublets
# sample cell reads with negative binomial
cell_nb_mu <- 5e3
cell_nb_size <- 1.5
cell_reads_min <- 200 # minimum number of reads per celll
cell_reads_max <- 100e3 # max number of reads per cell
# sample ambient reads with negative binomial
ambn_nb_mu <- 2
ambn_nb_size <- .1
# beta parameters for cell contamination
alpha_s1 <- 2
alpha_s2 <- 18
# beta parameters for FRIP in cell
frip_cell_s1 <- 4
frip_cell_s2 <- 6
# beta parameters for FRIP in ambient pool
frip_ambn_s1 <- 1
frip_ambn_s2 <- 9

# ==============================================================================
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    stop("arguments are required")
}

rna_wl_bc_fn <- args[1]
atac_wl_bc_fn <- args[2]
sample_fn <- args[3]
ct_fn <- args[4]
seedn <- as.numeric(args[5])
if (seedn == 0) seedn <- NULL
out_fn <- args[6]

rna_wl_bcs <- readLines(rna_wl_bc_fn)
rna_wl_bcs <- sapply(rna_wl_bcs, function(x) substring(x, 1, 16))
rna_wl_bcs <- unname(rna_wl_bcs)

atac_wl_bcs <- readLines(atac_wl_bc_fn)
atac_wl_bcs <- sapply(atac_wl_bcs, function(x) substring(x, 1, 16))
atac_wl_bcs <- unname(atac_wl_bcs)

samples <- read.table(sample_fn, header = FALSE, row.names = 1)
cts <- read.table(ct_fn, header = FALSE, row.names = 1)
n_cts <- nrow(cts)

wl_bcs <- data.frame("RNA_BC" = rna_wl_bcs, "ATAC_BC" = atac_wl_bcs)

n_amb <- nrow(wl_bcs) - n_nuc

n_dbl <- floor(dbl_rate * n_nuc)
n_sng <- n_nuc - n_dbl

set.seed(seedn)

gen_bc <- function(n, rmin, rmax) {
    # get sample and cell type IDs
    sam <- 0
    ct <- 0
    if (n == 0) {
        sam <- "NA"
        ct <- "NA"
    } else if (n == 1) {
        sam <- sample(rownames(samples), size = 1, prob = samples[, 1])
        ct <- sample(rownames(cts), size = 1, prob = cts[, 1])
    } else {
        # doublets are always heterotypic
        sam <- sample(rownames(samples), size = 1, prob = samples[, 1])
        sam_db <- samples[rownames(samples) != sam, , drop = FALSE]
        sam2 <- sample(rownames(sam_db), size = 1, prob = sam_db[, 1])
        sam <- paste(sam, sam2, sep = ",")
        ct <- sample(rownames(cts), size = 1, prob = cts[, 1])
        cts_db <- cts[rownames(cts) != ct, , drop = FALSE]
        ct2 <- sample(rownames(cts_db), size = 1, prob = cts_db[, 1])
        ct <- paste(ct, ct2, sep = ",")
    }

    # get read numbers
    if (n > 0) {
        n_rna_tot <- rnbinom(1, mu = cell_nb_mu, size = cell_nb_size)
        if (n_rna_tot < cell_reads_min) {
            n_rna_tot <- cell_reads_min
        }
        if (n_rna_tot > cell_reads_max) {
            n_rna_tot <- cell_reads_max
        }
        rna_alpha <- rbeta(1, shape1 = alpha_s1, shape2 = alpha_s2)
        rna_nr_a <- round(n_rna_tot * rna_alpha)
        rna_nr_c <- n_rna_tot - rna_nr_a

        n_atac_tot <- rnbinom(1, mu = cell_nb_mu, size = cell_nb_size)
        if (n_atac_tot < cell_reads_min) {
            n_atac_tot <- cell_reads_min
        }
        if (n_atac_tot > cell_reads_max) {
            n_atac_tot <- cell_reads_max
        }

        atac_frip_c <- rbeta(1, shape1 = frip_cell_s1, shape2 = frip_cell_s2)
        atac_nr_inpk <- round(n_atac_tot * atac_frip_c)
        atac_nr_outpk <- n_atac_tot - atac_nr_inpk

        atac_alpha <- rbeta(1, shape1 = alpha_s1, shape2 = alpha_s2)

        atac_nr_a_ip <- round(atac_nr_inpk * atac_alpha)
        atac_nr_c_ip <- atac_nr_inpk - atac_nr_a_ip
        atac_nr_a_op <- round(atac_nr_outpk * atac_alpha)
        atac_nr_c_op <- atac_nr_outpk - atac_nr_a_op
    } else {
        n_rna_tot <- rnbinom(1, mu = ambn_nb_mu, size = ambn_nb_size)
        rna_nr_a <- n_rna_tot
        rna_nr_c <- 0

        n_atac_tot <- rnbinom(1, mu = ambn_nb_mu, size = ambn_nb_size)

        atac_frip_a <- rbeta(1, shape1 = frip_ambn_s1, shape2 = frip_ambn_s2)
        atac_nr_inpk <- round(n_atac_tot * atac_frip_a)
        atac_nr_outpk <- n_atac_tot - atac_nr_inpk

        atac_alpha <- 1.0

        atac_nr_a_ip <- round(atac_nr_inpk * atac_alpha)
        atac_nr_c_ip <- atac_nr_inpk - atac_nr_a_ip
        atac_nr_a_op <- round(atac_nr_outpk * atac_alpha)
        atac_nr_c_op <- atac_nr_outpk - atac_nr_a_op
    }

    datf <- data.frame(
        n, sam, ct,
        rna_nr_c, rna_nr_a,
        atac_nr_c_ip, atac_nr_a_ip,
        atac_nr_c_op, atac_nr_a_op
    )
    return(datf)
}

sngs <- lapply(seq_len(n_sng), function(x) {
    gen_bc(1, cell_reads_shape, cell_reads_scale)
})
sngs <- do.call(rbind, sngs)

dbls <- lapply(seq_len(n_dbl), function(x) {
    gen_bc(2, cell_reads_shape, cell_reads_scale)
})
dbls <- do.call(rbind, dbls)

ambns <- lapply(seq_len(n_amb), function(x) {
    gen_bc(0, ambn_reads_shape, ambn_reads_scale)
})
ambns <- do.call(rbind, ambns)

datf <- do.call(rbind, list(sngs, dbls, ambns))
datf <- cbind(wl_bcs, datf)

write.table(datf, out_fn,
    row.names = FALSE,
    col.names = TRUE, quote = FALSE, sep = "\t"
)
ss <- sample(seq_len(nrow(datf)), size = 50e3, replace = FALSE)
datf_s <- datf[ss, , drop = FALSE]
write.table(datf_s, paste0(out_fn, ".small"),
    row.names = FALSE,
    col.names = TRUE, quote = FALSE, sep = "\t"
)

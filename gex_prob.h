
#ifndef GEX_PROB_H
#define GEX_PROB_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include <errno.h>
#include "g_list.h"
#include "gtf_anno.h"
#include "str_util.h"
#include "variants.h"
#include "rvd.h"
#include "bc_sim.h"

gex_prob_t *gex_prob_alloc();

void gex_prob_dstry(gex_prob_t *gp);

int gex_prob_load_gene_ids(gex_prob_t *gp, const char *fn);

int gex_prob_load_gtf(gex_prob_t *gp, const char *gtf_fn, int tx_basic);

int gex_prob_set_prob(gex_prob_t *gp, double *rho,
        int rho_ncol, int rho_nrow, 
        double *spl, int spl_len);

int gex_prob_load_prob(gex_prob_t *gp, const char *prob_fn, 
        const char *spl_fn);

// return 0 if match, -1 if no match or error
int gex_prob_check_num_genes(gex_prob_t *gp);

/* sample a gene
 *
 * @param gene Address of pointer. Should be NULL.
 *
 * @return gene index on success, -1 on error. Pointer of address in 
 *  gene is updated to point to the sampled gene_t. Do not free.
 */
int gex_prob_sample_gene(gex_prob_t *gp, uint16_t k, gene_t **gene);

/* sample a transcript
 *
 * Sample a transcript uniformly.
 *
 * @param iso Address of pointer. Should be NULL.
 *
 * @return transcript index on success, -1 on error. Pointer of address 
 * in iso is updated to point to the sampled isoform_t. Do not free.
 */
int gex_prob_sample_tx(gex_prob_t *gp, gene_t *gene, isoform_t **iso);

/* sample mature or pre
 *
 * Sample whether a gene/transcript is mature or unspliced
 *
 * @return 0 for pre-mRNA, 1 for mature mRNA, -1 on error.
 */
int gex_prob_sample_spl(gex_prob_t *gp, uint16_t k);

int isoform_mat_mrna_range(isoform_t *iso, int_ranges_t *ranges);

int isoform_pre_mrna_range(isoform_t *iso, int_ranges_t *ranges);

/* Sample an RNA read
 *
 * @param k cell type index, from 0, ..., K-1, K
 * @param rsam sample index
 * @param fa fasta object
 * @param read_len length of read
 * @param rna_read Pointer to rna_read object to store result in.
 * @return 0 on success, -1 on error.
 */
int gex_sample_read(sc_sim_t *sc_sim, uint16_t k, int rsam, rna_read_t *rna_read);

/*******************************************************************************
 * rna_read_t
 ******************************************************************************/

void rna_read_init(rna_read_t *read);
void rna_read_free(rna_read_t *read);

int rna_read_set_name(rna_read_t *rna_read, il_qname_t *names);

/* Set the seq and qual fields of RNA read 1 (BC+UMI) and 2 (cDNA)
 */
int rna_read_set_seq(rna_read_t *rna_read, const char *bc_name,
        size_t umi_len);

int rna_read_seq_error(sc_sim_t *sc_sim, rna_read_t *rna_read);

/* write out a fastq enetry for read1 (BC+UMI) or read2 (cDNA)
 */
int rna_read_write(rna_read_t *rna_read, BGZF *fs[2]);

/*******************************************************************************
 * test functions
 ******************************************************************************/

int test_sample_gex(gex_prob_t *gp, fa_seq_t *fa, int K, int n_sample, 
        uint32_t read_len, g_var_t *gv);

#endif // GEX_PROB_H


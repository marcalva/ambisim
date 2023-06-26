
#ifndef BC_SIM_H
#define BC_SIM_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include <errno.h>
#include "htslib/bgzf.h"
#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"
#include "g_list.h"
#include "gtf_anno.h"
#include "str_util.h"
#include "variants.h"
#include "rvd.h"

typedef struct {
    gene_anno_t *anno;
    uint16_t k; // number of cell types + ambient pool
    str_map *gene_ix; // gene names corresponding to rows of gene_probs array
    cat_ds_t **gene_probs; // array of gene expr. probs.
    cat_ds_t **spl_probs; // array of prob. of mature mRNA
    uint16_t probs_len;
} gex_prob_t;

typedef struct {
    seq_ranges_t *seq_ranges; // cDNA seq
    char *name[2]; // query names
    char *seq[2]; // seq 1: barcode+UMI, seq 2: cDNA
    mv_t(ui8v) qual[2];
    gene_t *gene;
    isoform_t *iso;
    int mat_rna;
} rna_read_t;

typedef struct {
    uint16_t k; // number of cell types + ambient pool
    iregs_t *peaks; // peaks
    cat_ds_t ** peak_probs; // array of peak probabilities per cell type.
    uint16_t probs_len;
} atac_prob_t;

typedef struct {
    seq_ranges_t *seq_ranges;
    mv_t(ui8v) qual;
    uint8_t peak; // flag indicating whether in peak or not
} atac_read_t;

typedef struct {
    atac_read_t pair[2];
    char *name[3]; // name 1, 2, 3
    char *seq[3];
    mv_t(ui8v) qual[3];
} atac_pair_t;

typedef struct {
    str_map *samples; // sample indexes
    cat_ds_t *sam_probs;
} sam_prob_t;

typedef struct {
    char *rna_barcode;
    char *atac_barcode;
    uint8_t n_cells; // num of cells in droplet, 0, 1, 2, ...
    mv_t(ui8v) samples; // index of sample ID corresponding to VCF
    mv_t(ui8v) cell_types; // index of cell type corresponding to cell_types
    uint32_t rna_nreads_cell; // number of reads from cell/nucleus.
    uint32_t rna_nreads_ambn;
    uint32_t atac_nreads_cell;
    uint32_t atac_nreads_ambn;
    uint32_t atac_npeak_cell;
    uint32_t atac_npeak_ambn;
} bc_sim_t;

mv_declare(bc_simv, bc_sim_t);

typedef struct {
    BGZF *rna_r1_fs, *rna_r2_fs;
    BGZF *atac_r1_fs, *atac_r2_fs, *atac_r3_fs;
} sc_streams_t;

/*! typedef
 * @abstract struct to hold all data for experiment simulation
 * @field barcodes Vector to store barcode information, including num. of reads.
 * @field rna_nreads Total number of RNA reads.
 * @field atac_nreads Total number of ATAC reads.
 * @field gex_prob Store expression probs for each gene.
 * @field atac_prob Store acccessibility probs for each peak.
 * @field K Number of clusters, including ambient cluster (K-1 is number of cell types).
 * @field sr Store vcf/bcf data.
 * @field vcf_hdr Store vcf/bcf header.
 * @field samples Sample string IDs.
 * @field cell_types Cell type string IDs.
 * @field gv Store variant data
 * @field rna_names RNA read names.
 * @field atac_names ATAC read names.
 * @field rna_rd_len Read length fieldeter for RNA.
 * @field atac_rd_len Read length fieldeter for ATAC.
 */
typedef struct {
    mv_t(bc_simv) barcodes;

    uint32_t rna_nreads, atac_nreads; // total
    
    gex_prob_t *gex_prob;
    atac_prob_t *atac_prob;
    sam_prob_t *bg_sam_prob; // background sample probs
    uint8_t has_bg_sam;

    uint16_t K; // number of cell types + ambient (K-1 is num of cell types)

    fa_seq_t *fa;

    bcf_srs_t *sr;
    bcf_hdr_t *vcf_hdr;
    str_map *samples; // sample IDs
    str_map *cell_types; // cell type IDs
    g_var_t *gv;

    str_map *chrms; // overlapping chromosomes

    il_qname_t rna_names;
    il_qname_t atac_names;

    uint32_t rna_rd_len; // length of RNA read
    uint32_t rna_umi_len; // length of RNA read
    uint32_t atac_rd_len; // length of ATAC read
    double seq_error; // parameter for sequencing error

    binom_ds_t rna_err_prob;
    binom_ds_t atac_err_prob;

    uint64_t rna_rd_counter; // counter for number of reads
    uint64_t atac_rd_counter; // counter for number of reads

    sc_streams_t scs; // output streams
    char *out; // output directory path
} sc_sim_t;

/*******************************************************************************
 * bc_sim_t
 ******************************************************************************/

void bc_sim_init(bc_sim_t *bc_sim);
void bc_sim_free(bc_sim_t *bc_sim);

mv_t(bc_simv) bc_sim_read_file(const char *file, str_map *samples,
        str_map *cell_types, int *ret);

/*******************************************************************************
 * sc_streams_t
 ******************************************************************************/

void sc_streams_init(sc_streams_t *scs);
void sc_streams_close(sc_streams_t *scs);
void sc_streams_free(sc_streams_t *scs);

/* Set and open streams for fastq files
 * @return 0 on success, -1 on error
 */
int sc_streams_set(sc_streams_t *scs, const char *dir);

/*******************************************************************************
 * sc_sim_t
 ******************************************************************************/

int sc_sim_init(sc_sim_t *sc_sim);
void sc_sim_free(sc_sim_t *sc_sim);

sc_sim_t *sc_sim_alloc();
void sc_sim_dstry(sc_sim_t *sc_sim);

int sc_sim_set_rd_len(sc_sim_t *sc_sim, uint32_t rna_rd_len,
        uint32_t rna_umi_len, uint32_t atac_rd_len);
int sc_sim_set_seq_error(sc_sim_t *sc_sim, double prob);

/* Load vcf file and variants.
 * Set the samples, if present
 * @return 0 on success, -1 on error.
 */
int sc_sim_load_vars(sc_sim_t *sc_sim, const char *vcf_fn, 
        const char *sample_fn);

int sc_sim_load_cell_types(sc_sim_t *sc_sim, const char *fn);

/* Load GTF file.
 * file is file path to GTF
 * tx_basic is flag (1 to read only basic tag transcripts.
 * @return 0 on success, -1 on error.
 */
int sc_sim_load_gtf(sc_sim_t *sc_sim, const char *file, int tx_basic);

/* Load fasta file.
 * Loads the index and sequence.
 * @return 0 on success, -1 on error.
 */
int sc_sim_load_fa(sc_sim_t *sc_sim, const char *file);

/* Load gex probs
 * rho_file contains the gene expression probabilities, while mmrna_file 
 * contains the splice probabilities. Each has K columns.
 *
 * @return number of cell types, or -1 on error.
 */
int sc_sim_load_gex(sc_sim_t *sc_sim, const char *rho_file, 
        const char *mmrna_file, const char *gene_ids_file);

/* Load background sample probs from file.
 * Returns 0 on success, -1 on error.
 */
int sc_sim_load_bg_sam_prob(sc_sim_t *sc_sim, const char *prob_file);

/* Load atac probs
 * prob_file has the peak probabilities per cell type.
 * peak_file has the peaks in bed format.
 * @return 0 on success, -1 on error.
 */
int sc_sim_load_atac(sc_sim_t *sc_sim, const char *prob_file, 
        const char *peak_file);

int sc_sim_read_file(sc_sim_t *sc_sim, const char *file);

/* Check K after reading rna and atac */
int sc_sim_check_k(sc_sim_t *sc_sim);

/* Intersect chromosomes in gex, atac, fasta, and vcf.
 * Store the results in the chrms field.
 * All modes must be initialized with their chromosome maps set.
 * Return -1 on error, 0 on success.
 */
int sc_sim_intrs_chrms(sc_sim_t *sc_sim);

int sc_sim_sample_gex(sc_sim_t *sc_sim, int n_reads, int read_len);
int sc_sim_sample_atac(sc_sim_t *sc_sim, int n_reads, int read_len);

/* Generate reads in fastq format
 */
int bc_sim_gen_reads(sc_sim_t *sc_sim, size_t bc_i);
int sc_sim_gen_reads(sc_sim_t *sc_sim);

void bc_sim_print(FILE *fs, bc_sim_t *bc_sim);
void sc_sim_print(FILE *fs, sc_sim_t *sc_sim);

#endif // BC_SIM_H 

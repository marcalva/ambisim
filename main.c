
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <getopt.h>
#include <errno.h>
#include "rvd.h"
#include "gex_prob.h"
#include "atac_prob.h"
#include "variants.h"
#include "bc_sim.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"

#define is_overflow(x) (errno == ERANGE)
#define is_underflow(x) (errno == ERANGE)

static void usage(FILE *fp, int exit_status){
    fprintf(fp, 
            "\n"
            "ambisim v0.1: single-cell multiome multiplex and ambient simulator\n"
            "\n"
            "Options:\n"
            "\n"
            "Input files\n"
            "\n"
            "  -v, --vcf                Indexed VCF file.\n" 
            "  -s, --samples            File containing sample IDs in VCF file to simulate genotypes.\n"
            "                           If this option is not given, all samples in VCF will be used.\n"
            "  -c, --cell-types         File containing cell type IDs. One ID per line.\n"
            "  -g, --gtf                Path to GTF file. Required to obtain gene, transcript, and\n"
            "                           exons coordinates.\n"
            "  -f, --fasta              FASTA file for genome, required for pulling the genome sequence.\n"
            "                           Note that the chromosome names must match with the VCF and GTF files.\n"
            "  -d, --drop-file          Path to droplet data table (see README for explanation).\n"
            "  -p, --expr-prob          File path to matrix with gene expression probabilities. Each column gives\n"
            "                           expression probabilities for each gene. The last column is ambient.\n"
            "                           No row names or header.\n"
            "  -G, --gene-ids           Path to file containing gene IDs. Each row contains the gene ID of the\n"
            "                           corresponding row in the 'expr-prob' file. Must have the same number of lines.\n"
            "  -m, --mmrna-prob         Path to file with splicing probabilities per cell type. Gives a numeric matrix\n"
            "                           with one row and K + 1 columns, where columns represent cell types.\n"
            "                           Last column is ambient. No row names or header.\n"
            "  -r, --peak-prob          Path to file with numeric matrix giving open chromatin peak accessibility probabilities.\n"
            "                           Rows represent peaks and columns represent cell types. Last column is ambient.\n"
            "                           No row names or header.\n"
            "  -b, --bg-prob-sam        Optional file containing the sample probabilities of the background. By default,\n"
            "                           the sample for each molecule is taken uniformly. First column contains the sample ID,\n"
            "                           second contains the probabilities. No header.\n"
            "  -o, --peaks              Path to BED file containing peaks. Number of peaks must match number of rows\n"
            "                           in peak probabilities file.\n"
            "  -O, --out                File path to output folder [sim].\n"
            "                           Subdirectories 'RNA' and 'ATAC' will be created, where fastq files will be output.\n"
            "\n"
            "Sequencing options\n"
            "\n"
            "  -e, --seq-error          Probability of a sequencing error [0.01].\n"
            "  -l, --rna-rd-len         Length of RNA read [91].\n"
            "  -U, --rna-umi-len        Length of RNA UMI [12].\n"
            "  -L, --atac-rd-len        Length of each ATAC read pair [50].\n"
            "\n"
            "  -S, --seed               Seed for random number generation.\n"
            "\n"
            "\n"
            "GTF options\n"
            "\n"
            "  -t, --tx-basic           Read only transcripts tagged as 'basic' in the GTF file.\n"
            "\n"
            "  -V, --verbose            Write status on output.\n"
            "      --help               Print this help screen.\n"
            "\n");
    exit(exit_status);
}

int main(int argc, char *argv[]){
    // test_chrm_seq();

    if (argc == 1) usage(stderr, EXIT_FAILURE);


    static const struct option loptions[] =
    {
        {"vcf", required_argument, NULL, 'v'},
        {"samples", required_argument, NULL, 's'},
        {"cell-types", required_argument, NULL, 'c'},
        {"gtf", required_argument, NULL, 'g'},
        {"fasta", required_argument, NULL, 'f'},
        {"drop-file", required_argument, NULL, 'd'},
        {"expr-prob", required_argument, NULL, 'p'},
        {"gene-ids", required_argument, NULL, 'G'},
        {"mmrna-prob", required_argument, NULL, 'm'},
        {"peak-prob", required_argument, NULL, 'r'},
        {"bg-prob-sam", required_argument, NULL, 'b'},
        {"peaks", required_argument, NULL, 'o'},
        {"out", required_argument, NULL, 'O'},
        {"seq-error", required_argument, NULL, 'e'},
        {"rna-rd-len", required_argument, NULL, 'l'},
        {"rna-umi-len", required_argument, NULL, 'U'},
        {"atac-rd-len", required_argument, NULL, 'L'},
        {"seed", required_argument, NULL, 'S'}, 
        {"tx-basic", no_argument, NULL, 't'}, 
        {"verbose", no_argument, NULL, 'V'},
        {"help", no_argument, NULL, '0'},
        {NULL, 0, NULL, 0}
    };

    // variables
    char *vcf_fn = NULL;
    char *sample_fn = NULL;
    char *cell_type_fn = NULL;
    char *gtf_fn = NULL;
    char *fasta_fn = NULL;
    char *bc_dat_fn = NULL;
    char *rho_fn = NULL;
    char *gene_ids_fn = NULL;
    char *mmrna_fn = NULL;
    char *bg_prob_sam_fn = NULL;
    char *peak_prob_fn = NULL;
    char *peaks_fn = NULL;
    char *out_fn = strdup("sim/");
    double seq_error = 0.01;
    uint32_t rna_rd_len = 91;
    uint32_t rna_umi_len = 12;
    uint32_t atac_rd_len = 50;
    int tx_basic = 0;
    int seed = -1;
    int verbose = 0;
    int ret = EXIT_SUCCESS;

    sc_sim_t *sc_sim = sc_sim_alloc();
    if (sc_sim == NULL){
        ret = err_msg(EXIT_FAILURE, 0, "failed to alloc sc_sim: %s",
                strerror(errno));
        goto cleanup;
    }
    il_qname_set_instr(&sc_sim->rna_names, "IL", 1, "FC");
    il_qname_set_instr(&sc_sim->atac_names, "IL", 2, "FC");


    // char *p_end = NULL;
    char *endptr;
    int option_index = 0;
    int cm, cret = 0;
    while ((cm = getopt_long_only(argc, argv, "v:s:c:g:f:d:p:G:m:r:b:o:O:e:l:U:L:S:tV0", 
                    loptions, &option_index)) != -1){
        switch(cm){
            case 'v': vcf_fn = strdup(optarg); 
                      if (vcf_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 's': sample_fn = strdup(optarg); 
                      if (sample_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'c': cell_type_fn = strdup(optarg); 
                      if (cell_type_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'g': gtf_fn = strdup(optarg); 
                      if (gtf_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'f': fasta_fn = strdup(optarg); 
                      if (fasta_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'd': bc_dat_fn = strdup(optarg); 
                      if (bc_dat_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'p': rho_fn = strdup(optarg); 
                      if (rho_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'G': gene_ids_fn = strdup(optarg); 
                      if (gene_ids_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'm':
                      mmrna_fn = strdup(optarg);
                      if (mmrna_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break;
            case 'r':
                      peak_prob_fn = strdup(optarg);
                      if (peak_prob_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break;
            case 'b':
                      bg_prob_sam_fn = strdup(optarg);
                      if (bg_prob_sam_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break;
            case 'o':
                      peaks_fn = strdup(optarg);
                      if (peaks_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break;
            case 'O':
                      free(out_fn);
                      out_fn = strdup(optarg);
                      break;
            case 'e':
                      seq_error = strtod(optarg, &endptr);
                      if ((seq_error < DBL_MIN || seq_error > -DBL_MIN) && 
                              optarg == endptr){
                          ret = err_msg(EXIT_FAILURE, 0, "no conversion of "
                                  "--seq-error %s", optarg);
                          goto cleanup;
                      }
                      if (is_overflow(seq_error)) {
                          ret = err_msg(EXIT_FAILURE, 0, "overflow of "
                                  "--seq-error %s", optarg);
                          goto cleanup;
                      }
                      if (is_underflow(seq_error)) {
                          ret = err_msg(EXIT_FAILURE, 0, "overflow of "
                                  "--seq-error %s", optarg);
                          goto cleanup;
                      }
                      break;
            case 'l': 
                      rna_rd_len = str2uint(optarg, &cret);
                      if (cret < 0){
                          ret = err_msg(EXIT_FAILURE, 0, "failed to convert "
                                  "%s to rna read length", optarg);
                          goto cleanup;
                      }
                      break;
            case 'U': 
                      rna_umi_len = str2uint(optarg, &cret);
                      if (cret < 0){
                          ret = err_msg(EXIT_FAILURE, 0, "failed to convert "
                                  "%s to rna umi length", optarg);
                          goto cleanup;
                      }
                      break;
            case 'L': 
                      atac_rd_len = str2uint(optarg, &cret);
                      if (cret < 0){
                          ret = err_msg(EXIT_FAILURE, 0, "failed to convert "
                                  "%s to atac read length", optarg);
                          goto cleanup;
                      }
                      break;
            case 'S':
                      errno = 0;
                      seed = strtol(optarg, &endptr, 10);
                      if (errno != 0) {
                          perror("strtol");
                          ret = EXIT_FAILURE;
                          goto cleanup;
                      }
                      if (endptr == optarg) {
                          ret = err_msg(EXIT_FAILURE, 0, "failed to convert "
                                  "%s to seed", optarg);
                          goto cleanup;
                      }
                      if (seed < 0) {
                          ret = err_msg(EXIT_FAILURE, 0, "seed='%s'"
                                  "must be >= 0", optarg);
                          goto cleanup;
                      }
                      break;
            case 't':
                      tx_basic = 1;
                      break;
            case 'V': verbose = 1;
                      break;
            case '0':
                      usage(stdout, EXIT_SUCCESS);
                      break;
            default: 
                      usage(stdout, EXIT_FAILURE);
        }
    }

    if (seed >= 0) srand(seed);
    else srand(time(NULL));

    // check required arguments
    if (vcf_fn == NULL) {
        ret = err_msg(EXIT_FAILURE, 0, "vcf file must be provided with --vcf");
        goto cleanup;
    }
    if (gtf_fn == NULL) {
        ret = err_msg(EXIT_FAILURE, 0, "gtf file must be provided with --gtf");
        goto cleanup;
    }
    if (cell_type_fn == NULL) {
        ret = err_msg(EXIT_FAILURE, 0, "cell types must be provided with --cell-types");
        goto cleanup;
    }
    if (fasta_fn == NULL) {
        ret = err_msg(EXIT_FAILURE, 0, "fasta file must be provided with --fasta");
        goto cleanup;
    }
    if (bc_dat_fn == NULL) {
        ret = err_msg(EXIT_FAILURE, 0, "bc data file must be provided with --drop-file");
        goto cleanup;
    }
    if (rho_fn == NULL) {
        ret = err_msg(EXIT_FAILURE, 0, "expression probs must be provided with --expr-prob");
        goto cleanup;
    }
    if (gene_ids_fn == NULL) {
        ret = err_msg(EXIT_FAILURE, 0, "gene IDs must be provided with --gene-ids");
        goto cleanup;
    }
    if (mmrna_fn == NULL) {
        ret = err_msg(EXIT_FAILURE, 0, "splice probs must be provided with --mmrna-fn");
        goto cleanup;
    }
    if (peak_prob_fn == NULL) {
        ret = err_msg(EXIT_FAILURE, 0, "peak probs must be provided with --peak-prob");
        goto cleanup;
    }
    if (peaks_fn == NULL) {
        ret = err_msg(EXIT_FAILURE, 0, "peak bed file must be provided with --peaks");
        goto cleanup;
    }

    if (sc_sim_set_rd_len(sc_sim, rna_rd_len, rna_umi_len, atac_rd_len) < 0) {
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if (sc_sim_set_seq_error(sc_sim, seq_error) < 0) {
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    sc_sim->out = strdup(out_fn);

    if (verbose) log_msg("loading variants");
    if (sc_sim_load_vars(sc_sim, vcf_fn, sample_fn) < 0) {
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (verbose) log_msg("loading GTF");
    if (sc_sim_load_gtf(sc_sim, gtf_fn, tx_basic) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (sc_sim_load_cell_types(sc_sim, cell_type_fn) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (verbose) log_msg("loading fasta");
    if (sc_sim_load_fa(sc_sim, fasta_fn) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (verbose) log_msg("loading gex");
    if (sc_sim_load_gex(sc_sim, rho_fn, mmrna_fn, gene_ids_fn) < 0) {
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (gex_prob_check_num_genes(sc_sim->gex_prob) < 0) {
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (verbose) log_msg("loading atac");
    if (sc_sim_load_atac(sc_sim, peak_prob_fn, peaks_fn) < 0) {
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (bg_prob_sam_fn != NULL) {
        if (verbose) log_msg("loading background sample probabilities");
        if (sc_sim_load_bg_sam_prob(sc_sim, bg_prob_sam_fn) < 0) {
            ret = EXIT_FAILURE;
            goto cleanup;
        }
    }

    if (sc_sim_check_k(sc_sim) < 0) {
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (sc_sim_intrs_chrms(sc_sim) < 0) {
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (sc_sim->chrms->n < 1) {
        err_msg(-1, 0, "no overlapping chromosomes, make sure naming conventions "
                "are consistent across input files");
    }
    if (verbose) {
        log_msg("%i overlapping chromosomes:", sc_sim->chrms->n);
        int i;
        for (i = 0; i < sc_sim->chrms->n; ++i)
            log_msg("\t%s", str_map_str(sc_sim->chrms, i));
    }

    if (verbose) log_msg("loading droplet data");
    if (sc_sim_read_file(sc_sim, bc_dat_fn) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if (verbose) {
        log_msg("generating %i RNA reads and %i ATAC reads",
                sc_sim->rna_nreads, sc_sim->atac_nreads);
    }

    sc_sim->rna_names.n_reads = sc_sim->rna_nreads;
    il_qname_set_step(&sc_sim->rna_names);
    sc_sim->atac_names.n_reads = sc_sim->atac_nreads;
    il_qname_set_step(&sc_sim->atac_names);
    
    if (sc_sim_gen_reads(sc_sim) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (verbose) log_msg("done");

cleanup:

    // free file names
    free(vcf_fn);
    free(sample_fn);
    free(cell_type_fn);
    free(gtf_fn);
    free(fasta_fn);
    free(bc_dat_fn);
    free(rho_fn);
    free(gene_ids_fn);
    free(peak_prob_fn);
    free(bg_prob_sam_fn);
    free(peaks_fn);
    free(mmrna_fn);
    free(out_fn);

    sc_sim_dstry(sc_sim);

    return ret;
}

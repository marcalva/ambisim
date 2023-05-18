
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include <errno.h>
#include <time.h>
#include "gtf_anno.h"
#include "str_util.h"
#include "gex_prob.h"
#include "array_util.h"
#include "variants.h"

gex_prob_t *gex_prob_alloc(){
    gex_prob_t *gp = (gex_prob_t *)malloc(sizeof(gex_prob_t));
    if (gp == NULL){
        err_msg(-1, 0, "gex_prob_alloc: %s", strerror(errno));
        return NULL;
    }
    gp->anno = NULL;
    gp->k = 0;
    gp->gene_ix = NULL;
    gp->gene_probs = NULL;
    gp->spl_probs = NULL;
    gp->probs_len = 0;
    return gp;
}

void gex_prob_dstry(gex_prob_t *gp){
    if (gp == NULL) return;
    gene_anno_dstry(gp->anno);
    destroy_str_map(gp->gene_ix);
    uint16_t i;
    for (i = 0; i < gp->probs_len; ++i){
        cat_ds_dstry(gp->gene_probs[i]);
    }
    free(gp->gene_probs);
    for (i = 0; i < gp->probs_len; ++i){
        cat_ds_dstry(gp->spl_probs[i]);
    }
    free(gp->spl_probs);
    free(gp);
}

int gex_prob_load_gtf(gex_prob_t *gp, const char *gtf_fn, int tx_basic){
    if (gp == NULL || gtf_fn == NULL) return -1;

    gp->anno = read_from_gtf(gtf_fn, tx_basic);
    if (gp->anno == NULL)
        return err_msg(-1, 0, "gex_prob_load_gtf: failed to read GTF file");

    // fill in coding sequence lengths
    if (gene_anno_cds_len(gp->anno) < 0) return -1;

    return 0;
}

int gex_prob_load_gene_ids(gex_prob_t *gp, const char *fn) {
    if (gp == NULL || fn == NULL)
        return err_msg(-1, 0, "gex_prob_load_gene_ids: argument is null");

    str_map *map = read_str_map(fn);
    if (map == NULL)
        return -1;

    if (gp->gene_ix != NULL){
        destroy_str_map(gp->gene_ix);
        gp->gene_ix = NULL;
    }

    gp->gene_ix = map;

    return 0;
}

int gex_prob_set_prob(gex_prob_t *gp, double *rho, 
        int rho_ncol, int rho_nrow, 
        double *spl, int spl_len){
    if (gp == NULL || rho == NULL || spl ==  NULL)
        return -1;

    assert(rho_ncol >= 0);
    assert(rho_nrow >= 0);
    assert(spl_len >= 0);

    gp->k += rho_ncol;

    uint16_t i;
    uint16_t init_len = gp->probs_len;

    // realloc gene_probs and spl_probs
    gp->probs_len += rho_ncol;
    gp->gene_probs = realloc(gp->gene_probs, gp->probs_len * sizeof(cat_ds_t *));
    if (gp->gene_probs == NULL)
        return err_msg(-1, 0, "gex_prob_set_prob: %s", strerror(errno));
    gp->spl_probs = realloc(gp->spl_probs, gp->probs_len * sizeof(cat_ds_t *));
    if (gp->spl_probs == NULL)
        return err_msg(-1, 0, "gex_prob_set_prob: %s", strerror(errno));
    for (i = init_len; i < gp->probs_len; ++i) gp->gene_probs[i] = NULL;
    for (i = init_len; i < gp->probs_len; ++i) gp->spl_probs[i] = NULL;

    // add probs
    double spl_probs[2];
    for (i = 0; i < rho_ncol; ++i){
        // gene probs index
        int p_ix = i + init_len;

        // add gene probs
        gp->gene_probs[p_ix] = cat_ds_alloc();
        if (cat_ds_set_p(gp->gene_probs[p_ix], rho + (i * rho_nrow), rho_nrow) < 0)
            return -1;

        // add spl
        int spl_ix = i % spl_len; // recycle spl values.
        if (spl[spl_ix] < 0 || spl[spl_ix] > 1)
            return err_msg(-1, 0, "gex_prob_set_prob: splice prob. %f "
                    "must be within [0,1]", spl[spl_ix]);
        spl_probs[0] = 1.0 - spl[spl_ix]; // index 0 is pre-mRNA
        spl_probs[1] = spl[spl_ix]; // index 1 is mature mRNA
        gp->spl_probs[p_ix] = cat_ds_alloc();
        if (cat_ds_set_p(gp->spl_probs[p_ix], spl_probs, 2) < 0)
            return -1;
    }

    return rho_ncol;
}

int gex_prob_load_prob(gex_prob_t *gp, const char *prob_fn, 
        const char *spl_fn){
    if (gp == NULL || prob_fn == NULL || spl_fn == NULL)
        return err_msg(-1, 0, "gex_prob_load_prob: argument is null");

    int rret = 0;
    double **arr = NULL;
    char **rownames = NULL, **colnames = NULL;
    char *delim = "\t", newline = '\n';

    // get the gene expression probabilities
    int nrow = 0, ncol = 0, rowcol = 0, header = 0;
    rret = read_matrix_double(prob_fn, &arr, rowcol, header, 
            &rownames, &nrow, &colnames, &ncol, 
            delim, newline);
    if (rret < 0)
        return err_msg(-1, 0, "gex_prob_load_prob: failed to read matrix from file");

    double *rho = malloc(nrow * ncol * sizeof(double));
    if (rho == NULL)
        return err_msg(-1, 0, "gex_prob_load_prob: %s", strerror(errno));

    int i, j;
    for (i = 0; i < nrow; ++i){
        for (j = 0; j < ncol; ++j)
            rho[CMI(i,j,nrow)] = arr[i][j];
    }

    int rho_ncol = ncol, rho_nrow = nrow;

    for (i = 0; i < nrow; ++i){
        free(arr[i]);
        if (rowcol) free(rownames[i]);
    }
    for (j = 0; j < ncol; ++j){
        if (header) free(colnames[j]);
    }
    free(arr); free(rownames); free(colnames);
    arr = NULL; rownames = NULL; colnames = NULL;

    // get the splice probabilities
    nrow = 0; ncol = 0; rowcol = 0; header = 0;
    rret = read_matrix_double(spl_fn, &arr, rowcol, header, 
            &rownames, &nrow, &colnames, &ncol, 
            delim, newline);
    if (rret < 0)
        return err_msg(-1, 0, "gex_prob_load_prob: failed to read matrix from file");
    if (ncol <= 0)
        return err_msg(-1, 0, "gex_prob_load_prob: no data in spl file %s", 
                spl_fn);

    double *spl = malloc(ncol * sizeof(double));
    if (spl == NULL)
        return err_msg(-1, 0, "gex_prob_load_prob: %s", strerror(errno));

    i = 0;
    for (j = 0; j < ncol; ++j)
        spl[CMI(i,j,nrow)] = arr[i][j];

    for (i = 0; i < nrow; ++i){
        free(arr[i]);
        if (rowcol) free(rownames[i]);
    }
    for (j = 0; j < ncol; ++j){
        if (header) free(colnames[j]);
    }
    free(arr); free(rownames); free(colnames);
    arr = NULL; rownames = NULL; colnames = NULL;

    int spl_len = ncol;

    // add the probs to gex_prob
    rret = gex_prob_set_prob(gp, rho, rho_ncol, rho_nrow, spl, spl_len);
    if (rret < 0) return -1;

    free(rho);
    free(spl);

    return rho_ncol;
}

int gex_prob_check_num_genes(gex_prob_t *gp) {
    if (gp == NULL)
        return err_msg(-1, 0, "gex_prob_check_num_genes: argument is null");

    if (gp->gene_probs == NULL)
        return err_msg(-1, 0, "gex_prob_check_num_genes: gene probs hasn't been set");
    if (gp->gene_ix == NULL)
        return err_msg(-1, 0, "gex_prob_check_num_genes: gene ixs hasn't been set");
    uint64_t gene_ix_un = (uint64_t)gp->gene_ix->n;
    if (gene_ix_un != gp->gene_probs[0]->n)
        return err_msg(-1, 0, "gex_prob_check_num_genes: gene numbers from 'gene-ids' "
                "and 'expr-prob' don't match");

    return 0;
}

int gex_prob_sample_gene(gex_prob_t *gp, uint16_t k, gene_t **gene){
    if (gp == NULL || gene == NULL) return -1;

    // check that k in [0, probs_len)
    if (k >= gp->probs_len)
        return err_msg(-1, 0, "gex_prob_sample_gene: k=%u must be < %u", 
                k, gp->probs_len);

    // get random gene index
    int rv_ret = 0;
    uint64_t r_ix = cat_ds_rv(gp->gene_probs[k], &rv_ret);
    if (rv_ret < 0)
        return -1;

    char *gene_id = str_map_str(gp->gene_ix, r_ix);
    assert(gene_id);

    int g_ix = str_map_ix(gp->anno->gene_ix, gene_id);
    if (g_ix < 0)
        return err_msg(-1, 0, "gex_prob_sample_gene: gene '%s' given in prob "
                "file but not present in GTF", gene_id);

    // get gene from index
    *gene = mv_i(&gp->anno->gix2gene, g_ix);
    if (*gene == NULL)
        return err_msg(-1, 0, "gex_prob_sample_gene: *gene is null");

    return g_ix;
}

int gex_prob_sample_tx(gex_prob_t *gp, gene_t *gene, isoform_t **iso){
    if (gp == NULL || gene == NULL || iso == NULL)
        return err_msg(-1, 0, "gex_prob_sample_tx: argument is null");

    // make sure n_tx > 0
    int n_tx = gene->isoforms_n;
    assert(n_tx > 0);

    int tx_ix = cat_ds_uni_rand(0, n_tx);

    kbtree_t(kb_iso) *bt = gene->bt_isoforms;
    assert(bt);
    int kb_ix = 0;
    kbitr_t itr;
    kb_itr_first(kb_iso, bt, &itr);
    for (; kb_itr_valid(&itr); kb_itr_next(kb_iso, bt, &itr)){
        *iso = &kb_itr_key(isoform_t, &itr);
        if (kb_ix == tx_ix) break;
        ++kb_ix;
    }
    assert(kb_ix >= 0);
    assert(kb_ix < n_tx);
    assert(kb_itr_valid(&itr));
    return tx_ix;
}

int gex_prob_sample_spl(gex_prob_t *gp, uint16_t k){
    if (gp == NULL)
        return err_msg(-1, 0, "gex_prob_sample_spl: argument is null");

    if (k >= gp->probs_len)
        return err_msg(-1, 0, "gex_prob_sample_spl: k=%u must be < %u", 
                k, gp->probs_len);

    // get random gene index
    int rv_ret = 0;
    uint64_t spl = cat_ds_rv(gp->spl_probs[k], &rv_ret);
    if (rv_ret < 0)
        return -1;

    return spl;
}

int isoform_mat_mrna_range(isoform_t *iso, int_ranges_t *ranges){
    // return the genomic position ranges of the coding sequence
    if (iso == NULL || ranges == NULL)
        return err_msg(-1, 0, "isoform_mrna_range: argument is null");
    
    ml_node_t(exon_list) *e_node;
    for (e_node = ml_begin(&iso->exons); 
            e_node; 
            e_node = ml_node_next(e_node)){
        exon_t exon = ml_node_val(e_node);
        int_range_t range;
        range.beg = exon.beg;
        range.end = exon.end;
        if (int_ranges_add_range(ranges, range) < 0)
            return -1;
    }
    return 0;
}

int isoform_pre_mrna_range(isoform_t *iso, int_ranges_t *ranges){
    // return the genomic position ranges of the coding sequence
    if (iso == NULL || ranges == NULL)
        return err_msg(-1, 0, "isoform_mrna_range: argument is null");
    
    int_range_t range;
    range.beg = iso->beg;
    range.end = iso->end;
    if (int_ranges_add_range(ranges, range) < 0)
        return -1;

    return 0;
}

int gex_sample_read(sc_sim_t *sc_sim, uint16_t k, int rsam, rna_read_t *rna_read){
    if (sc_sim == NULL || rna_read == NULL)
        return err_msg(-1, 0, "gex_sample_read: argument is null");

    gex_prob_t *gp = sc_sim->gex_prob;
    g_var_t *gv = sc_sim->gv;
    fa_seq_t *fa = sc_sim->fa;
    assert(gp && gv && fa);
    uint32_t read_len = sc_sim->rna_rd_len;

    // sample a gene and make sure read can fit
    gene_t *gene = NULL;
    int gene_ix;
    uint32_t gene_cds_len = 0;
    do {
        if ((gene_ix = gex_prob_sample_gene(gp, k, &gene)) < 0)
            return -1;
        gene_cds_len = gene->min_cds_len;
        assert(gene_cds_len > 0);
    } while (gene_cds_len < read_len);

    // sample an isoform and make sure read can fit
    isoform_t *iso = NULL;
    int iso_ix;
    uint32_t iso_cds_len = 0;
    do {
        if ((iso_ix = gex_prob_sample_tx(gp, gene, &iso)) < 0)
            return -1;
        iso_cds_len = iso->cds_len;
        assert(iso_cds_len > 0);
    } while (iso_cds_len < read_len);

    // sample whether mature RNA or premature (unspliced RNA)
    int mat_rna = gex_prob_sample_spl(gp, k);
    if (mat_rna < 0)
        return -1;

    // get the position range(s) of the rna
    int_ranges_t ranges;
    int_ranges_init(&ranges);
    // if mature
    if (mat_rna > 0){
        if (isoform_mat_mrna_range(iso, &ranges) < 0)
            return -1;
    } else {
        if (isoform_pre_mrna_range(iso, &ranges) < 0)
            return -1;
    }

    // get the sequence from the range
    // first get the chromosome name of the gene from the GTF file
    int rid = gene->chrm;
    const char *c_name = str_map_str(gp->anno->chrm_ix, rid);
    assert(c_name);
    seq_ranges_t *seq_ranges = NULL;
    if (fa_seq_seq_ranges(fa, c_name, ranges, &seq_ranges) < 0)
        return -1;
    if (seq_ranges == NULL)
        return err_msg(-1, 0, "gex_sample_read: failed to get range from fasta");

    // sample a read of length read_len uniformly
    // get sample position
    int q_len = seq_ranges->len - (int)read_len + 1;
    assert(q_len > 0);
    int pos_sample = rand() % q_len;
    assert(pos_sample < q_len);

    // get subset of gene from sampled position
    seq_ranges_t *subset = malloc(sizeof(seq_ranges_t));
    if (subset == NULL)
        return err_msg(-1, 0, "gex_sample_read: %s", strerror(errno));
    seq_ranges_init(subset);

    int sret = seq_ranges_subset(seq_ranges, subset, pos_sample, read_len);
    if (sret < 0) return -1;

    // get overlapping variants of the subset
    if (seq_ranges_var(subset, c_name, gv) < 0)
        return -1;

    // sample allele
    int vret = seq_ranges_sample_allele(subset, gv->vcf_hdr, rsam);
    if (vret < 0)
        return -1;

    // set allele in sequence
    vret = seq_ranges_set_allele_seq(subset, gv->vcf_hdr);
    if (vret < 0)
        return -1;

    // add sequencing error
    size_t p_ix, p_num = mv_size(&subset->rv);
    for (p_ix = 0; p_ix < p_num; ++p_ix) {
        seq_range_t *sr = &mv_i(&subset->rv, p_ix);
        if (seq_range_seq_error(sr, sc_sim->seq_error) < 0)
            return -1;
    }

    rna_read->seq_ranges = subset;
    rna_read->gene = gene;
    rna_read->iso = iso;
    rna_read->mat_rna = mat_rna;

    int_ranges_free(&ranges);
    seq_ranges_free(seq_ranges);
    free(seq_ranges);

    return 0;
}

void rna_read_init(rna_read_t *read){
    if (read == NULL) return;
    // seq_ranges_init(read->seq_ranges);
    read->seq_ranges = NULL;
    int i;
    for (i = 0; i < 2; ++i) {
        read->name[i] = NULL;
        read->seq[i] = NULL;
        mv_init(&read->qual[i]);
    }
    read->gene = NULL;
    read->iso = NULL;
    read->mat_rna = 0;
}

void rna_read_free(rna_read_t *read) {
    if (read == NULL) return;

    seq_ranges_free(read->seq_ranges);
    free(read->seq_ranges);
    read->seq_ranges = NULL;
    int i;
    for (i = 0; i < 2; ++i) {
        free(read->name[i]);
        read->name[i] = NULL;
        free(read->seq[i]);
        read->seq[i] = NULL;
        mv_free(&read->qual[i]);
    }
}

int rna_read_set_name(rna_read_t *rna_read, il_qname_t *names){
    if (rna_read == NULL)
        return err_msg(-1, 0, "rna_read_set_name: argument is null");
    if (names == NULL)
        return err_msg(-1, 0, "rna_read_set_name: names is null");
    
    char i5[] = "GTAACATGCG";
    char i7[] = "AGGTAACACT";

    char *base_name = il_qname_get_name(names);
    if (base_name == NULL)
        return -1;

    int i;
    int s_len[2] = {0, 0};
    size_t name_size[2] = {1000, 1000};
    char *name_a[2] = {NULL, NULL};
    for (i = 0; i < 2; ++i) {
        name_a[i] = calloc(name_size[i], sizeof(char));
        if (name_a[i] == NULL)
            return err_msg(-1, 0, "rna_read_set_name: %s", strerror(errno));

        // set read name in name_a
        s_len[i] = snprintf(name_a[i], name_size[i], "%s %i:N:0:%s+%s",
                base_name, i+1, i5, i7);
        if (s_len[i] < 0)
            return err_msg(-1, 0, "rna_read_set_name: failed to convert string");
        assert((size_t)s_len[i] < name_size[i]);
        name_a[i] = realloc(name_a[i], (s_len[i] + 1) * sizeof(char));
        rna_read->name[i] = name_a[i];
    }
    free(base_name);

    return 0;
}

int rna_read_set_seq(rna_read_t *rna_read, const char *bc_name,
        size_t umi_len) {
    if (rna_read == NULL || bc_name == NULL)
        return err_msg(-1, 0, "rna_read_set_seq: argument is null");

    size_t i;
    
    // get a random UMI sequence
    char *umi_seq = seq_rand(umi_len);
    if (umi_seq == NULL)
        return err_msg(-1, 0, "rna_read_set_seq: %s", strerror(errno));

    // cat BC + UMI
    rna_read->seq[0] = strcat2(bc_name, umi_seq);
    free(umi_seq);
    if (rna_read->seq[0] == NULL)
        return err_msg(-1, 0, "rna_read_set_seq: failed to cat strings");

    // allocate the sequence, even if empty
    size_t seq2_len = rna_read->seq_ranges->len;
    rna_read->seq[1] = malloc((seq2_len + 1) * sizeof(char));
    if (rna_read->seq[1] == NULL){
        err_msg(-1, 0, "rna_read_set_seq: %s", strerror(errno));
        return -1;
    }

    size_t n_seq = mv_size(&rna_read->seq_ranges->rv);
    for (i = 0; i <= seq2_len; ++i) rna_read->seq[1][i] = '\0';
        
    size_t src_tot_len = 0;
    for (i = 0; i < n_seq; ++i){
        seq_range_t *sr = &mv_i(&rna_read->seq_ranges->rv, i);
        char *src_seq = sr->seq;
        size_t src_len = (size_t)(sr->range.end - sr->range.beg);
        strncat(rna_read->seq[1], src_seq, src_len);
        src_tot_len += src_len;
    }
    assert(src_tot_len == seq2_len);

    if (rna_read->gene->strand == '-')
        seq_rev_cpl(rna_read->seq[1]);

    // set base qualities
    for (i = 0; i < 2; ++i) {
        size_t j, seq_len = strlen(rna_read->seq[i]);
        if (mv_resize(ui8v, &rna_read->qual[i], seq_len) < 0)
            return err_msg(-1, 0, "rna_read_set_seq: failed to set qual");

        for (j = 0; j < seq_len; ++j)
            mv_i(&rna_read->qual[i], j) = 37;
        mv_size(&rna_read->qual[i]) = seq_len;
    }

    return 0;
}

int rna_read_seq_error(sc_sim_t *sc_sim, rna_read_t *rna_read) {
    if (sc_sim == NULL || rna_read == NULL)
        return err_msg(-1, 0, "rna_read_seq_error: argument is null");

    if (rna_read->seq[1] == NULL)
        return err_msg(-1, 0, "rna_read_seq_error: seq not set");

    int b, arr_len;
    int *arr = binom_ds_sample(&sc_sim->rna_err_prob, &arr_len);
    if (arr_len < 0)
        return -1;

    for (b = 0; b < arr_len; ++b)
        rna_read->seq[1][arr[b]] = base_rand();

    free(arr);

    return 0;
}

int rna_read_write(rna_read_t *rna_read, BGZF *fs[2]) {
    if (rna_read == NULL || fs == NULL)
        return err_msg(-1, 0, "rna_read_write_r1: argument is null");

    int wret, i;

    for (i = 0; i < 2; ++i) {
        if (fs[i] == NULL)
            return err_msg(-1, 0, "rna_read_write: fastq stream %i "
                    "is null", i);
    }

    for (i = 0; i < 2; ++i) {
        size_t j, n_bp = mv_size(rna_read->qual + i);
        char *qual_s = calloc(n_bp + 1, sizeof(char));
        if (qual_s == NULL)
            return err_msg(-1, 0, "rna_read_write: %s", strerror(errno));
        for (j = 0; j < n_bp; ++j)
            qual_s[j] = mv_i(rna_read->qual + i, j) + 33;
        qual_s[j] = '\0';

        int ix, sa_len = 6;
        char *sa[6];
        sa[0] = rna_read->name[i];
        sa[1] = "\n";
        sa[2] = rna_read->seq[i];
        sa[3] = "\n+\n";
        sa[4] = qual_s;
        sa[5] = "\n";

        size_t str_len = 0;
        for (ix = 0; ix < sa_len; ++ix)
            str_len += strlen(sa[ix]);

        char *str = calloc(str_len + 1, sizeof(char));
        if (str == NULL)
            return err_msg(-1, 0, "rna_read_write: %s", strerror(errno));
        str[0] = '\0';
        for (ix = 0; ix < sa_len; ++ix)
            strcat(str, sa[ix]);

        wret = bgzf_write(fs[i], str, str_len);

        if (wret < 0)
            return err_msg(-1, 0, "rna_read_write: failed to write to fastq "
                    "stream");

        free(qual_s);
        free(str);
    }

    return 0;
}

/*
   int test_sample_gex(gex_prob_t *gp, fa_seq_t *fa, int K, int n_sample, 
   uint32_t read_len, g_var_t *gv){
   if (gp == NULL || fa == NULL) return -1;
   if (n_sample < 0) n_sample = 10;

   srand(time(NULL));
   int vcf_n_sam = bcf_hdr_nsamples(gv->vcf_hdr);
   int i;
   for (i = 0; i < n_sample; ++i){
   rna_read_t read;
   rna_read_init(&read);
   int k = rand() % K;
   int rsam = rand() % vcf_n_sam;
   int sret = gex_sample_read(gp, (uint16_t)k, rsam, fa, read_len, 
   &read, gv);
   if (sret < 0) return -1;
   const char bc[9] = "ATCGATCG";
   size_t umi_len = 8;
   if (rna_read_set_seq(&read, bc, umi_len) < 0) return -1;
   if (read.seq_ranges->n_var > 0){
   fprintf(stdout, "################### sample\n");
   fprintf(stdout, "%s\n", read.seq[1]);
   fprintf(stdout, "\nseq len=%i", read.seq_ranges->len);
   fprintf(stdout, "\tk=%i", k);
   fprintf(stdout, "\tsample=%s", gv->vcf_hdr->samples[rsam]);
   fprintf(stdout, "\tgene=%s", read.gene->id);
   fprintf(stdout, "\tisoform=%s", read.iso->id);
   fprintf(stdout, "\tspl=%i", read.mat_rna);
   fprintf(stdout, "\tn_var=%i\n", read.seq_ranges->n_var);
   }

   if (sret < 0) return -1;
   rna_read_free(&read);
   }

   return 0;
   }
   */

// given an isoform_t object, return a splice range vector
// given a splice range vector, subset given a position and a length
// given a splice range vector, return the genomic sequence

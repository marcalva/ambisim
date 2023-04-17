
#include "atac_prob.h"
#include "array_util.h"


void atac_prob_init(atac_prob_t *ap){
    if (ap == NULL) return;
    
    ap->k = 0;
    ap->peaks = iregs_init();
    ap->peak_probs = NULL;
    ap->probs_len = 0;
}

void atac_prob_free(atac_prob_t *ap){
    if (ap == NULL) return;
    
    iregs_dstry(ap->peaks);
    ap->peaks = NULL;
    uint16_t i;
    for (i = 0; i < ap->probs_len; ++i)
        cat_ds_dstry(ap->peak_probs[i]);
    free(ap->peak_probs);
    ap->peak_probs = NULL;
    ap->probs_len = 0;
}

atac_prob_t *atac_prob_alloc(){
    atac_prob_t *ap = malloc(sizeof(atac_prob_t));
    if (ap == NULL){
        err_msg(-1, 0, "atac_prob_alloc: %s", strerror(errno));
        return NULL;
    }
    atac_prob_init(ap);
    return ap;
}

void atac_prob_dstry(atac_prob_t *ap){
    if (ap == NULL) return;
    atac_prob_free(ap);
    free(ap);
}

int atac_prob_load_peaks(atac_prob_t *ap, const char *peaks_fn){
    if (ap == NULL || peaks_fn == NULL)
        return err_msg(-1, 0, "atac_load_peaks: argument is null");

    if (iregs_add_bed(ap->peaks, peaks_fn) < 0) return -1;
    if (iregs_parse_bed(ap->peaks) < 0) return -1;

    return 0;
}

int atac_prob_add_prob(atac_prob_t *ap, double *probs, int ncol, int nrow){
    if (ap == NULL || probs == NULL)
        return err_msg(-1, 0, "atac_prob_add_prob: argument is null");

    assert(ncol > 0);
    assert(nrow > 0);

    ap->k = ncol; // number of cell types + ambient pool

    uint16_t i;
    // allocate probs array
    ap->probs_len = ncol;
    ap->peak_probs = malloc(ap->probs_len * sizeof(cat_ds_t *));
    if (ap->peak_probs == NULL)
        return err_msg(-1, 0, "atac_prob_add_prob: %s", strerror(errno));
    for (i = 0; i < ap->probs_len; ++i) ap->peak_probs[i] = NULL;

    // add peak probs
    for (i = 0; i < ncol; ++i){
        ap->peak_probs[i] = cat_ds_alloc();
        if (cat_ds_set_p(ap->peak_probs[i], probs + (i * nrow), nrow) < 0)
            return -1;
    }

    return ncol;
}

int atac_prob_load_file(atac_prob_t *ap, const char *prob_fn){
    if (ap == NULL || prob_fn == NULL)
        return err_msg(-1, 0, "atac_prob_load_file: argument is null");

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
        return err_msg(-1, 0, "atac_prob_load_file: failed to read matrix from file");

    double *probs = malloc(nrow * ncol * sizeof(double));
    if (probs == NULL)
        return err_msg(-1, 0, "atac_prob_load_file: %s", strerror(errno));

    int i, j;
    for (i = 0; i < nrow; ++i){
        for (j = 0; j < ncol; ++j)
            probs[CMI(i,j,nrow)] = arr[i][j];
    }

    // free array
    for (i = 0; i < nrow; ++i){
        free(arr[i]);
    }
    free(arr);

    rret = atac_prob_add_prob(ap, probs, ncol, nrow);
    if (rret < 0) return -1;

    free(probs);

    return 0;
}

int atac_prob_sample_peak(atac_prob_t *ap, uint16_t k, g_region *reg){
    if (ap == NULL)
        return err_msg(-1, 0, "atac_prob_sample_peak: argument is null");

    if (k >= ap->probs_len)
        return err_msg(-1, 0, "atac_prob_sample_peak: k=%u must be < %u", 
                k, ap->probs_len);

    // get random peak index
    int rv_ret = 0;
    uint64_t p_ix = cat_ds_rv(ap->peak_probs[k], &rv_ret);
    if (rv_ret < 0)
        return -1;

    // get region from index
    *reg = ap->peaks->reg[p_ix];

    return 0;
}

int peak_range(g_region *reg, int_ranges_t *ranges){
    if (reg == NULL || ranges == NULL)
        return err_msg(-1, 0, "peak_range: argument is null");

    int_range_t range;
    range.beg = reg->start;
    range.end = reg->end;
    if (int_ranges_add_range(ranges, range) < 0)
        return -1;

    return 0;
}

int atac_prob_sample_read(sc_sim_t *sc_sim, uint16_t k, int peak, int rsam, atac_pair_t *pair){
    // check input
    if (sc_sim == NULL || pair == NULL)
        return err_msg(-1, 0, "atac_prob_sample_read: argument is null");

    atac_prob_t *ap = sc_sim->atac_prob;
    fa_seq_t *fa = sc_sim->fa;
    g_var_t *gv = sc_sim->gv;
    assert(ap && fa && gv);
    uint32_t read_len = sc_sim->atac_rd_len;

    if (k >= ap->probs_len)
        return err_msg(-1, 0, "atac_prob_sample_read: k=%u must be < %u", 
                k, ap->probs_len);

    // TODO: create a better random insert size sample
    int32_t ins_size = (rand() % (150 + read_len)) - read_len;
    uint32_t read_pair_len = 2 * read_len + ins_size;

    // get sampling region
    g_region pk_reg;
    char *c_name;
    if (peak){
        // if read is in peak, sample a peak
        uint32_t pk_len;
        do {
            if (atac_prob_sample_peak(ap, k, &pk_reg) < 0)
                return -1;
            pk_len = pk_reg.end - pk_reg.start;
            c_name = str_map_str(ap->peaks->chr_map, pk_reg.rid);
        } while (pk_len < read_pair_len);
    } else {
        // otherwise, sample a region from the genome. No 'N' bp allowed
        // TODO sample peak size
        int n_n;
        do {
            if (fa_seq_rand_range(fa, read_pair_len + 600, '.', &pk_reg) < 0)
                return -1;
            c_name = str_map_str(fa->c_names, pk_reg.rid);
            n_n = fa_seq_n_n(fa, c_name, pk_reg.start, pk_reg.end);
        } while (n_n > 0);
        
    }

    // sample a read pair, each of length read_len, separated by 
    // sample the read inside the peak region
    int pk_len = pk_reg.end - pk_reg.start;
    int q_len = pk_len - (int)read_pair_len + 1;
    assert(q_len > 0);
    int pos_sample = rand() % q_len;
    assert(pos_sample < q_len);
    int beg1 = pk_reg.start + pos_sample, end1 = beg1 + read_len;
    int beg2 = end1 + ins_size, end2 = beg2 + read_len;

    // get int range subset of the sampled read
    int_range_t rd_rng[2];
    rd_rng[0].beg = beg1;
    rd_rng[0].end = end1;
    rd_rng[1].beg = beg2;
    rd_rng[1].end = end2;

    int_range_t rd1_rng, rd2_rng;
    int_range_init(&rd1_rng);
    int_range_init(&rd2_rng);
    rd1_rng.beg = beg1;
    rd1_rng.end = end1;
    rd2_rng.beg = beg2;
    rd2_rng.end = end2;
    assert( (rd2_rng.end - rd1_rng.beg) == (int)read_pair_len );

    int i, n_read = 2;
    for (i = 0; i < n_read; ++i) {
        int_ranges_t read_rng;
        int_ranges_init(&read_rng);
        if (int_ranges_add_range(&read_rng, rd_rng[i]) < 0)
            return -1;
        seq_ranges_t *read_seq = NULL;
        // printf("seq for %s:%i-%i\n", c_name, rd_rng[i].beg, rd_rng[i].end);
        if (fa_seq_seq_ranges(fa, c_name, read_rng, &read_seq) < 0)
            return -1;
        if (read_seq == NULL)
            return err_msg(-1, 0, "atac_prob_sample_read: failed to get range from fasta");

        // get overlapping variants of the subset
        if (seq_ranges_var(read_seq, c_name, gv) < 0)
            return -1;

        // sample allele
        int vret = seq_ranges_sample_allele(read_seq, gv->vcf_hdr, rsam);
        if (vret < 0)
            return -1;

        // set allele in sequence
        vret = seq_ranges_set_allele_seq(read_seq, gv->vcf_hdr);
        if (vret < 0)
            return -1;

        // add sequencing error
        size_t p_ix, p_num = mv_size(&read_seq->rv);
        for (p_ix = 0; p_ix < p_num; ++p_ix) {
            seq_range_t *sr = &mv_i(&read_seq->rv, p_ix);
            if (seq_range_seq_error(sr, sc_sim->seq_error) < 0)
                return -1;
        }

        pair->pair[i].seq_ranges = read_seq;
        pair->pair[i].peak = peak ? 1 : 0;
        uint32_t bp_i, rlen = read_seq->len;
        if (mv_resize(ui8v, &pair->pair[i].qual, rlen) < 0)
            return -1;
        for (bp_i = 0; bp_i < rlen; ++bp_i)
            mv_push(ui8v, &pair->pair[i].qual, 37);

        int_ranges_free(&read_rng);
    }

    // reverse complement one of the reads
    int pair_ix = rand() % 2;
    seq_range_t *seq_r = &mv_i(&pair->pair[pair_ix].seq_ranges->rv, 0);
    if (seq_rev_cpl(seq_r->seq) < 0)
        return -1;
    seq_r->rc = 1;

    return 0;
}

/*
int test_sample_atac(atac_prob_t *ap, g_var_t *gv, fa_seq_t *fa, 
        int n_reads, uint32_t read_len){
    if (ap == NULL) return -1;

    srand(time(NULL));
    int vcf_n_sam = bcf_hdr_nsamples(gv->vcf_hdr);
    int i;
    for (i = 0; i < n_reads; ++i){
        atac_pair_t atac_pair;
        atac_pair_init(&atac_pair);
        int k = rand() % ap->k, pk = rand() % 2;
        int rsam = rand() % vcf_n_sam;
        g_region reg;
        int pret = atac_prob_sample_peak(ap, k, &reg);
        if (pret < 0) return -1;
        fprintf(stdout, "k=%i %i:%i-%i len=%i\n", k, reg.rid, reg.start, reg.end, 
                reg.end - reg.start);
        if (atac_prob_sample_read(ap, k, pk, rsam, fa, read_len, gv, &atac_pair) < 0)
            return -1;
        char *seq1 = mv_i(&atac_pair.pair[0].seq_ranges->rv, 0).seq;
        char *seq2 = mv_i(&atac_pair.pair[1].seq_ranges->rv, 0).seq;
        fprintf(stdout, "################### sample\n");
        fprintf(stdout, "%s\n", seq1);
        fprintf(stdout, "%s\n", seq2);
        fprintf(stdout, "pk=%i\n", pk);
        atac_pair_free(&atac_pair);
    }

    return 0;
}
*/

/*******************************************************************************
 * atac_read_t
 ******************************************************************************/

void atac_read_init(atac_read_t *atac_read){
    if (atac_read == NULL) return;
    atac_read->seq_ranges = NULL;
    mv_init(&atac_read->qual);
    atac_read->peak = 0;
}

void atac_read_free(atac_read_t *atac_read){
    if (atac_read == NULL) return;
    seq_ranges_free(atac_read->seq_ranges);
    free(atac_read->seq_ranges);
    atac_read->seq_ranges = NULL;
    mv_free(&atac_read->qual);
    atac_read->peak = 0;
}

void atac_pair_init(atac_pair_t *atac_pair) {
    if (atac_pair == NULL) return;
    atac_read_init(atac_pair->pair);
    atac_read_init(atac_pair->pair + 1);
    int i;
    for (i = 0; i < 3; ++i) {
        atac_pair->name[i] = NULL;
        atac_pair->seq[i] = NULL;
        mv_init(atac_pair->qual + i);
    }
}

void atac_pair_free(atac_pair_t *atac_pair) {
    if (atac_pair == NULL) return;
    atac_read_free(atac_pair->pair);
    atac_read_free(atac_pair->pair + 1);
    int i;
    for (i = 0; i < 3; ++i) {
        free(atac_pair->name[i]);
        atac_pair->name[i] = NULL;
        free(atac_pair->seq[i]);
        atac_pair->seq[i] = NULL;
        mv_free(atac_pair->qual + i);
    }
}

int atac_pair_set_name(atac_pair_t *atac_pair, il_qname_t *names) {
    if (atac_pair == NULL)
        return err_msg(-1, 0, "atac_pair_set_name: argument is null");
    if (names == NULL)
        return err_msg(-1, 0, "rna_read_set_name: names is null");

    char i7[] = "AAACGGCG";

    char *base_name = il_qname_get_name(names);
    if (base_name == NULL)
        return -1;

    int i;
    int s_len[3] = {0, 0, 0}; // length of string in name_a
    size_t name_size[3] = {1000, 1000, 1000}; // allocated size in name_a
    char *name_a[3] = {NULL, NULL, NULL};
    for (i = 0; i < 3; ++i) {
        // allocate array
        name_a[i] = calloc(name_size[i], sizeof(char));
        if (name_a[i] == NULL)
            return err_msg(-1, 0, "atac_pair_set_name: %s", strerror(errno));

        // set read name in name_a
        s_len[i] = snprintf(name_a[i], name_size[i], "%s %i:N:0:%s",
                base_name, i+1, i7);
        if (s_len[i] < 0)
            return err_msg(-1, 0, "atac_pair_set_name: failed to convert string");
        assert((size_t)s_len[i] < name_size[i]);
        name_a[i] = realloc(name_a[i], (s_len[i] + 1) * sizeof(char));
        atac_pair->name[i] = name_a[i];
    }
    free(base_name);

    return 0;
}

int atac_pair_set_seq(atac_pair_t *atac_pair, const char *bc_name) {
    if (atac_pair == NULL || bc_name == NULL)
        return err_msg(-1, 0, "atac_pair_set_seq: argument is null");

    size_t seq_len, i;
    char *seq_c[2] = {NULL, NULL};
    mv_t(ui8v) *qual_c[2] = {&atac_pair->qual[0], &atac_pair->qual[2]};

    // read 1 and read 3
    for (i = 0; i < 2; ++i) {
        atac_read_t *atac_read = atac_pair->pair + i;
        seq_len = atac_read->seq_ranges->len;
        seq_c[i] = malloc((seq_len + 1) * sizeof(char));
        if (seq_c[i] == NULL)
            return err_msg(-1, 0, "atac_pair_set_seq: %s", strerror(errno));
        size_t j;
        for (j = 0; j <= seq_len; ++j) seq_c[i][j] = '\0';
        size_t src_tot_len = 0, n_seq = mv_size(&atac_read->seq_ranges->rv);
        for (j = 0; j < n_seq; ++j) {
            seq_range_t *sr = &mv_i(&atac_read->seq_ranges->rv, j);
            assert(sr);
            assert(sr->seq);
            char *src_seq = sr->seq;
            size_t src_len = (size_t)(sr->range.end - sr->range.beg);
            strncat(seq_c[i], src_seq, src_len);
            src_tot_len += src_len;
        }
        assert(src_tot_len == seq_len);

        // read is already rev complemented so don't need to flip strand

        // set qual
        if (mv_resize(ui8v, qual_c[i], seq_len) < 0)
            return err_msg(-1, 0, "atac_pair_set_seq: failed to set qual");
        for (j = 0; j < seq_len; ++j)
            mv_i(qual_c[i], j) = 37;
        mv_size(qual_c[i]) = seq_len;
    }
    atac_pair->seq[0] = seq_c[0];
    atac_pair->seq[2] = seq_c[1];


    // read 2 (barcode)
    i = 1;
    atac_pair->seq[i] = strdup(bc_name);
    if (atac_pair->seq[i] == NULL)
        return err_msg(-1, 0, "atac_pair_set_seq: %s", strerror(errno));
    // barcode is rev complement from sequencing
    seq_rev_cpl(atac_pair->seq[i]);

    size_t bc_len = strlen(atac_pair->seq[i]);
    if (mv_resize(ui8v, &atac_pair->qual[1], seq_len) < 0)
        return err_msg(-1, 0, "atac_pair_set_seq: failed to set qual");

    for (i = 0; i < bc_len; ++i)
        mv_i(&atac_pair->qual[1], i) = 37;
    mv_size(&atac_pair->qual[1]) = bc_len;

    return 0;
}

int atac_pair_seq_error(sc_sim_t *sc_sim, atac_pair_t *atac_pair) {
    if (sc_sim == NULL || atac_pair == NULL)
        return err_msg(-1, 0, "atac_pair_seq_error: argument is null");

    if (atac_pair->seq[0] == NULL || atac_pair->seq[2] == NULL)
        return err_msg(-1, 0, "atac_pair_seq_error: seq not set");

    int p, rd_ix[2] = {0,2};

    for (p = 0; p < 2; ++p) {
        int p_ix = rd_ix[p];

        int b, arr_len;
        int *arr = binom_ds_sample(&sc_sim->atac_err_prob, &arr_len);
        if (arr_len < 0)
            return -1;

        for (b = 0; b < arr_len; ++b)
            atac_pair->seq[p_ix][arr[b]] = base_rand();
        free(arr);
    }

    return 0;
}

int atac_pair_write(atac_pair_t *atac_pair, BGZF *fs[3]) {
    if (atac_pair == NULL || fs == NULL)
        return err_msg(-1, 0, "atac_pair_write: argument is null");

    int wret, i;
    
    for (i = 0; i < 3; ++i) {
        if (fs[i] == NULL)
            return err_msg(-1, 0, "atac_pair_write: fastq stream %i "
                    "is null", i);
    }

    for (i = 0; i < 3; ++i) {
        size_t j, n_bp = mv_size(atac_pair->qual + i);
        char *qual_s = calloc(n_bp + 1, sizeof(char));
        if (qual_s == NULL)
            return err_msg(-1, 0, "atac_pair_write: %s", strerror(errno));
        for (j = 0; j < n_bp; ++j)
            qual_s[j] = mv_i(atac_pair->qual + i, j) + 33;
        qual_s[j] = '\0';

        int ix, sa_len = 6;
        char *sa[6];
        sa[0] = atac_pair->name[i];
        sa[1] = "\n";
        sa[2] = atac_pair->seq[i];
        sa[3] = "\n+\n";
        sa[4] = qual_s;
        sa[5] = "\n";

        size_t str_len = 0;
        for (ix = 0; ix < sa_len; ++ix)
            str_len += strlen(sa[ix]);

        char *str = calloc(str_len + 1, sizeof(char));
        if (str == NULL)
            return err_msg(-1, 0, "atac_pair_write: %s", strerror(errno));
        str[0] = '\0';
        for (ix = 0; ix < sa_len; ++ix)
            strcat(str, sa[ix]);

        wret = bgzf_write(fs[i], str, str_len);

        if (wret < 0)
            return err_msg(-1, 0, "atac_pair_write: failed to write to fastq "
                    "stream");
        
        free(qual_s);
        free(str);
    }

    return 0;
}


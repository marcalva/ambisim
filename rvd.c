
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <inttypes.h>
#include <errno.h>
#include <assert.h>
#include <float.h>
#include "kavl.h"
#include "htslib/faidx.h"
#include "htslib/vcf.h"
#include "g_list.h"
#include "region.h"
#include "rvd.h"
#include "str_util.h"
#include "bits.h"

uint64_t get_rand(uint64_t max){
    uint64_t r = (uint64_t)rand();
    r = r % max;
    return(r);
}

int seq_rev_cpl(char *seq){
    if (seq == NULL)
        return err_msg(-1, 0, "seq_rev_cpl: argument is null");

    int rc_err = 0;
    div_t r;
    long int seq_len = (long int)strlen(seq);
    r = div(seq_len, 2);
    long int i, mid = r.quot, extra = r.rem;
    for (i = 0; i < mid; ++i){
        long int j = seq_len - 1 - i;
        char tmp = seq[j];
        seq[j] = base_rev_cmp(seq[i], &rc_err);
        if (rc_err)
            return err_msg(-1, 0, "seq_rev_cpl: base '%c' "
                    "could not be reverse complemented", seq[j]);
        seq[i] = base_rev_cmp(tmp, &rc_err);
        if (rc_err)
            return err_msg(-1, 0, "seq_rev_cpl: base '%c' "
                    "could not be reverse complemented", seq[i]);
    }
    if (extra) seq[i] = base_rev_cmp(seq[i], &rc_err);
    if (rc_err)
        return err_msg(-1, 0, "seq_rev_cpl: base '%c' could not be "
                "reverse complemented", seq[i]);
    return 0;
}

char *seq_rand(size_t len) {
    char bases[4] = "ATCG";
    char *seq = pcalloc(len+1, sizeof(char), "seq_rand");
    if (seq == NULL)
        return NULL;
    size_t i;
    for (i = 0; i < len; ++i)
        seq[i] = bases[rand() % 4];
    seq[len] = '\0';
    return seq;
}

/*******************************************************************************
 * read name
 ******************************************************************************/

void il_qname_set_step(il_qname_t *il_qname) {
    if (il_qname == NULL) return;
    uint64_t rds_per_tile = il_qname->n_reads / il_qname->n_tiles + 1;
    assert(il_qname->spots_per_tile >= rds_per_tile);
    float fstp = (float)il_qname->spots_per_tile;
    float frpt = (float)rds_per_tile;
    il_qname->step = (uint64_t)sqrt(floorf(fstp / frpt));
}

void il_qname_init(il_qname_t *il_qname) {
    if (il_qname == NULL) return;

    il_qname->instr = NULL;
    il_qname->rn = -1;
    il_qname->fc = NULL;

    il_qname->min_lane = 1;
    il_qname->max_lane = 1;
    il_qname->min_t1 = 1;
    il_qname->max_t1 = 2;
    il_qname->min_t2 = 1;
    il_qname->max_t2 = 6;
    il_qname->min_t3 = 1;
    il_qname->max_t3 = 78;
    il_qname->min_spot = 1000;
    il_qname->max_spot = 20000;

    il_qname->lane = il_qname->min_lane;
    il_qname->tile1 = il_qname->min_t1;
    il_qname->tile2 = il_qname->min_t2;
    il_qname->tile3 = il_qname->min_t3;
    il_qname->x = il_qname->min_spot;
    il_qname->y = il_qname->min_spot;

    il_qname->n_tiles = (il_qname->max_lane - il_qname->min_lane + 1) *
        (il_qname->max_t1 - il_qname->min_t1 + 1) *
        (il_qname->max_t2 - il_qname->min_t2 + 1) *
        (il_qname->max_t3 - il_qname->min_t3 + 1);

    il_qname->spots_per_tile = (il_qname->max_spot - il_qname->min_spot + 1) *
        (il_qname->max_spot - il_qname->min_spot + 1);

    il_qname->n_reads = 1e9;
    il_qname_set_step(il_qname);
}

void il_qname_free(il_qname_t *il_qname) {
    if (il_qname == NULL) return;
    free(il_qname->instr);
    free(il_qname->fc);
}

char *il_qname_get_name(il_qname_t *il_qname) {
    if (il_qname == NULL){
        err_msg(-1, 0, "il_qname_get_name: argument is null");
        return NULL;
    }

    char *base_name = pmalloc(1000 * sizeof(char), "il_qname_get_name");
    if (base_name == NULL)
        return NULL;
    int base_len = snprintf(base_name, 1000, "@%s:%i:%s:%i:%i%i%02i:%i:%i",
                il_qname->instr, il_qname->rn, il_qname->fc,
                il_qname->lane,
                il_qname->tile1, il_qname->tile2, il_qname->tile3,
                il_qname->x, il_qname->y);
    if (base_len < 0){
        err_msg(-1, 0, "il_qname_get_name: failed to convert to string");
        return NULL;
    }
    assert(base_len < 1000);
    base_name = prealloc(base_name, (base_len + 1) * sizeof(char), "il_qname_get_name");
    if (base_name == NULL)
        return NULL;

    // increment counter
    il_qname->y += il_qname->step;
    if (il_qname->y > il_qname->max_spot) {
        il_qname->y = il_qname->min_spot;
        il_qname->x += il_qname->step;
    }
    if (il_qname->x > il_qname->max_spot) {
        il_qname->x = il_qname->min_spot;
        il_qname->tile3 += 1;
    }
    if (il_qname->tile3 > il_qname->max_t3) {
        il_qname->tile3 = il_qname->min_t3;
        il_qname->tile2 += 1;
    }
    if (il_qname->tile2 > il_qname->max_t2) {
        il_qname->tile2 = il_qname->min_t2;
        il_qname->tile1 += 1;
    }
    if (il_qname->tile1 > il_qname->max_t1) {
        il_qname->tile1 = il_qname->min_t1;
        il_qname->lane += 1;
    }
    if (il_qname->lane > il_qname->max_lane){
        err_msg(-1, 0, "il_qname_get_name: overrran lanes, likely a bug");
        return NULL;
    }
    base_name = prealloc(base_name, (base_len + 1) * sizeof(char), "il_qname_get_name");
    if (base_name == NULL)
        return NULL;

    return base_name;
}

void il_qname_set_instr(il_qname_t *il_qname, const char *instr,
        int rn, const char *fc) {
    if (il_qname == NULL)
        return;

    il_qname->instr = strdup(instr);
    il_qname->fc = strdup(fc);
    il_qname->rn = rn;
}

/*******************************************************************************
 * math
 ******************************************************************************/

int factorial(int x) {
    if (x < 0)
        return err_msg(-1, 0, "factorial: x=%i < 0", x);
    if (x > 12)
        return err_msg(-1, 0, "factorial: x=%i too large (overflow)", x);
    int y = 1;
    for (; x > 0; --x) y *= x;
    return y;
}

double binom_choose(double n, double k) {
    if (n < 0 || k < 0)
        return err_msg(-1, 0, "binom_choose: arguments are negative");

    if (k <= DBL_EPSILON)
        return 1;

    if (k > (n/2)) k = n - k;
    double x = 1;
    int i;
    for (i = 0; i < k; ++i) {
        double m = (n-i) / (k - i);
        if (x > (DBL_MAX / m))
            return err_msg(-1, 0, "binom_choose: choose(%i,%i) too large "
                    "(overflow)", (int)n, (int)k);
        x *= m;
    }
    return x;
}

/*******************************************************************************
 * Categorical distribution
 ******************************************************************************/

void cat_ds_init(cat_ds_t *d) {
    if (d == NULL) return;
    d->root = NULL;
    d->n = 0;
}

cat_ds_t *cat_ds_alloc(){
    cat_ds_t *d = pmalloc(sizeof(cat_ds_t), "cat_ds_alloc");
    if (d == NULL)
        return NULL;
    cat_ds_init(d);
    return(d);
}

void cat_ds_free(cat_ds_t *d){
    if (d == NULL) return;
    if (d->root) {
        kavl_itr_t(k_av_ix) itr;
        kavl_itr_first(k_av_ix, d->root, &itr);  // place at first
        do {                             // traverse
            const struct av_ix_t *p = kavl_at(&itr);
            free((void*)p);                // free node
        } while (kavl_itr_next(k_av_ix, &itr));
    }
    d->root = NULL;
    d->n = 0;
}

void cat_ds_dstry(cat_ds_t *d){
    if (d == NULL) return;
    cat_ds_free(d);
    free(d);
}

int cat_ds_set_p(cat_ds_t *d, double *p, uint64_t n){
    if (d == NULL){
        fprintf(stderr, "cat_ds_set_p: d is NULL\n");
        return -1;
    }

    if (n == 0)
        return 0;

    double p_total = 0.0;
    uint64_t i;
    for (i = 0; i < n; ++i) {
        if (p[i] < 0)
            return err_msg(-1, 0, "cat_ds_set_p: p[%i]=%f must be non-negative",
                    i, p[i]);
        p_total += p[i];
    }
    if (p_total <= 0.0)
        return err_msg(-1, 0, "cat_ds_set_p: all values in p are 0");

    // add to avl tree
    double p_this = 0.0;
    for (i = 0; i < n; ++i){
        // 0 probability values are skipped
        double prob = p[i] / p_total;
        if (prob <= DBL_EPSILON)
            continue;
        av_ix_t *s, *t = pmalloc(sizeof(av_ix_t), "cat_ds_set_p");
        if (t == NULL)
            return -1;
        p_this += prob;
        t->ix = i;
        t->val = p_this;
        t->prob = prob;
        // return t if not present, otherwise return duplicate node
        s = kavl_insert(k_av_ix, &d->root, t, NULL);
        // sanity check
        // if already present
        if (t != s) {
            fprintf(stderr, "prob=%f\n", prob);
            fprintf(stderr, "t->val = %f\ns->val = %f\n", t->val, s->val);
            free(t);
            return err_msg(-1, 0, "cat_ds_set_p: trying to set same prob value twice");
        }
    }

    d->n = n;

    return 0;
}

uint64_t cat_ds_rv(cat_ds_t *d, int *ret){
    *ret = 0;
    if (d == NULL) {
        *ret = 1;
        return err_msg(0, 0, "cat_ds_rv: argument is null");
        return 0;
    }

    if (d->n == 0 || d->root == NULL) {
        *ret = -1;
        return err_msg(0, 0, "cat_ds_rv: no elements in distribution");
    }
    if (kavl_size(head, d->root) == 0) {
        *ret = -1;
        return err_msg(0, 0, "cat_ds_rv: no elements in distribution");
    }

    int rmax = RAND_MAX;
    uint64_t r = get_rand(rmax);
    double rd = (double)r / (double)rmax;

    kavl_itr_t(k_av_ix) itr;
    av_ix_t pf;
    pf.val = rd;
    // place itr so it is greater than or equal to value of pf
    uint64_t ix;
    kavl_itr_find(k_av_ix, d->root, &pf, &itr);
    const av_ix_t *qf = kavl_at(&itr);
    // return last index if nothing greater. Shouldn't happen, maybe in edge 
    // case with rounding error
    if (qf == NULL) {
        kavl_itr_first(k_av_ix, d->root, &itr);
        do {
            const av_ix_t *qtmp = kavl_at(&itr);
            ix = qtmp->ix;
        } while (kavl_itr_next(k_av_ix, &itr));
    } else {
        ix = qf->ix;
    }
    return ix;
}

int cat_ds_uni_rand(int min, int max){
    int len = max - min;
    int r = rand() % len;
    r += min;
    return r;
}

/*******************************************************************************
 * Binomial distribution
 ******************************************************************************/

void binom_ds_init(binom_ds_t *x) {
    if (x == NULL) return;
    x->prob = -1;
    x->size = -1;
    cat_ds_init(&x->catd);
}

binom_ds_t *binom_ds_alloc() {
    binom_ds_t *x = pmalloc(sizeof(binom_ds_t), "binom_ds_alloc");
    binom_ds_init(x);
    return x;
}

void binom_ds_free(binom_ds_t *x) {
    if (x == NULL) return;
    cat_ds_free(&x->catd);
}

void binom_ds_dstry(binom_ds_t *x) {
    if (x == NULL) return;
    binom_ds_free(x);
    free(x);
}

int binom_ds_set(binom_ds_t *x, double prob, int size) {
    if (x == NULL)
        return err_msg(-1, 0, "binom_ds_set: argument is null");
    if (prob < 0 || prob > 1)
        return err_msg(-1, 0, "binom_ds_set: prob=%f must be in [0,1]",
                prob);

    if (size < 0)
        return err_msg(-1, 0, "binom_ds_set: size=%i must be >= 0",
                size);

    ++size; // include 0, ..., size
    x->prob = prob;
    x->size = size;

    int k = 0;
    double *t = pmalloc(size * sizeof(double), "binom_ds_set");
    if (t == NULL)
        return -1;
    for (k = 0; k < size; ++k) {
        t[k] = binom_choose(size, k) * pow(prob, k) * pow(1 - prob, size - k);
        if (t[k] < 0)
            return err_msg(-1, 0, "binom_ds_set: failed to get prob");
    }

    if (cat_ds_set_p(&x->catd, t, size) < 0)
        return err_msg(-1, 0, "binom_ds_set: failed to set parameters");

    free(t);

    return 0;
}

int binom_ds_rv(binom_ds_t *x) {
    if (x == NULL)
        return err_msg(-1, 0, "binom_ds_rv: argument is null");

    if (x->size == -1)
        return err_msg(-1, 0, "binom_ds_rv: distribution not initialized");

    int ret;
    int r = cat_ds_rv(&x->catd, &ret);
    if (ret < 0)
        return err_msg(-1, 0, "failed to sample random variable");

    return r;
}

int *binom_ds_sample(binom_ds_t *x, int *len) {
    *len = -1;
    if (x == NULL) {
        err_msg(-1, 0, "binom_ds_sample: argument is null");
        return NULL;
    }

    if (x->size < 0) {
        err_msg(-1, 0, "binom_ds_sample: distribution not initialized");
        return NULL;
    }

    int i, rv = binom_ds_rv(x);
    *len = rv;
    if (rv == 0)
        return NULL;

    int *p = pmalloc(rv * sizeof(int), "binom_ds_sample");
    if (p == NULL) {
        err_msg(-1, 0, "binom_ds_sample: failed to allocate memory");
        return NULL;
    }

    bflg_t flgs;
    bflg_init_empty(&flgs);
    bflg_init(&flgs, x->size);

    int xs = x->size-1;
    for (i = 0; i < rv; ++i) {
        int r_ix;
        do {
            r_ix = rand() % xs;
        } while (bflg_get(&flgs, r_ix));
        bflg_set(&flgs, r_ix);
        p[i] = r_ix;
    }
    *len = rv;
    bflg_free(&flgs);

    return p;
}

/*******************************************************************************
 * seq_allele_t
 ******************************************************************************/

void seq_allele_init(seq_allele_t *sa){
    if (sa == NULL) return;
    sa->bcf1 = NULL;
    sa->allele = 0;
}

/*******************************************************************************
 * genomic position/sequence range
 ******************************************************************************/

void int_range_init(int_range_t *range){
    if (range == NULL) return;
    range->beg = range->end = -1;
}

int int_range_valid(int_range_t range) {
    return range.beg >= 0 && range.beg <= range.end;
}

void seq_range_init(seq_range_t *seq_rng){
    if (seq_rng == NULL) return;
    int_range_init(&seq_rng->range);
    seq_rng->seq = NULL;
    mv_init(&seq_rng->av);
    seq_rng->rc = 0;
}

int int_range_subset(int_range_t *int_rng, int_range_t *out, size_t pos, 
        size_t len){
    if (int_rng == NULL || out == NULL)
        return err_msg(-1, 0, "int_range_subset: argument is null");

    size_t q_beg = pos, q_end = pos + len, q_len = len;
    int ir_len = int_rng->end - int_rng->beg;
    assert(ir_len > 0);
    size_t ir_ulen = (size_t)ir_len;

    // check subset fits within range
    if (q_len > ir_ulen)
        return err_msg(-1, 0, "int_range_subset: requested subset [%u,%u) "
                "extends beyond end of range %u", q_beg, q_end, ir_ulen);
    assert(q_beg < ir_ulen);

    out->beg = int_rng->beg + q_beg;
    out->end = out->beg + len;

    return 0;
}

void seq_range_free(seq_range_t *seq_rng){
    if (seq_rng == NULL) return;
    int_range_init(&seq_rng->range);
    free(seq_rng->seq);
    mv_free(&seq_rng->av);
    seq_rng->rc = 0;
}

int seq_range_subset(seq_range_t *seq_rng, seq_range_t *out, 
        size_t pos, size_t len){
    if (seq_rng == NULL || out == NULL)
        return err_msg(-1, 0, "seq_range_subset: argument is null");

    size_t q_beg = pos, q_end = pos + len, q_len = len;
    int sr_len = seq_rng->range.end - seq_rng->range.beg;
    assert(sr_len > 0);
    size_t sr_ulen = (size_t)sr_len;
    assert(strlen(seq_rng->seq) == sr_ulen);
    // fprintf(stdout, "seq_range_subset: q_beg=%zu q_len=%zu sr_len=%zu\n", q_beg, q_len, sr_ulen);
    
    if (q_end > sr_ulen)
        return err_msg(-1, 0, "seq_range_subset: requested subset [%u,%u) "
                "extends beyond end of range %u", q_beg, q_end, sr_ulen);
    assert(q_beg < sr_ulen);

    int_range_t r_out = {0,0};
    r_out.beg = seq_rng->range.beg + q_beg;
    r_out.end = r_out.beg + q_len;

    char *seq = pcalloc(q_len + 1, sizeof(char), "seq_range_subset");
    if (seq == NULL)
        return -1;

    memcpy(seq, seq_rng->seq + q_beg, q_len * sizeof(char));
    seq[q_len] = '\0';

    out->range = r_out;
    out->seq = seq;

    // add variants
    size_t i, n_var = mv_size(&seq_rng->av);
    for (i = 0; i < n_var; ++i){
        seq_allele_t sa = mv_i(&seq_rng->av, i);
        int var_pos = sa.bcf1->pos;
        if (var_pos >= r_out.beg && var_pos < r_out.end)
            mv_push(seq_allelev, &out->av, sa);
    }
    
    return 0;
}

// get overlapping variants
int seq_range_get_vars(seq_range_t *seq_rng, const char *chr, g_var_t *gv){
    if (seq_rng == NULL || gv == NULL)
        return err_msg(-1, 0, "seq_range_get_vars: argument is null");

    int n_vars = 0;
    ml_t(vcfr_list) vars;
    ml_init(vcfr_list, &vars);

    int32_t beg = seq_rng->range.beg, end = seq_rng->range.end;
    int rret = g_var_get_region_vars(gv, chr, beg, end, &vars);
    // if chromosome is not found, skip without warning
    if (rret == -1) {
        ml_free(vcfr_list, &vars);
        return n_vars;
    }
    if (rret < -1) {
        ml_free(vcfr_list, &vars);
        return err_msg(-1, 0, "seq_range_get_vars: failed to get regional variants");
    }

    // add seq_allele_t objects of overlapping variants
    ml_node_t(vcfr_list) *node;
    for (node = ml_begin(&vars); node; node = ml_node_next(node)){
        var_t var = ml_node_val(node);
        if (is_biallelic_snp(var.b) < 0){
            continue; // only support biallelic snps for now.
        }
        seq_allele_t sa;
        seq_allele_init(&sa);
        assert(var.b != NULL); // make sure we don't add an empty variant
        sa.bcf1 = var.b;
        if (mv_push(seq_allelev, &seq_rng->av, sa) < 0)
            return err_msg(-1, 0, "seq_range_get_vars: failed to push seq_allele "
                    "to range");
        ++n_vars;
    }
    ml_free(vcfr_list, &vars);

    return n_vars;
}

// sample genotypes from an individual
int seq_range_sample_allele(seq_range_t *seq_rng, bcf_hdr_t *vcf_hdr, int rsam){
    if (seq_rng == NULL || vcf_hdr == NULL)
        return err_msg(-1, 0, "seq_range_sample_allele: argument is null");

    int nsam = bcf_hdr_nsamples(vcf_hdr);
    if (rsam < 0 || rsam > nsam)
        return err_msg(-1, 0, "seq_range_sample_allele: rsam=%i must be "
                "within [0,%i)", rsam, nsam);

    size_t i, n_seq_vars = mv_size(&seq_rng->av);
    for (i = 0; i < n_seq_vars; ++i){
        seq_allele_t *sa = &mv_i(&seq_rng->av, i);
        assert(sa->bcf1);

        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt;
        ngt = bcf_get_genotypes(vcf_hdr, sa->bcf1, &gt_arr, &ngt_arr);
        int ploidy = ngt / nsam;
        if (ploidy != 2) err_msg(-1, 0, "ploidy=%i!!!", ploidy);

        // if rsam == nsam, then sample is ambient and get random sample
        if (rsam == nsam) rsam = rand() % nsam;
        int32_t *ptr = gt_arr + rsam * ploidy;

        // sample a homologous chromosome
        int chr_rand = rand() % ploidy;
        int is_miss = bcf_gt_is_missing(ptr[chr_rand]) || 
                (ptr[chr_rand]==bcf_int32_vector_end);

        // if allele is missing, sample an allele from another sample
        int n_tries = 0;
        while (is_miss){
            if (n_tries > 1000){
                free(gt_arr);
                return err_msg(-1, 0, "could not sample allele to replace "
                        "missing genotype for sample %s", 
                        vcf_hdr->samples[rsam]);
            }
            rsam = rand() % nsam;
            ptr = gt_arr + rsam * ploidy;
            chr_rand = rand() % ploidy;
            is_miss = bcf_gt_is_missing(ptr[chr_rand]) || 
                (ptr[chr_rand]==bcf_int32_vector_end);
            ++n_tries;
        }
        int32_t allele_index = bcf_gt_allele(ptr[chr_rand]);

        // type check
        if (allele_index < 0){
            free(gt_arr);
            return err_msg(-1, 0, "allele index < 0");
        }
        if (allele_index > (UINT8_MAX-1)){
            free(gt_arr);
            return err_msg(-1, 0, "allele index too large");
        }

        // set allele seq_allele
        uint8_t allele_uix = (uint8_t)allele_index;
        sa->allele = allele_uix;

        free(gt_arr);
    }

    return 0;
}

int seq_range_set_allele_seq(seq_range_t *seq_rng, bcf_hdr_t *vcf_hdr){
    if (seq_rng == NULL || vcf_hdr == NULL)
        return err_msg(-1, 0, "seq_range_set_allele_seq: argument is null");

    // make sure av is properly initialized
    size_t i, n_seq_vars = mv_size(&seq_rng->av);
    for (i = 0; i < n_seq_vars; ++i){
        seq_allele_t *sa = &mv_i(&seq_rng->av, i); // take a copy since same pointer
        assert(sa->bcf1);

        int var_pos = sa->bcf1->pos;
        assert(var_pos >= seq_rng->range.beg && var_pos < seq_rng->range.end);

        int pos_ix = var_pos - seq_rng->range.beg;
        assert((size_t)pos_ix < strlen(seq_rng->seq)); // TODO: remove after checking since slow

        char *allele = sa->bcf1->d.allele[sa->allele];
        assert(strlen(allele) == 1);
        seq_rng->seq[pos_ix] = allele[0];
    }
    return 0;
}

int seq_range_seq_error(seq_range_t *seq_rng, double prob) {
    if (seq_rng == NULL)
        return err_msg(-1, 0, "seq_range_seq_error: argument is null");

    if (seq_rng->seq == NULL)
        return err_msg(-1, 0, "seq_range_seq_error: seq_rng->seq is null");

    if (prob < 0 || prob > 1)
        return err_msg(-1, 0, "seq_range_seq_error: prob=%f must be within "
                "[0,1]", prob);

    char *p, base;
    for (p = seq_rng->seq; *p != '\0'; ++p)
        if ((base = base_seq_error(prob)) != 0) *p = base;

    return 0;
}

/*******************************************************************************
 * genomic position/sequence ranges vector
 ******************************************************************************/

void int_ranges_init(int_ranges_t *ranges){
    if (ranges == NULL) return;
    mv_init(&ranges->rv);
    ranges->len = 0;
}

void int_ranges_free(int_ranges_t *ranges){
    if (ranges == NULL) return;
    mv_free(&ranges->rv);
    ranges->len = 0;
}

int int_ranges_add_range(int_ranges_t *ranges, int_range_t range){
    if (ranges == NULL)
        return err_msg(-1, 0, "int_ranges_add_range: argument is null");

    int range_len = range.end - range.beg;
    assert(range_len >= 0);

    if (mv_push(int_rangev, &ranges->rv, range) < 0)
        return err_msg(-1, 0, "int_ranges_add_range: failed to push range");

    ranges->len += range_len;

    return 0;
}

int int_ranges_subset(int_ranges_t *int_rngs, int_ranges_t *out, 
        int pos, size_t len){
    if (int_rngs == NULL || out == NULL)
        return err_msg(-1, 0, "int_ranges_subset: argument is null");

    assert(pos >= 0);
    size_t q_beg = (size_t)pos;
    size_t q_len = len;
    size_t q_end = q_beg + q_len;

    if (q_end > (size_t)int_rngs->len)
        return err_msg(-1, 0, "int_ranges_subset: end of subset is greater "
                "than the requested length");

    size_t i, ir_beg = 0, ir_end = 0, n_rngs = mv_size(&int_rngs->rv);
    for (i = 0; i < n_rngs; ++i){
        int_range_t this_ir = mv_i(&int_rngs->rv, i);
        int ir_len = this_ir.end - this_ir.beg;
        assert(ir_len >= 0);
        size_t ir_ulen = (size_t)ir_len;
        ir_end += ir_ulen;

        if (q_beg >= ir_end){
            ir_beg += ir_ulen;
            continue;
        }

        if (ir_beg >= q_end)
            break;

        size_t put_ulen = ir_end < q_end ? (ir_end - q_beg) : (q_end - q_beg);
        assert(put_ulen <= len);
        assert(q_beg >= ir_beg);
        size_t ir_q_beg_ix = q_beg - ir_beg;
        // get the subset for this range
        int_range_t r_sbs;
        int_range_init(&r_sbs);
        int iret = int_range_subset(&this_ir, &r_sbs, ir_q_beg_ix, put_ulen);
        if (iret < 0) return -1;
        iret = int_ranges_add_range(out, r_sbs);
        if (iret < 0) return -1;

        q_beg += put_ulen;
        q_len -= put_ulen;
        ir_beg += ir_ulen;
    }
    
    return 0;
}

void seq_ranges_init(seq_ranges_t *seq_rngs){
    if (seq_rngs == NULL) return;

    mv_init(&seq_rngs->rv);
    seq_rngs->len = 0;
    seq_rngs->n_var = 0;
}

void seq_ranges_free(seq_ranges_t *seq_rngs){
    if (seq_rngs == NULL) return;
    // free the allocated char arrays in each seq_range
    int i;
    int n_ranges = mv_size(&seq_rngs->rv);
    for (i = 0; i < n_ranges; ++i){
        seq_range_free(&mv_i(&seq_rngs->rv, i));
    }
    // free the nodes in the list
    mv_free(&seq_rngs->rv);
    seq_rngs->len = 0;
    seq_rngs->n_var = 0;
}

void seq_ranges_dstry(seq_ranges_t *seq_rngs){
    if (seq_rngs == NULL) return;
    seq_ranges_free(seq_rngs);
    free(seq_rngs);
}

int seq_ranges_add_range(seq_ranges_t *seq_rngs, seq_range_t seq_rng){
    if (seq_rngs == NULL)
        return err_msg(-1, 0, "seq_ranges_add_range: argument is null");

    int range_len = seq_rng.range.end - seq_rng.range.beg;
    assert(range_len >= 0);

    if (mv_push(seq_rangev, &seq_rngs->rv, seq_rng) < 0)
        return err_msg(-1, 0, "seq_ranges_add_range: failed to push range");

    seq_rngs->len += range_len;

    return 0;
}

int seq_ranges_subset(seq_ranges_t *seq_rngs, seq_ranges_t *out,
        int pos, size_t len){
    if (seq_rngs == NULL || out == NULL)
        return err_msg(-1, 0, "seq_ranges_subset: argument is null");

    assert(pos >= 0);
    size_t q_beg = (size_t)pos;
    size_t q_len = len;
    size_t q_end = q_beg + q_len;
    size_t r_len = (size_t)(seq_rngs->len);

    if (q_end > r_len) {
        fprintf(stderr, "seq_ranges_subset: invalid subset of length %zu "
                "at position %i in sequence of length %zu\n",
                len, pos, r_len);
        return -1;
    }

    int n_vars = 0;
    size_t i;
    size_t sr_beg = 0, sr_end = 0; // begin and end index of sequence range
    for (i = 0; i < mv_size(&seq_rngs->rv); ++i){
        seq_range_t sr = mv_i(&seq_rngs->rv, i);
        int sr_len = sr.range.end - sr.range.beg;
        assert(sr_len >= 0);
        size_t sr_ulen = (size_t)sr_len;
        sr_end += sr_ulen;

        if (q_beg >= sr_end){
            sr_beg += sr_ulen;
            continue;
        }

        if (sr_beg >= q_end)
            break;

        // q_beg < sr_end
        // get end of range
        size_t put_ulen = sr_end < q_end ? (sr_end - q_beg) : (q_end - q_beg);
        assert(put_ulen <= len);
        assert(q_beg >= sr_beg);
        size_t sr_q_beg_ix = q_beg - sr_beg;
        // get the subset for this range
        seq_range_t r_sbs;
        seq_range_init(&r_sbs);
        int sret = seq_range_subset(&sr, &r_sbs, sr_q_beg_ix, put_ulen);
        if (sret < 0) return -1;
        // note char array is shallow copied, don't free after
        sret = seq_ranges_add_range(out, r_sbs);
        if (sret < 0) return -1;

        q_beg += put_ulen;
        q_len -= put_ulen;
        sr_beg += sr_ulen;
        n_vars += mv_size(&r_sbs.av);
    }
    out->n_var = n_vars;

    return 0;
}

int seq_ranges_var(seq_ranges_t *seq_rngs, const char *chr, g_var_t *gv){
    if (seq_rngs == NULL || gv == NULL)
        return err_msg(-1, 0, "seq_ranges_var: argument is null");

    int n_vars = 0;
    size_t i, n_rng = mv_size(&seq_rngs->rv);
    for (i = 0; i < n_rng; ++i){
        int vret = seq_range_get_vars(&mv_i(&seq_rngs->rv, i), chr, gv);
        if (vret < 0)
            return -1;
        n_vars += vret;
    }
    seq_rngs->n_var = n_vars;
    return 0;
}

int seq_ranges_sample_allele(seq_ranges_t *seq_rngs, bcf_hdr_t *vcf_hdr, 
        int rsam){
    if (seq_rngs == NULL || vcf_hdr == NULL)
        return err_msg(-1, 0, "seq_ranges_sample_allele: argument is null");

    size_t i, n_rng = mv_size(&seq_rngs->rv);
    for (i = 0; i < n_rng; ++i){
        int vret = seq_range_sample_allele(&mv_i(&seq_rngs->rv, i), vcf_hdr, 
                rsam);
        if (vret < 0)
            return -1;
    }
    return 0;
}

int seq_ranges_set_allele_seq(seq_ranges_t *seq_rngs, bcf_hdr_t *vcf_hdr){
    if (seq_rngs == NULL || vcf_hdr == NULL)
        return err_msg(-1, 0, "seq_ranges_set_allele_seq: argument is null");

    size_t i, n_rng = mv_size(&seq_rngs->rv);
    for (i = 0; i < n_rng; ++i){
        int vret = seq_range_set_allele_seq(&mv_i(&seq_rngs->rv, i), vcf_hdr);
        if (vret < 0)
            return -1;
    }
    return 0;
}

int seq_ranges_seq_error(seq_ranges_t *seq_rngs, double prob) {
    if (seq_rngs == NULL)
        return err_msg(-1, 0, "seq_ranges_seq_error: argument is null");

    if (prob < 0 || prob > 1)
        return err_msg(-1, 0, "seq_ranges_seq_error: prob=%f must be within "
                "[0,1]", prob);

    size_t i, n_v = mv_size(&seq_rngs->rv);
    for (i = 0; i < n_v; ++i) {
        if (seq_range_seq_error(&mv_i(&seq_rngs->rv, i), prob) < 0)
            return -1;
    }

    return 0;
}

int seq_ranges_check_len(seq_ranges_t *seq_rngs) {
    if (seq_rngs == NULL)
        return err_msg(-1, 0, "seq_ranges_check_len: argument is null");
    size_t i, n_v = mv_size(&seq_rngs->rv), seq_len = 0, str_len = 0;
    for (i = 0; i < n_v; ++i) {
        seq_range_t sr = mv_i(&seq_rngs->rv, i);
        int sr_len = sr.range.end - sr.range.beg;
        seq_len += sr_len;
        str_len += strlen(sr.seq);
    }
    if (seq_len != str_len)
        return 0;
    else
        return 1;
}

/*******************************************************************************
 * chromosome
 ******************************************************************************/

int str_num_n(const char *str, size_t len) {
    int num_n = 0;

    size_t ix;
    for (ix = 0; ix < len && *str; ++ix, ++str) {
        if (*str == 'N' || *str == 'n')
            ++num_n;
    }
    if (ix != len)
        return err_msg(-1, 0, "str_num_n: str '%s' has length %zu "
                "but param len=%zu", ix, len);

    return num_n;
}

fa_seq_t *fa_seq_alloc(){
    fa_seq_t *cs = pmalloc(sizeof(fa_seq_t), "fa_seq_alloc");
    if (cs == NULL)
        return NULL;
    cs->fai = NULL;
    cs->chrm_prob = cat_ds_alloc();
    if (cs->chrm_prob == NULL) {
        free(cs);
        return NULL;
    }
    mv_init(&cs->seqs);
    cs->c_names = init_str_map();
    if (cs->c_names == NULL) {
        cat_ds_dstry(cs->chrm_prob);
        free(cs);
        return NULL;
    }

    cs->seq_n_bp = NULL;
    return cs;
}

void fa_seq_dstry(fa_seq_t *cs){
    if (cs == NULL) return;
    fai_destroy(cs->fai);
    cat_ds_dstry(cs->chrm_prob);
    size_t i = 0;
    for (i = 0; i < mv_size(&cs->seqs); ++i)
        free(mv_i(&cs->seqs, i));
    mv_free(&cs->seqs);
    destroy_str_map(cs->c_names);
    iregs_dstry(cs->seq_n_bp);
    free(cs);
}

int fa_seq_add_fai(fa_seq_t *cs, const char *fa_fn){
    if (cs == NULL) return(-1);

    // load sequence
    cs->fai = fai_load(fa_fn);
    if (cs->fai == NULL)
        return err_msg(-1, 0, "fa_seq_add_fai: failed to load index "
                "from '%s'", fa_fn);
    int n_seq = faidx_nseq(cs->fai);
    if (n_seq == 0)
        return err_msg(-1, 0, "fa_seq_add_fai: no chromosomes "
                "found in '%s'", fa_fn);

    double *probs = pmalloc(n_seq * sizeof(double), "fa_seq_add_fai");
    if (probs == NULL)
        return -1;
    int i;
    for (i = 0; i < n_seq; ++i){
        const char *seqn = faidx_iseq(cs->fai, i);
        probs[i] = (double)faidx_seq_len(cs->fai, seqn);
        if (probs[i] < 0)
            return err_msg(-1, 0, "fa_seq_add_fai: seq len return "
                    "negative length for contig %i", i);
        int found = 0;
        // add chromosome name
        if (add2str_map(cs->c_names, seqn, &found) < 0)
            return -1;
        if (found)
            return err_msg(-1, 0, "fa_seq_add_fai: contig '%s' found twice", 
                    seqn);
    }
    // set probabilities by length
    int cret = cat_ds_set_p(cs->chrm_prob, probs, (uint64_t)n_seq);
    if (cret < 0)
        return -1;
    free(probs);

    return 0;
}

int fa_seq_load_seq(fa_seq_t *cs){
    if (cs == NULL)
        return err_msg(-1, 0, "fa_seq_load_seq: argument is null");

    int i, n_seq = faidx_nseq(cs->fai);
    for (i = 0; i < n_seq; ++i){
        const char *seqn = faidx_iseq(cs->fai, i);
        int chrm_len = faidx_seq_len(cs->fai, seqn);
        int end = chrm_len - 1, qlen = 0;
        char *seq = faidx_fetch_seq(cs->fai, seqn, 0, end, &qlen);
        if (qlen == -2)
            return err_msg(-1, 0, "fa_seq_load_seq: contig '%s' not found",
                   seqn);
        if (qlen == -1)
            return err_msg(-1, 0, "fa_seq_load_seq: could not fetch sequence "
                    "for contig '%s'", seqn);
        if (seq == NULL)
            return err_msg(-1, 0, "fa_seq_load_seq: %s", strerror(errno));
        if (mv_push(seq_v, &cs->seqs, seq) < 0)
            return err_msg(-1, 0, "fa_seq_load_seq: failed to add sequence to vec");
    }

    return 0;
}

int fa_seq_get_seq_n_bp(fa_seq_t *fs) {
    if (fs == NULL)
        return err_msg(-1, 0, "fa_seq_get_seq_n_bp: argument is null");

    if (mv_size(&fs->seqs) == 0)
        return err_msg(-1, 0, "fa_seq_get_seq_n_bp: no sequences loaded");

    fs->seq_n_bp = fa_seq_char_ranges(fs, 'N', 1, -1, 0);
    if (fs->seq_n_bp == NULL)
        return -1;

    return 0;
}

int fa_seq_rand_chr(fa_seq_t *fa, const char **seqn, int *chrm_len) {
    if (fa == NULL || seqn == NULL || chrm_len == NULL)
        return err_msg(-1, 0, "fa_seq_rand_chr: argument is null");

    if (fa->chrm_prob == NULL || fa->chrm_prob->n == 0)
        return err_msg(-1, 0, "fa_seq_rand_chr: no chromosomes added");

    // pick a random chromosome
    int dret = 0;
    int c_ix = cat_ds_rv(fa->chrm_prob, &dret);
    if (dret < 0)
        return -1;

    // get chromosome name
    *seqn = faidx_iseq(fa->fai, c_ix);
    if (*seqn == NULL)
        return err_msg(-1, 0, "fa_seq_rand_chr: chromosome for "
                              "index %i not found", c_ix);

    // get chromosome length
    *chrm_len = faidx_seq_len(fa->fai, *seqn);
    if (*chrm_len < 0)
        return err_msg(-1, 0, "fa_seq_rand_chr: seq len for contig "
                              "%s not found", *seqn);

    return c_ix;
}

int fa_seq_rand_range(fa_seq_t *fa, int len, char strand, g_region *reg){
    if (fa == NULL || reg == NULL)
        return err_msg(-1, 0, "fa_seq_rand_range: argument is null");

    if (fa->chrm_prob == NULL || fa->chrm_prob->n == 0)
        return err_msg(-1, 0, "fa_seq_rand_range: no chromosomes added");

    // pick a random chromosome
    int n_tries = 0, max_tries = 30;
    int c_ix;
    int chrm_len = 0;
    const char *seqn;
    while (chrm_len <= 0 && n_tries < max_tries) {
        c_ix = fa_seq_rand_chr(fa, &seqn, &chrm_len);
        if (c_ix < 0)
            return -1;
        chrm_len -= len;
        ++n_tries;
    }
    if (n_tries == max_tries)
        return err_msg(-1, 0, "fa_seq_rand_range: could not get a chromosome "
                "long enough to fit range of length %i", len);

    // get a random range in the chromosomoe
    int pos = rand() % chrm_len;
    reg->strand = strand;
    reg->rid = c_ix;
    reg->start = pos;
    reg->end = pos + len;
    return 0;
}

int fa_seq_range_valid(fa_seq_t *fa, const char *c_name, int beg, int end){
    if (fa == NULL || c_name == NULL)
        return err_msg(-1, 0, "fa_seq_range_valid: argument is null");

    if (beg < 0 || end < 0)
        return err_msg(-1, 0, "fa_seq_range_valid: beg=%i end=%i "
                "is negative", beg, end);

    // if chromosome not present
    if (faidx_has_seq(fa->fai, c_name) == 0)
        return 0;

    // if beg or end outside of chromosome length
    int chrm_len = faidx_seq_len(fa->fai, c_name);

    if (beg >= chrm_len || end >= chrm_len)
        return 0;

    return 1;
}

iregs_t *fa_seq_char_ranges(const fa_seq_t *fs, char base, int nmin, int nmax,
                       int case_sensitive) {
    iregs_t *iregs = iregs_init_empty();
    if (iregs == NULL)
        return NULL;

    base = case_sensitive ? base : toupper(base);

    // loop over chromosomes in `fs`
    int i;
    for (i = 0; i < fs->c_names->n; ++i) {
        const char *seqn = faidx_iseq(fs->fai, i);
        const char *seq = mv_i(&fs->seqs, i);
        int chrm_len = faidx_seq_len(fs->fai, seqn);
        hts_pos_t beg, end;
        int j = 0;
        while (j < chrm_len) {
            char b = case_sensitive ? seq[j] : toupper(seq[j]);
            beg = (hts_pos_t)j;
            end = beg - 1;
            while (b == base) {
                end = (hts_pos_t)j;
                ++j;
                b = case_sensitive ? seq[j] : toupper(seq[j]);
            };
            ++j;
            int match_len = end - beg + 1;
            if (match_len >= nmin) {
                if (nmax <= 0 || match_len <= nmax) {
                    int ret;
                    char *seqn_beg = (char *)seqn;
                    size_t slen = strlen(seqn) - 1;
                    char *seqn_end = seqn_beg + slen;
                    ret = iregs_push_reg(iregs, seqn_beg, seqn_end, beg, end);
                    if (ret < 0)
                        return NULL;
                }
            }
        }
    }

    // add to hash
    if (iregs_parse_reghash(iregs) < 0) {
        iregs_dstry(iregs);
        return NULL;
    }

    return iregs;
}

int fa_seq_n_n(fa_seq_t *fa, const char *c_name, int beg, int end){
    if (fa == NULL || c_name == NULL)
        return err_msg(-1, 0, "fa_seq_n_n: argument is null");

    if (fa->seq_n_bp == NULL)
        return err_msg(-1, 0, "fa_seq_n_n: no N sequences found, "
                       "call `fa_seq_get_seq_n_bp`");

    // make sure range is valid
    int rval = 0;
    rval = fa_seq_range_valid(fa, c_name, beg, end);
    if (rval < 0)
        return -1;

    if (rval == 0)
        return err_msg(-1, 0, "fa_seq_n_n: range '%s:%i-%i' is invalid", 
                c_name, beg, end);

    int n_ovrlp = iregs_has_overlap(fa->seq_n_bp, c_name, beg, end);
    if (n_ovrlp < 0)
        return -1;
    if (n_ovrlp > 0) {
        int c_ix = str_map_ix(fa->c_names, (char *)c_name);
        size_t i, len = end - beg;
    }
    return n_ovrlp;

    int c_ix = str_map_ix(fa->c_names, (char *)c_name);

    size_t len = end - beg;
    int n_n = str_num_n(&mv_i(&fa->seqs, c_ix)[beg], len);
    if (n_n < 0)
        return -1;

    return n_n;
}

char *fa_seq_rand_seq(fa_seq_t *fa, int len, int *chrm_ix, int *beg){
    if (fa == NULL)
        return(NULL);

    // pick a random chromosome
    g_region reg;
    init_g_region(&reg);
    char *seq = NULL;
    int n_N = 1;
    do {
        int ret = fa_seq_rand_range(fa, len, '+', &reg);
        if (ret < 0)
            return NULL;
        free(seq);
        seq = strndup(mv_i(&fa->seqs, reg.rid) + reg.start, len);

        n_N = str_num_n(seq, len);
        if (n_N < 0) {
            free(seq);
            return NULL;
        }

    } while (n_N > 0);
    *chrm_ix = reg.rid;
    *beg = reg.start;
    return seq;
}

int fa_seq_seq_range(fa_seq_t *fa, const char *c_name, int_range_t range, 
        seq_range_t **seq_rng){
    if (fa == NULL || c_name == NULL || seq_rng == NULL)
        return err_msg(-1, 0, "fa_seq_seq_range: argument is null");

    // sanity check on input
    if (!int_range_valid(range))
        return err_msg(-1, 0, "fa_seq_seq_range: range [%i,%i) is invalid",
                range.beg, range.end);

    // check that range.end < chromosome length
    int64_t chrm_len = faidx_seq_len64(fa->fai, c_name);
    if (chrm_len < 0)
        return err_msg(-1, 0, "fa_seq_seq_range: chromosome %s not found in "
                "fasta", c_name);
    int chrm_ix = str_map_ix(fa->c_names, (char *)c_name);

    if (range.end >= chrm_len)
        return err_msg(-1, 0, "fa_seq_seq_range: range end %i is larger than "
                "chromosome '%i' size", range.end, chrm_len);

    // duplicate fasta sequence at [beg,end)
    int qlen = range.end - range.beg;
    char *qseq = calloc(qlen + 1, sizeof(char));
    if (qseq == NULL)
        return err_msg(-1, 0, "fa_seq_seq_range: %s", strerror(errno));

    const char *fa_seq = mv_i(&fa->seqs, chrm_ix);
    if (fa_seq == NULL)
        return err_msg(-1, 0, "fa_seq_seq_range: contig '%s' is null", c_name);
    memcpy(qseq, fa_seq + range.beg, qlen);
    qseq[qlen] = '\0';
    if (qseq == NULL) {
        err_msg(0, 1, "fa_seq_seq_range: query sequence has length 0 "
                "for region '%s:%i-%i'. Fasta contig length is '%zu",
                c_name, range.beg, range.end, strlen(fa_seq));
    }

    // create output
    seq_range_t *sr_out = malloc(sizeof(seq_range_t));
    if (sr_out == NULL)
        return err_msg(-1, 0, "fa_seq_seq_range: %s", strerror(errno));
    seq_range_init(sr_out);

    sr_out->range = range;
    sr_out->seq = qseq;

    // put range in output
    *seq_rng = sr_out;

    return 0;
}

int fa_seq_seq_ranges(fa_seq_t *fa, const char *c_name, int_ranges_t ranges, 
        seq_ranges_t **seq_rngs){
    if (fa == NULL || seq_rngs == NULL)
        return err_msg(-1, 0, "fa_seq_seq_ranges: argument is null");

    // allocate output
    seq_ranges_t *l = pmalloc(sizeof(seq_ranges_t), "fa_seq_seq_ranges");
    if (l == NULL)
        return -1;
    seq_ranges_init(l);

    // loop over each range and add to seq_rngs
    int n_ranges = (int)mv_size(&ranges.rv);
    int i;
    for (i = 0; i < n_ranges; ++i){
        int_range_t range = mv_i(&ranges.rv, i);
        seq_range_t *seq_rng = NULL;
        if (fa_seq_seq_range(fa, c_name, range, &seq_rng) < 0)
            return -1;
        if (seq_ranges_add_range(l, *seq_rng) < 0)
            return -1;
        free(seq_rng);
    }

    *seq_rngs = l;

    return 0;
}

/*******************************************************************************
 ******************************************************************************/

int test_chrm_seq(){
    srand(time(NULL));
    cat_ds_t *ds = cat_ds_alloc();
    assert(ds != NULL);
    uint64_t i, n = 11;
    double *ps = malloc(n * sizeof(ps));
    for (i = 0; i < n; ++i){
        int r = rand() % 10;
        ps[i] = (double)r;
        if (i) ps[i] = .01;
        else ps[i] = .9;
    }
    cat_ds_set_p(ds, ps, n);

    kavl_itr_t(k_av_ix) itr;
    kavl_itr_first(k_av_ix, ds->root, &itr);
    do {
        const struct av_ix_t *p = kavl_at(&itr);
        printf("i=%lu, cm_p=%f\n", p->ix, p->val);
    } while (kavl_itr_next(k_av_ix, &itr));

    printf("RAND MAX=%i\n", RAND_MAX);

    for (i = 0; i < 20; ++i){
        int dret = 0;
        uint64_t ri = cat_ds_rv(ds, &dret);
        printf("ri=%lu\n", ri);
    }

    free(ps);
    cat_ds_dstry(ds);


    // fasta file
    const char fa_fn[] = "/u/project/pajukant/malvarez/"
        "kobs_scrna/vat_multiome/data/ref/gencode41/GRCh38.primary_assembly.genome.fa";

    faidx_t *fai = fai_load(fa_fn);
    int n_seq = faidx_nseq(fai);
    printf("n_seq=%i\n", n_seq);
    int ii;
    for (ii = 0; ii < n_seq; ++ii){
        const char *seqn = faidx_iseq(fai, ii);
        int seq_len = faidx_seq_len(fai, seqn);
        printf("seq[%i] = %s [n=%i]\n", ii, seqn, seq_len);
        int qlen;
        char *seq = faidx_fetch_seq(fai, seqn, 10e3, 10e3 + 1, &qlen);
        printf("\tseq=%s\n", seq);
        free(seq);
    }
    // 248956422

    fai_destroy(fai);

    // test fai
    fa_seq_t *chrms = fa_seq_alloc();
    if (chrms == NULL) return -1;
    if (fa_seq_add_fai(chrms, fa_fn) < 0) return -1;
    printf("loading sequence\n");
    if (fa_seq_load_seq(chrms) < 0) return -1;

    printf("getting seqs\n");
    for (ii = 0; ii < 100; ++ii){
        int c_ix, pos, qlen = 20;
        char *seq = fa_seq_rand_seq(chrms, qlen, &c_ix, &pos);
        const char *seqn = faidx_iseq(chrms->fai, c_ix);
        fprintf(stdout, "%s:%i-%i [%s]\n", seqn, pos, pos + qlen, seq);
        free(seq);
    }

    printf("destroying sequence\n");
    fa_seq_dstry(chrms);

    return 0;
}


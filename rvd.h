
#ifndef RVD_H
#define RVD_H

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <errno.h>
#include <assert.h>
#include "kavl.h"
#include "htslib/faidx.h"
#include "g_list.h"
#include "region.h"
#include "variants.h"

typedef struct av_ix_t {
    uint64_t ix;
    double val, prob;
    KAVL_HEAD(struct av_ix_t) head;
} av_ix_t;
#define dcmp(p, q) (((q)->val < (p)->val) - ((p)->val < (q)->val))
KAVL_INIT2(k_av_ix, static, struct av_ix_t, head, dcmp)

uint64_t get_rand(uint64_t max);

static inline char base_rev_cmp(char n, int *err){
    *err = 0;
    switch(n){
        case 'A':
            return 'T';
        case 'a':
            return 'T';
        case 'T':
            return 'A';
        case 't':
            return 'A';
        case 'C':
            return 'G';
        case 'c':
            return 'G';
        case 'G':
            return 'C';
        case 'g':
            return 'C';
        case 'N':
            return 'N';
        case 'n':
            return 'N';
        default:
            *err = 1;
            return n;
    }
}

int seq_rev_cpl(char *seq);

char *seq_rand(size_t len);

static inline char base_rand() {
    int base = rand() % 4;
    switch (base) {
        case 0:
            return 'A';
        case 1:
            return 'T';
        case 2:
            return 'C';
        case 3:
            return 'G';
        default:
            return 'N';
    }
}

static inline char base_seq_error(double prob) {
    double r = (double)rand() / (double)RAND_MAX;
    if (r < prob)
        return base_rand();
    else
        return 0;
}

/*******************************************************************************
 * read name
 ******************************************************************************/

typedef struct {
    char *instr;
    int rn;
    char *fc;

    // counters
    int lane;
    // tile is e.g. 2403
    int tile1; // [1,2]
    int tile2; // [1,6]
    int tile3; // [1,78]
    int x; // [1000,20000]
    int y; // [1000,20000]

    // min and max parameters
    int min_lane, max_lane;
    int min_t1, max_t1;
    int min_t2, max_t2;
    int min_t3, max_t3;
    int min_spot, max_spot;

    // parameters
    uint64_t n_tiles;
    uint64_t spots_per_tile;
    uint64_t n_reads; // rough number of reads to ensure even distribution
                 // across tiles
    uint64_t step; // step per spot (x and y)
} il_qname_t;

void il_qname_set_step(il_qname_t *il_qname);
void il_qname_init(il_qname_t *il_qname);
void il_qname_free(il_qname_t *il_qname);
char *il_qname_get_name(il_qname_t *il_qname);
void il_qname_set_instr(il_qname_t *il_qname, const char *instr,
        int rn, const char *fc);

/*******************************************************************************
 * math
 ******************************************************************************/

/* compute and return x! factorial.
 * Only [0, 12] values are valid, otherwise the factorial will overflow.
 */
int factorial(int x);

/* Compute the binomial choose coefficient (n, k).
 * Negative @p or @k inputs return an error.
 * An overflow will return an error.
 * @return 0 on success, -1 on error.
 */
double binom_choose(double n, double k);

/*******************************************************************************
 * Categorical distribution
 ******************************************************************************/

/*! typdef
 * @abstract Categorical distribution
 * Stores probs as CDF. If x1=.4 and x2=.6, then t1=.4 and t2=1.
 * @p n might be greater than the number of items in the avl tree when some
 * values are 0.
 */
typedef struct {
    av_ix_t *root; // holds the parameter values
    uint64_t n;
} cat_ds_t;

/* Initialize a cat_ds_t struct */
void cat_ds_init(cat_ds_t *d);

/* Allocate and initialize a cat_ds_t struct */
cat_ds_t *cat_ds_alloc();

/* Free the memory and re-initialize fields in @p d. */
void cat_ds_free(cat_ds_t *d);

/* Free all memory in @p d. */
void cat_ds_dstry(cat_ds_t *d);

/* Set probability values in @p d.
 * The @p p array stores the non-negative probability values. The values are
 * normalized to sum to 1. The sum of values in @p p must be greater than 0.
 * @p d Pointer to cat_ds_t struct.
 * @p p Array of doubles containing the probability values of each category.
 * @p n Number of elements in @p p.
 */
int cat_ds_set_p(cat_ds_t *d, double *p, uint64_t n);

/* return a random index ix with probability p[ix]. Possible ix values are 
 * 0, ..., n-1.
 * @param d The categorical distribution.
 * @param ret Store return value.
 * @return Set @p ret to -1 on error, 0 on success.
 */
uint64_t cat_ds_rv(cat_ds_t *d, int *ret);

/* Sample from categorical with uniform distribution.
 * This returns a number from [min,max) with equal probability.
 * The minimum value @p min must be < the maximum value @p max. Otherwise,
 * the behaviour is undefined.
 */
int cat_ds_uni_rand(int min, int max);

/*******************************************************************************
 * Binomial distribution
 ******************************************************************************/

typedef struct {
    double prob;
    int size;
    cat_ds_t catd;
} binom_ds_t;

// init, alloc, free, destroy empty binom_ds_t objects.
void binom_ds_init(binom_ds_t *x);
binom_ds_t *binom_ds_alloc();
void binom_ds_free(binom_ds_t *x);
void binom_ds_dstry(binom_ds_t *x);

// set binomial parameters for this object.
// samples from 0, 1, ..., size.
int binom_ds_set(binom_ds_t *x, double prob, int size);

// randomly sample a variable from binomial distribution.
// returns an int from 0, 1, ... size.
int binom_ds_rv(binom_ds_t *x);

// sample elements from the vector (0, 1, ..., size - 1) without replacement.
// The number of elements is drawn from the binomial distribution in x.
// The size parameter is given by the 'n' field of @p x.
// Returns the samples elements in an array of size @p len.
// If x = 0 elements are sampled from the binomial, return NULL and len = 0.
// Sets len = -1 on error.
int *binom_ds_sample(binom_ds_t *x, int *len);

/*******************************************************************************
 * range types
 ******************************************************************************/

typedef struct {
    bcf1_t *bcf1;
    uint8_t allele;
} seq_allele_t;

// vector of seq_allele_t
mv_declare(seq_allelev, seq_allele_t);

typedef struct {
    int beg, end;
} int_range_t;

// vector of int_range_t
mv_declare(int_rangev, int_range_t);

typedef struct {
    mv_t(int_rangev) rv; // vector of int_range_t
    int len; // total length of vector
} int_ranges_t;

typedef struct {
    int_range_t range;
    char *seq;
    mv_t(seq_allelev) av;
    uint8_t rc; // 1 for reverse complement
} seq_range_t;

// vector of seq_range_t
mv_declare(seq_rangev, seq_range_t);

typedef struct {
    mv_t(seq_rangev) rv; // vector of seq_range_t
    int len; // total length of vector
    int n_var;
} seq_ranges_t;

/*******************************************************************************
 * seq_allele_t
 ******************************************************************************/

void seq_allele_init(seq_allele_t *sa);

/*******************************************************************************
 * int_range_t
 ******************************************************************************/

// Initialize int_range_t to -1 values.
void int_range_init(int_range_t *range);

// return 1 if range is valid, -1 if invalid
int int_range_valid(int_range_t range);

// subset the range by index
int int_range_subset(int_range_t *int_rng, int_range_t *out, size_t pos, 
        size_t len);

/*******************************************************************************
 * seq_range_t
 ******************************************************************************/

void seq_range_init(seq_range_t *seq_rng);

void seq_range_free(seq_range_t *seq_rng);

/* subset a seq_range.
 * Update the objected pointer to by 'out' with the subset.
 * The seq char array is allocated and must be freed by the caller.
 * @return -1 on error, 0 on success.
 */
int seq_range_subset(seq_range_t *seq_rng, seq_range_t *out, 
        size_t pos, size_t len);

// get overlapping variants
int seq_range_get_vars(seq_range_t *seq_rng, const char *chr, g_var_t *gv);

// sample an allele for sample 'sam'
int seq_range_sample_allele(seq_range_t *seq_rng, bcf_hdr_t *vcf_hdr, int rsam);

// set allele in sequence
int seq_range_set_allele_seq(seq_range_t *seq_rng, bcf_hdr_t *vcf_hdr);

/* Sample sequencing error(s)
 *
 * Expects seq field to be non-null, and for seq to be a null-byte terminated 
 * char array.
 * Modifies the seq field to contain random sequencing errors.
 *
 * @param prob Probability value within [0,1].
 * @return 0 on success, -1 on error.
 */
int seq_range_seq_error(seq_range_t *seq_rng, double prob);

/*******************************************************************************
 * int_ranges_t
 ******************************************************************************/

// Initialize empty int_ranges_t to empty
void int_ranges_init(int_ranges_t *ranges);

// Free int_ranges_t object
void int_ranges_free(int_ranges_t *ranges);

// add range
int int_ranges_add_range(int_ranges_t *ranges, int_range_t range);

// subset the range by index
int int_ranges_subset(int_ranges_t *int_rngs, int_ranges_t *out, 
        int pos, size_t len);

/*******************************************************************************
 * seq_ranges_t
 ******************************************************************************/

// Initialize empty seq_ranges_t to empty
void seq_ranges_init(seq_ranges_t *seq_rngs);

// Free seq_ranges_t object
void seq_ranges_free(seq_ranges_t *seq_rngs);
void seq_ranges_dstry(seq_ranges_t *seq_rngs);

/* add range
 *
 * Note that seq_rng is added by a shallow copy, it is not duplicated.
 * It is an error if seq_rngs is null.
 *
 * @return -1 on error, 0 on success.
 */
int seq_ranges_add_range(seq_ranges_t *seq_rngs, seq_range_t seq_rng);

/* subset seq_ranges_t
 * @param pos 0-based index of the start of the subset. This is not the 
 *  position in beg and end, but rather the index of the positions.
 * @return -1 on error, 0 on success.
 */
int seq_ranges_subset(seq_ranges_t *seq_rngs, seq_ranges_t *out,
        int pos, size_t len);

/* Get overlapping variants.
 *
 * Store the overlapping variants as bcf1_t pointers in the individual 
 * seq_range_t objects.
 *
 * @return 0 on success, -1 on error.
 */
int seq_ranges_var(seq_ranges_t *seq_rngs, const char *chr, g_var_t *gv);

int seq_ranges_sample_allele(seq_ranges_t *seq_rngs, bcf_hdr_t *vcf_hdr, 
        int rsam);

int seq_ranges_set_allele_seq(seq_ranges_t *seq_rngs, bcf_hdr_t *vcf_hdr);

int seq_ranges_seq_error(seq_ranges_t *seq_rngs, double prob);

int seq_ranges_check_len(seq_ranges_t *seq_rngs);

/*******************************************************************************
 * chromosome
 ******************************************************************************/

/* Get the number of occurences of 'N' or 'n' in a sequence string.
 * Return number of n occurrences, or -1 if lengths don't match. */
int str_num_n(const char *str, size_t len);

// Create different type for prob.
// given range vector, return the sequence spliced togather.

// vector of sequences
mv_declare(seq_v, char *);

typedef struct fa_seq_t {
    faidx_t *fai;
    cat_ds_t *chrm_prob;
    mv_t(seq_v) seqs;
    str_map *c_names;
    iregs_t *seq_n_bp; // genomic ranges of 'N' base pairs in fasta
} fa_seq_t;

/* Allocate memory and initialize all members to empty. */
fa_seq_t *fa_seq_alloc();

/* Destory fa_seq_t object and free all associated memory. */
void fa_seq_dstry(fa_seq_t *cs);

/* Add a fasta index from file @p fa_fn.
 * This reads the chromosome names and calculated the percent, or probability,
 * of each chromosome with respect to the length of the genome.
 * @return -1 on error, 0 on success.
 */
int fa_seq_add_fai(fa_seq_t *cs, const char *fa_fn);

/* Read fasta sequence into memory and store in @p cs.
 * @return -1 on error, 0 on success.
 */
int fa_seq_load_seq(fa_seq_t *cs);

/* Get the sequences of 'N' in the fasta.
 * Fill the `seq_n_bp` field.
 * The sequences must be filled with `fa_seq_load_seq`.
 * @return -1 on error, 0 on success.
 */
int fa_seq_get_seq_n_bp(fa_seq_t *fs);

/* Select a random chromosome and return the index number.
 * Return -1 on error.
 */
int fa_seq_rand_chr(fa_seq_t *fa, const char **seqn, int *chrm_len);

/* Get a random range of length @p len from the genome.
 * Store the result in @p reg. Store the strand from the parameter value
 * in @p strand.
 * @return -1 on error, 0 on success.
 */
int fa_seq_rand_range(fa_seq_t *fa, int len, char strand, g_region *reg);

/* return 1 if range is valid, 0 if invalid, -1 on error. */
int fa_seq_range_valid(fa_seq_t *fa, const char *c_name, int beg, int end);

/* Get the genomic ranges that contain a contiguous sequence of base pairs,
 *  such as 'N'.
 *
 * @param fs A fa_seq_t instance containing the genomic sequence.
 * @param base The character base to get.
 * @param nmin Return regions with at least this many contiguous characters.
 * @param nmax Return regions with at most this many contiguous characters.
 *  If a region has a longer contiguous sequence, it is not split and none of
 *  the region is returned. Set to <= 0 to ignore.
 * @param case_sensitive Whether to ignore the case (=0) or consider exact
 *  matches (=1).
 * @return An allocated iregs_t object that contains the base regions.
 */
iregs_t *fa_seq_char_ranges(const fa_seq_t *fs, char base, int nmin, int nmax,
                       int case_sensitive);

/* Get number of N bases in range, or return -1 on error. */
int fa_seq_n_n(fa_seq_t *fa, const char *c_name, int beg, int end);

char *fa_seq_rand_seq(fa_seq_t *fa, int len, int *chrm_ix, int *beg);

/* given contig and int_range_t, return the sequence in the fasta.
 *
 * The chromosome name in @p c_name is used to search for the sequence 
 * in the fasta.
 *
 * @param seq_rng Address of pointer. The pointer itself should be null.
 *
 * @return -1 on error, 0 on success.
 *  The value returned in @p seq_rng must be freed and destroyed.
 */
int fa_seq_seq_range(fa_seq_t *fa, const char *c_name, int_range_t range, 
        seq_range_t **seq_rng);
int fa_seq_seq_ranges(fa_seq_t *fa, const char *c_name, int_ranges_t ranges, 
        seq_ranges_t **seq_rngs);

int test_chrm_seq();

#endif // RVD_h


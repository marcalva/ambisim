
#ifndef REGION_H 
#define REGION_H

#include <stdlib.h>
#include "str_util.h"
#include "g_list.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "htslib/regidx.h"

#define count_t uint16_t
#define n_bases 16 // 2^4

/*******************************************************************************
 * g_region
 ******************************************************************************/

/*! @typedef
 * @abstrac Structure to hold read region
 * Start position must never be greater than end position.
 *
 * @field rid reference contid int ID of read pair 1
 * @field start start position (0-based) of the first base that consumes the query
 * @field end end position (0-based) of the last base that consumes the query
 * @field strand Character giving strand. '+' if plus strand, '-' if minus 
 *  strand, and '.' if missing.
 *
 * @note start and end are 0-based and half open-ended. Given the start and 
 * end, the bases that are part of the read are [start, end).
 */
typedef struct g_region {
    char strand;
    int32_t rid;
    int32_t start;
    int32_t end;
} g_region;

// hash function for a g_region
// kh_g_region_hf
static inline khint_t kh_g_region_hf(g_region reg) {
    khint_t h = (khint_t)kh_int_hash_func(reg.rid);
    khint_t sh = kh_int_hash_func(reg.start);
    khint_t eh = kh_int_hash_func(reg.end);
    h = ((h << 5) - h) + (khint_t)sh;
    h = ((h << 5) - h) + (khint_t)eh;
    return(h);
}

static inline khint_t kh_g_region_he(g_region reg1, g_region reg2){
    if (reg1.rid != reg2.rid) return(0);
    if (reg1.start != reg2.start) return(0);
    if (reg1.end != reg2.end) return(0);
    if (reg1.strand != reg2.strand) return(0);
    return(1);
}

// collision check
KHASH_INIT(kh_reg, g_region, int, 1, kh_g_region_hf, kh_g_region_he);

/* initialize members of a g_region object
 * Sets the rid, start, and end members to -1.
 */
void init_g_region(g_region *reg);

/* Set the coordinates of a g_region object.
 * The parameter start must be less than or equal to parameter end.
 */
void set_region(g_region *reg, int32_t rid, int32_t start, int32_t end, 
        char strand);

/* Compare two regions for equality.
 *
 * Comparison of two regions is done by comparing the individual members.
 * The comparison order is rid, start, end, and strand.
 *
 * @return <0 if r1 < r2, 0 if r1 = r2, and >0 if r1 > r2.
 *  The magnitude is meaningless.
 */
int regioncmp(g_region r1, g_region r2);

/* Print a g_region */
void print_g_region(FILE *f, g_region g);

/*******************************************************************************
 * g_pos
 ******************************************************************************/

typedef struct g_pos {
    char strand;
    int32_t rid;
    int32_t pos;
} g_pos;

void init_g_pos(g_pos *p);

/* Compare positions from two g_pos objects.
 *
 * Return negative value if p1 < p2, positive value if p1 > p2, and 
 * 0 if p1 = p2. Strands are compared for equality as well.
 *
 * @param p1 First g_pos object.
 * @param p2 Second g_pos object.
 *
 * @return 0 if p1 = p2.
 */
int poscmp(g_pos p1, g_pos p2);

/* Returns a string representation of the position.
 * Format is (Chr ID) %ID (Pos) %POS (Strand) %STRAND
 * Returned string must be freed by caller.
 */
char *str_g_pos(g_pos p);

/* Print g_pos object
 *
 * Prints (Chr ID) %ID (Pos) %POS (Strand) %STRAND
 */
void fprint_g_pos(FILE *f, g_pos p);

/*******************************************************************************
 * g_reg_pair
 ******************************************************************************/

typedef struct g_reg_pair {
    g_region r1;
    g_region r2;
} g_reg_pair;

void init_reg_pair(g_reg_pair *rp);

/* Create a g_reg_pair object from two g_region objects. */
g_reg_pair get_reg_pair(g_region r1, g_region r2);

/* Hash function for g_reg_pair object.
 *
 * @param p A g_reg_pair object.
 * @return khint_t from hash function.
 */
khint_t kh_reg_pair_hash(g_reg_pair p);

/* Test for equality between two g_reg_pair objects.
 *
 * The objects are the same if the first region is the same 
 * and if the second region is the same between p1 and p2.
 *
 * @param p1 g_reg_pair object.
 * @param p2 g_reg_pair object.
 *
 * @return 1 if equal, 0 if not.
 */
int kh_reg_pair_equal(g_reg_pair p1, g_reg_pair p2);

/*******************************************************************************
 * indexed region
 ******************************************************************************/

/*! @typedef
 * @abstract a compact struct to hold indices.
 *
 * @field ix An array of integers providing the index of the region.
 * @field n The number of elements in ix.
 * @field m The max numver of elements in ix.
 *
 * @note The index in ix corresponds to the index provided in iregs_t
 */
typedef struct iregn_t {
    int *ix;
    size_t n, m;
} iregn_t;

// vector of ints (for indices)
// mv_t(int_vec)
mv_declare(int_vec, int);

/* Add ix to iregn_t object.
 *
 * @return 0 on success, -1 on error.
 */
int iregn_add_ix(iregn_t *ireg, int ix);

/*! @typedef
 * @abstract Store regions accessible by overlap query or index query.
 *
 * Store g_region objects accessible by index. Fast retrieval of its 
 * index given the g_region.
 * Indexed overlap using htslib regidx.
 *
 * Regions are stored 0-based, closed (so end-inclusive). In constrast,
 * BED files are 0-based, half-open (so end-exclusive).
 *
 * @field idx Pointer to a regidx_t object from htslib. This stores 
 *  regions read from a BED file or other tab-delimited region file.
 * @field itr Pointer to a regitr_t object from htslib. This is an 
 *  iterator used internally by iregs_t to retrieve overlapping regions.
 * @field reg This stores an array of regions in g_region format.
 * @field hash Pointer to a khash_t(kh_reg) object. This maps the region 
 *  to the index of the region in reg.
 * @note the coordinates are always stored as 0-based inclusive [beg,end].
 */
typedef struct iregs_t {
    regidx_t *idx;
    regitr_t *itr;
    g_region *reg;
    int n, m;

    khash_t(kh_reg) *hash;

    str_map *chr_map;
} iregs_t;

/* Initialize or destroy an iregs_t object.
 * init returns NULL on failure. Sets fields to 0/NULL.
 */
iregs_t *iregs_init();
void iregs_dstry(iregs_t *iregs);

/* Initialize an iregs_t object. Also initialize the
 *  regidx_t `idx` field empty. This is to manually
 *  add regions, in contrast to adding from a BED file.
 *
 *  @return The allocated instance, or NULL on failure.
 */
iregs_t *iregs_init_empty();

/* return the ith region in reg.
 *
 * @return A g_region object.
 */
#define iregs_regi(iregs, i) ((iregs)->reg)[(i)]

/* Add a region to iregs_t
 *
 * If iregs is null, do nothing and return 0.
 * Adds a ireg to reg array, and updates the hash table to return the index.
 * beg and end are 0-based inclusive [beg,end]
 *
 * @return 0 on success, -1 on error.
 */
int iregs_add2reghash(iregs_t *iregs, const char *chr, int32_t beg, 
        int32_t end, char strand);

/* Add iregs from a BED file
 *
 * If iregs is null, do nothing and return 0.
 *
 * @param iregs Pointer to initialized iregs_t object.
 * @param fn Pointer to character string.
 * @return 0 on success, -1 on error.
 */
int iregs_add_bed(iregs_t *iregs, const char *fn);

/* Parse the iregs into a reghash
 *
 * @param iregs Pointer to iregs_t object.
 * @return 0 on success, -1 if no regions were added, -2 on error.
 */
int iregs_parse_reghash(iregs_t *iregs);

/* Add a region manually to the regidx `idx` instance.
 * If @p iregs is null or iregs->idx index is null, return 0.
 *
 * chr_beg points to the beginning of the character array containing
 *  the chromosome name, chr_end points to the last character of the array.
 *
 * The region is 0-based and closed, end-inclusive.
 *
 * @return 0 on success, -1 on error
 */
int iregs_push_reg(iregs_t *iregs, char *chr_beg, char *chr_end,
                   hts_pos_t beg, hts_pos_t end);

/* Return overlapping iregs given a region.
 *
 * All overlapping regions with at least one base pair overlap are returned
 * in @p overlaps, which contains the indexes of the overlapping regions 
 * after successful call.
 * The parameters @p beg and @p end are 0-based inclusive [beg,end].
 * An overlap is defined if >= 1 base pair from the query region is 
 * overlapping 
 *
 * @return The number of overlapping regions (stored in @p n as well), or
 *  -1 on error
 */
int iregs_overlap(iregs_t *iregs, const char *chr, int32_t beg, int32_t end, 
        mv_t(int_vec) *overlaps);

/* Test if a region has an overlapping region in iregs.
 *
 * @return -1 on error, 0 if no overlap, > 0 if overlap
 */
static inline 
int iregs_has_overlap(iregs_t *iregs, const char *chr, 
        int32_t beg, int32_t end){
    mv_t(int_vec) v;
    mv_init(&v);
    int ret = iregs_overlap(iregs, chr, beg, end, &v);
    mv_free(&v);
    return(ret);
}

/* write regions from iregs
 * 
 * @return 0 on success, -1 on error
 */
int iregs_write(iregs_t *pks, BGZF *fp);

#endif // REGION_H

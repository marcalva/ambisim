
#ifndef ANNO_H
#define ANNO_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "kbtree.h"
#include "str_util.h"
#include "bins.h"
#include "g_list.h"

/* gtf file parsing */

enum spl {SPLICE, UNSPLICE, AMBIG, SNA};
#define N_SPL 3

/* columns */
enum {GTF_SEQNAME, 
    GTF_SOURCE, 
    GTF_FEAT, 
    GTF_START, 
    GTF_END, 
    GTF_SCORE, 
    GTF_STRAND, 
    GTF_ATTR};

#define FEAT "transcript"
#define GENE_NAME "gene_name"
#define GENE_ID "gene_id"
#define GENE_TYPE "gene_type"
#define TX_NAME "transcript_name"
#define TX_ID "transcript_id"
#define TX_TYPE "transcript_type"
#define EXON_ID "exon_id"
#define GENE "gene"
#define TX "transcript"
#define EXON "exon"
#define UTR "utr"

/* all @p beg and @p end are [beg,end) 
 * all coordinates are 0-based coordinates */

/* exon_t structure */
typedef struct exon_t {
    int beg;
    int end;
} exon_t;

#define exon_t_cmp(p, q) ( ((q).beg < (p).beg) - ((p).beg < (q).beg) )
ml_declare(exon_list, exon_t, exon_t_cmp);

/* isoform_t structure */
typedef struct isoform_t {
    char *id;
    int beg;
    int end;
    ml_t(exon_list) exons;
    ml_t(exon_list) utrs;
    int cds_len;
} isoform_t;

#define iso_cmp(p, q) (strcmp(((p).id), ((q).id)))
KBTREE_INIT(kb_iso, isoform_t, iso_cmp);

/* gene_t structure */
typedef struct gene_t {
    char *id; // gene_id attribute in GTF
    char *name; // gene_name attribute in GTF
    char *type; // gene_type attribute in GTF
    int beg; // 0-based start (inclusive)
    int end; // 0-based end (exclusive)
    char strand;
    int chrm; // chromosome index
    int bin;
    kbtree_t(kb_iso) *bt_isoforms;
    int isoforms_n; // number of isoforms in @field isoforms
    int min_cds_len;
    int max_cds_len;
} gene_t;

/* gene_t struct 
 * use ml_type(gl)
 * */
static inline int gtf_gene_cmp(gene_t *p, gene_t *q){
    if (p->beg < q->beg) return(-1);
    if (p->beg > q->beg) return(1);
    if (p->end < q->end) return(-1);
    if (p->end > q->end) return(1);

    return( strcmp(p->id, q->id) );
}

ml_declare(gl, gene_t *, gtf_gene_cmp);

mv_declare(gv, gene_t *);

/* Chromosome structure */
typedef struct chr_genes_t {
    ml_t(gl) bins[MAX_BIN]; // array of pointers to gene_t objects. NULL if empty.
} chr_genes_t;

mv_declare(cv, chr_genes_t *);

/* main struct to hold gene-transcript-exon info
 * Given a chromosome string, get the index in @field chrms with kh_get and kh_val.
 */
typedef struct {

    str_map *gene_ix; // gene to index. g->ix has same memory as key here
    mv_t(gv) gix2gene; // gene index to gene. Points to gene in chrms lists
                       // free only the vector, not dereferenced genes

    str_map *chrm_ix; // string to index in chrms array
    mv_t(cv) chrms; // array of pointers to chr_genes_t objects.

} gene_anno_t;

/*******************************************************************************
 * GTF line struct
 ******************************************************************************/

/* Structure to hold gtf information
 * The fields attr_tag and attr_val hold linked lists of the 
 * tag and value members in attribute field. These fields 
 * are populated from @field attribute after calling parse_gtf_attr.
 */
typedef struct {
    kstring_t chrname;
    kstring_t source;
    kstring_t feature;
    int start;
    int end;
    int score;
    char strand;
    int frame;
    kstring_t attribute; // holds the space delimited attribute field from a GTF line.
    kstr_node *attr_tag; // tag strings of each attribute
    kstr_node *attr_val; // value strings of each attribute. Match @field attr_tag.
    int n_attr;
} gtf1_t;

/****************
 * Functions
 ****************/

/*******************************************************************************
 * exon_t
 ******************************************************************************/

/* Initialize an exon to -1, -1 coordinates */
void exon_init(exon_t *e);

/*******************************************************************************
 * isoform_t
 ******************************************************************************/

/* Initialize an isoform to empty data */
int isoform_init(isoform_t *iso);

/* Free the memory underlying an isoform */
void isoform_free(isoform_t *iso);

/* get coding sequence length (length of exons)
 * @return -1 on error, or length of coding sequence on success. */
int isoform_get_cds_len(isoform_t *iso);

/*******************************************************************************
 * gene_t
 ******************************************************************************/

/* Initialize an empty gene */
void gene_init(gene_t *gene);

/* free memory underlying a gene */
void gene_free(gene_t *gene);

/* get min and max coding sequence from the isoforms
 * @return 0 on success, -1 on error.
 */
int gene_get_cds_len(gene_t *gene);

/*******************************************************************************
 * chr_genes
 ******************************************************************************/

void chr_genes_init(chr_genes_t *cg);
/* Return pointer to allocated object, or NULL on error */
chr_genes_t *chr_genes_alloc();

void chr_genes_free(chr_genes_t *cg);
void chr_genes_dstry(chr_genes_t *cg);

/* fill in cds fields for genes and isoforms in chr_genes_t
 * @return 0 on success, -1 on error.
 */
int chr_genes_get_cds_len(chr_genes_t *cg);

/*******************************************************************************
 * gene_anno_t
 ******************************************************************************/

/* Initialize/allocate a gene_anno object */
// return -1 on error, 0 on success
int gene_anno_init(gene_anno_t *gene_anno);

/* Return pointer to allocated object, or NULL on error */
gene_anno_t *gene_anno_alloc();

void gene_anno_free(gene_anno_t *gene_anno);
void gene_anno_dstry(gene_anno_t *gene_anno);

/* Return a pointer to the gene_t object given a gene_id string.
 * Returns null on error.
 */
gene_t *gene_anno_get_gene(gene_anno_t *gene_anno, char *gene_id);

int isoform_collapse_utrs(isoform_t *iso);
int gene_anno_collapse_utrs(gene_anno_t *a);

int gene_anno_cds_len(gene_anno_t *a);

/****************************
 * gene_anno_t structure
 ****************************/

/* Add chromosome to the annotation object.
 *
 * @param a pointer to annotation object.
 * @param c character array of chromosome name.
 * @return index of the chromosome in @p a->chrm_ix, or -1 on error.
 */
int add_chrm(gene_anno_t *a, char *c);

/****************************
 * GTF processing
 ****************************/

/* Initialize gtf line object */
gtf1_t *gtf1_alloc();

/* Clears the memory in gtf1_t g and reset values */
void gtf1_free(gtf1_t *g);

/* Clear memory and free gtf1_t g object */
void gtf1_dstry(gtf1_t *g);

/* Add gene, isoform, or exon from GTF line to anno
 *
 * @param a pointer to annotation object
 * @param gl pointer to gtf1_t object
 * @return 0 on success, -1 on failure.
 */
int gtf_gene2anno(gene_anno_t *a, gtf1_t *gl);
int gtf_iso2anno(gene_anno_t *a, gtf1_t *gl);
int gtf_exon2anno(gene_anno_t *a, gtf1_t *gl, int is_utr);
int gtf1_to_anno(gene_anno_t *a, gtf1_t *gl);

/*******************************************************************************
 * GTF file parsing
 ******************************************************************************/

/* Parse GTF line in string and place data into gtf1_t
 *
 * The function makes successive calls to strtok_r on the 
 * string @p line. The corresponding GTF fields are populated .
 *
 * @param line string that contains the GTF line.
 * @param g pointer to gtf1_t to populate data.
 * @return 0 on success, -1 on error.
 */
int parse_gtf1(kstring_t *line, gtf1_t *g);

/* Parse GTF line attributes.
 *
 * The attributes in @p g are stored in a single string. 
 * This function tokenizes the string and stores the key-value 
 * pairs as strings.
 *
 * @return 0 on success
 */
int parse_gtf_attr(gtf1_t *g);

/* Get attribute value of key from gtf line.
 * Searches for the first occurence of key, and returns the 
 * corresponding value in the gtf line.
 *
 * @param g pointer to gtf1_t object
 * @param key char array of key of GTF attribute
 * @return char array of the first occurence of the key's value.
 *  NULL if the attribute is not found.
 */
char *get_attr_val(gtf1_t *g, char *key);

/* Test if key-value pair is present in gtf line
 *
 * Returns 1 if found, 0 if not found.
 */
int has_key_val(gtf1_t *g, char *key, char *val);

/* Read GTF annotations into gene_anno_t object.
 *
 * The GTF file must have attributes 'gene_id' and 'transcript_id'. The 
 * exons of a transcript must be listed after the transcript, and the 
 * transcripts must be listed after the genes. The basic parameter allows 
 * to filter the GTF to include only transcripts/exons that are tagged as 
 * basic.
 *
 * @param file char array containing file name of GTF.
 * @param basic specify whether to use only isoforms that are tagged as basic (basic=1).
 * @return pointer to gene_anno_t object, or NULL if failure.
 *
 * Returned object must be freed with destroy_anno()
 */
gene_anno_t *read_from_gtf(const char *file, int basic);

/* Write out gene summary data
 *
 * Takes the genes in gene_ix and writes summary data in that order.
 * This writes out chr, start, end, strand, type, name and key.
 *
 * @return 0 on success, -1 on error.
 */
int write_gene_data(BGZF *fp, gene_anno_t *anno, str_map *gene_ix);

/****************************
 * Region overlap
 ****************************/

/* Get features that overlap given region [beg, end).
 *
 * genes must point to a valid ml_t(gl) object.
 * Free the list by calling ml_free(gl, genes) only.
 *
 * @param a gene_anno_t object to retrieve features from.
 * @param ref Reference sequence name (chromosome) in character array.
 * @param beg 0-based start position of region.
 * @param end 0-based position the base pair after the end. Open end position.
 * @param stranded 0 if given region is unstranded, any other value if the region is stranded.
 * @param strand Character of strand, one of '+' or '-'. Ignored if @p stranded is 0.
 * @param genes array of pointers to gene_t objects, for overlapping genes.
 * @param genes_len Pointer to integer that contains the length of the genes array.
 * @param p Minimum fraction of the region that overlaps.
 *
 * @return number of overlapping features returned, -1 on error.
 *
 * The function will reallocate the @p genes array to fit the number of overlapping genes. 
 * If @p ref chromosome is not found in @p a, then 0 is returned and @p genes is unchanged.
 * If no genes overlap, then 0 is returned and @p genes is unchanged.
 * The object @p genes must be allocated before the function is called, and freed after the call.
 */
int feats_from_region_p(const gene_anno_t *a, const char* ref, 
        int32_t beg, int32_t end, uint8_t stranded, char strand, 
        ml_t(gl) *genes, double p);

/* Get features that overlap completely with @p set to 1 in the above company */
static inline int feats_from_region(const gene_anno_t *a, const char* ref, int32_t beg, 
        int32_t end, uint8_t stranded, char strand, ml_t(gl) *genes){
    return feats_from_region_p(a, ref, beg, end, stranded, strand, genes, 1.0);
}

/****************************
 * Summary functions
 ****************************/

int n_feat(gene_anno_t *a, int *n_gene, int *n_iso, int *n_exon);

int test_anno(char *file);

#endif //ANNO_H


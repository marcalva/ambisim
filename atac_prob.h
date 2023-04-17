
#ifndef ATAC_PROB_H
#define ATAC_PROB_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include <errno.h>
#include "gtf_anno.h"
#include "str_util.h"
#include "bc_sim.h"
#include "variants.h"
#include "region.h"
#include "rvd.h"

void atac_prob_init(atac_prob_t *ap);
void atac_prob_free(atac_prob_t *ap);

atac_prob_t *atac_prob_alloc();
void atac_prob_dstry(atac_prob_t *ap);

int atac_prob_load_peaks(atac_prob_t *sp, const char *peaks_fn);

int atac_prob_add_prob(atac_prob_t *ap, double *probs, int ncol, int nrow);
int atac_prob_load_file(atac_prob_t *ap, const char *prob_fn);

int peak_range(g_region *reg, int_ranges_t *ranges);

int atac_prob_sample_peak(atac_prob_t *ap, uint16_t k, g_region *reg);

int atac_prob_sample_read(sc_sim_t *sc_sim, uint16_t k, int peak, int rsam, atac_pair_t *pair);

/*******************************************************************************
 * atac_read_t
 ******************************************************************************/

void atac_read_init(atac_read_t *atac_read);

void atac_read_free(atac_read_t *atac_read);

int test_sample_atac(atac_prob_t *ap, g_var_t *gv, fa_seq_t *fa, 
        int n_reads, uint32_t read_len);

void atac_pair_init(atac_pair_t *atac_pair);

void atac_pair_free(atac_pair_t *atac_pair);

int atac_pair_set_name(atac_pair_t *atac_pair, il_qname_t *names);

int atac_pair_set_seq(atac_pair_t *atac_pair, const char *bc_name);

int atac_pair_seq_error(sc_sim_t *sc_sim, atac_pair_t *atac_pair);

int atac_pair_write(atac_pair_t *atac_pair, BGZF *fs[3]);

#endif // ATAC_PROB_H


#ifndef SAM_PROB_H
#define SAM_PROB_H

#include "bc_sim.h"

// allocate and destroy sam_prob_t structs
sam_prob_t *sam_prob_alloc();
void sam_prob_dstry(sam_prob_t *sp);

/* Load sample probabilities from a file.
 * First column contains sample IDs, second contains probabilities.
 * Returns 0 on success, -1 on error.
 */
int sam_prob_load_probs(sam_prob_t *sp, const char *fn);

/* Sample a sample from a sam_prob_t struct.
 * Returns the sample index on success (0-based), or -1 on error.
 */
int sam_prob_sample_sam(sam_prob_t *sp);

#endif // SAM_PROB_H

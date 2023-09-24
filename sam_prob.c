

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include <errno.h>
#include <time.h>
#include "sam_prob.h"
#include "array_util.h"

sam_prob_t *sam_prob_alloc() {
    sam_prob_t *sp = malloc(sizeof(sam_prob_t));
    if (sp == NULL) {
        err_msg(-1, 0, "sam_prob_alloc: %s", strerror(errno));
        return NULL;
    }
    sp->samples = NULL;
    sp->sam_probs = NULL;
    return sp;
}

void sam_prob_dstry(sam_prob_t *sp) {
    if (sp == NULL)
        return;

    destroy_str_map(sp->samples);
    cat_ds_dstry(sp->sam_probs);
    free(sp);
}

int sam_prob_load_probs(sam_prob_t *sp, const char *fn) {
    if (sp == NULL || fn == NULL)
        return err_msg(-1, 0, "sam_prob_load_probs: argument is null");

    int rret = 0;
    double **arr = NULL;
    char **rownames = NULL, **colnames = NULL;
    char *delim = "\t", newline = '\n';

    int nrow = 0, ncol = 0, rowcol = 1, header = 0;
    rret = read_matrix_double(fn, &arr, rowcol, header, 
            &rownames, &nrow, &colnames, &ncol, 
            delim, newline);
    if (rret < 0)
        return err_msg(-1, 0, "sam_prob_load_probs: failed to read matrix from file");

    sp->samples = init_str_map_array(rownames, nrow);
    if (sp->samples == NULL)
        return -1;

    // flatten array
    double *arr_f = malloc(nrow * sizeof(double));
    int i;
    for (i = 0; i < nrow; ++i)
        arr_f[i] = arr[i][0];

    sp->sam_probs = cat_ds_alloc();
    if (sp->sam_probs == NULL)
        return -1;

    if (cat_ds_set_p(sp->sam_probs, arr_f, nrow) < 0)
        return -1;

    for (i = 0; i < nrow; ++i)
        free(arr[i]);
    free(arr);
    for (i = 0; i < nrow; ++i)
        free(rownames[i]);
    free(rownames);
    free(arr_f);
    
    return 0;
}

int sam_prob_sample_sam(sam_prob_t *sp) {
    if (sp == NULL)
        return err_msg(-1, 0, "sam_prob_sample_sam: argument is null");

    if (sp->sam_probs == NULL || sp->samples == NULL)
        return err_msg(-1, 0, "sam_prob_sample_sam: sam_prob_t is unitialized");

    int ret = 0;
    uint64_t s = cat_ds_rv(sp->sam_probs, &ret);
    if (ret < 0)
        return -1;

    return s;
}

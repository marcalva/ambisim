
#ifndef BITS_H
#define BITS_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>

static inline size_t _bsize(size_t x){
    x = x >> 3;
    ++x;
    return(x);
}

typedef struct {
    uint8_t *f; // array of size _bsize(n) bytes
    size_t n; // number of bits
    size_t n_buckets; // length of f
} bflg_t;

#define bflg_init_empty(flg) ( (flg)->f = NULL, (flg)->n = 0, (flg)->n_buckets = 0 ) 

#define bflg_free(flg) ( free((flg)->f), (flg)->n = 0, (flg)->n_buckets = 0 )
#define bflg_unset(flg, ix) ( (flg)->f[(ix) >> 3] &= ~(1UL << ((ix) & 0x7U)) )
#define bflg_set(flg, ix) ( (flg)->f[(ix) >> 3] |= (1UL << ((ix) & 0x7U)) )
#define bflg_get(flg, ix) ( ( (flg)->f[(ix) >> 3] >> ((ix) & 0x7U) ) & 1 ) 
#define bflg_size(flg) (flg)->n

static inline int bflg_resize(bflg_t *bflg, size_t n){
    if (n <= bflg->n) return(0);
    size_t m = _bsize(n);
    bflg->f = (uint8_t *)realloc(bflg->f, m * sizeof(uint8_t));
    if (bflg->f == NULL){
        fprintf(stderr, "error: bflg_init: %s", strerror(errno));
        return(-1);
    }
    bflg->n_buckets = m;
    // set allocated bits to 0
    size_t i;
    for (i = bflg->n; i < n; ++i) bflg_unset(bflg, i);
    bflg->n = n;
    return(0);
}

static inline int bflg_init(bflg_t *bflg, size_t n){
    bflg->n = n;
    size_t m = _bsize(n);
    bflg->f = (uint8_t *)realloc(bflg->f, m * sizeof(uint8_t));
    if (bflg->f == NULL){
        fprintf(stderr, "error: bflg_init: %s", strerror(errno));
        return(-1);
    }
    memset(bflg->f, 0, m * sizeof(uint8_t));
    bflg->n_buckets = m;
    return(0);
}

#endif // BITS_H


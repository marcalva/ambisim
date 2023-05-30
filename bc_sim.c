
#include "htslib/vcf.h"
#include "bc_sim.h"
#include "g_list.h"
#include "gex_prob.h"
#include "atac_prob.h"

#define ROUND_2_INT(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))

void bc_sim_init(bc_sim_t *bc_sim){
    if (bc_sim == NULL) return;

    bc_sim->rna_barcode = NULL;
    bc_sim->atac_barcode = NULL;
    bc_sim->n_cells = 0;
    mv_init(&bc_sim->samples);
    mv_init(&bc_sim->cell_types);
    bc_sim->rna_nreads_cell = 0;
    bc_sim->rna_nreads_ambn = 0;
    bc_sim->atac_nreads_cell = 0;
    bc_sim->atac_nreads_ambn = 0;
    bc_sim->atac_npeak_cell = 0;
    bc_sim->atac_npeak_ambn = 0;
}

void bc_sim_free(bc_sim_t *bc_sim){
    if (bc_sim == NULL) return;

    free(bc_sim->rna_barcode);
    free(bc_sim->atac_barcode);
    mv_free(&bc_sim->samples);
    mv_free(&bc_sim->cell_types);
    bc_sim_init(bc_sim);
}

mv_t(bc_simv) bc_sim_read_file(const char *file, str_map *samples,
        str_map *cell_types, int *ret){
    *ret = 0;
    mv_t(bc_simv) bc_sim_v;
    mv_init(&bc_sim_v);
    if (file == NULL){
        *ret = err_msg(-1, 0, "bc_sim_read_file: argument is null");
        return bc_sim_v;
    }
    // resize to help allocations
    if (mv_resize(bc_simv, &bc_sim_v, 2e6) < 0){
        *ret = err_msg(-1, 0, "bc_sim_read_file: failed to allocate vector");
        return bc_sim_v;
    }

    BGZF *fp = bgzf_open(file, "r");
    if (fp == 0){
        *ret = err_msg(-1, 0, "bc_sim_read_file: failed to open file %s", file);
        return bc_sim_v;
    }

    // local variables
    unsigned int ui;
    char newline = '\n', **tokens = NULL, delims[3] = "\t ";
    int slen = 0, sm = 0;
    char **subtokens = NULL, subdelim[2] = ",";
    int sub_slen = 0, sub_sm = 0;
    int bret = 0;

    // read out header
    kstring_t kstr = KS_INITIALIZE;
    bret = bgzf_getline(fp, newline, &kstr);

    while ((bret = bgzf_getline(fp, newline, &kstr)) >= 0){
        // tokenize
        if (split_line(kstr.s, &tokens, delims, &slen, &sm) < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: failed to tokenize '%s'", kstr.s);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }

        // check all columns are present
        if (slen < 11){
            *ret = err_msg(-1, 0, "bc_sim_read_file: line '%s' has less than 11 "
                    "columns [%i]. Make sure columns are separated by whitespace.",
                    kstr.s, slen);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }

        bc_sim_t bc_sim;
        bc_sim_init(&bc_sim);

        // get barcodes
        bc_sim.rna_barcode = strdup(tokens[0]);
        bc_sim.atac_barcode = strdup(tokens[1]);

        // get num. of cells
        int cret = 0;
        unsigned int n_cells = str2uint(tokens[2], &cret);
        if (cret < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: n_cells entry %s must be unsigned int", 
                    tokens[2]);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        bc_sim.n_cells = n_cells;

        // get sample indices
        if (split_line(tokens[3], &subtokens, subdelim, &sub_slen, &sub_sm) < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: failed to tokenize '%s'", tokens[3]);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        unsigned int sub_ulen = (unsigned int)sub_slen;
        if (n_cells && sub_ulen != n_cells){
            *ret = err_msg(-1, 0, "bc_sim_read_file: number of samples in '%s' "
                    "does not agree with number of cells '%u'", tokens[3], n_cells);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        // only add if not empty
        for (ui = 0; n_cells && ui < sub_ulen; ++ui){
            int sam_ix = str_map_ix(samples, subtokens[ui]);
            if (sam_ix < 0){
                *ret = err_msg(-1, 0, "bc_sim_read_file: could not find sample '%s' "
                        "in VCF", subtokens[ui]);
                free(tokens);
                ks_free(&kstr);
                return bc_sim_v;
            }
            uint8_t sam_uix = (uint8_t)sam_ix;
            if (mv_push(ui8v, &bc_sim.samples, sam_uix) < 0){
                *ret = err_msg(-1, 0, "bc_sim_read_file: failed to push to bc_sim.sample");
                free(tokens);
                ks_free(&kstr);
                return bc_sim_v;
            }
        }

        // get cell type indices
        // subtokenize
        if (split_line(tokens[4], &subtokens, subdelim, &sub_slen, &sub_sm) < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: failed to tokenize '%s'", tokens[4]);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        sub_ulen = (unsigned int)sub_slen;
        if (n_cells && sub_ulen != n_cells){
            *ret = err_msg(-1, 0, "bc_sim_read_file: number of cell types in '%s' "
                    "does not agree with number of cells '%u'", tokens[4], n_cells);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        // only add if not empty
        for (ui = 0; n_cells && ui < sub_ulen; ++ui){
            int k = str_map_ix(cell_types, subtokens[ui]);
            if (k < 0){
                *ret = err_msg(-1, 0, "bc_sim_read_file: could not find "
                        "cell type '%s'", subtokens[ui]);
                free(tokens);
                ks_free(&kstr);
                return bc_sim_v;
            }
            uint8_t uk = (uint16_t)k;
            if (mv_push(ui8v, &bc_sim.cell_types, uk) < 0){
                *ret = err_msg(-1, 0, "bc_sim_read_file: failed to push to bc_sim.cell_types");
                free(tokens);
                ks_free(&kstr);
                return bc_sim_v;
            }
        }

        // get rna nreads cell
        unsigned int rna_nreads_cell = str2uint(tokens[5], &cret);
        if (cret < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: rna_nreads_cell entry %s must be unsigned int", 
                    tokens[5]);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        bc_sim.rna_nreads_cell = rna_nreads_cell;

        // get rna_nreads_ambn
        unsigned int rna_nreads_ambn = str2uint(tokens[6], &cret);
        if (cret < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: rna_nreads_ambn entry %s must be unsigned int", 
                    tokens[6]);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        bc_sim.rna_nreads_ambn = rna_nreads_ambn;

        // get atac_nreads_cell
        unsigned int atac_nreads_cell = str2uint(tokens[7], &cret);
        if (cret < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: atac_nreads_cell entry %s must be unsigned int", 
                    tokens[7]);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        bc_sim.atac_nreads_cell = atac_nreads_cell;

        // get atac_nreads_ambn
        unsigned int atac_nreads_ambn = str2uint(tokens[8], &cret);
        if (cret < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: atac_nreads_ambn entry %s must be unsigned int", 
                    tokens[8]);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        bc_sim.atac_nreads_ambn = atac_nreads_ambn;

        // get atac_npeak_cell
        unsigned int atac_npeak_cell = str2uint(tokens[9], &cret);
        if (cret < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: atac_npeak_cell entry %s must be unsigned int", 
                    tokens[9]);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        bc_sim.atac_npeak_cell = atac_npeak_cell;

        // get atac_npeak_ambn
        unsigned int atac_npeak_ambn = str2uint(tokens[10], &cret);
        if (cret < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: atac_npeak_ambn entry %s must be unsigned int", 
                    tokens[10]);
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }
        bc_sim.atac_npeak_ambn = atac_npeak_ambn;

        // add to vector
        if (mv_push(bc_simv, &bc_sim_v, bc_sim) < 0){
            *ret = err_msg(-1, 0, "bc_sim_read_file: failed to push bc_sim to vec");
            free(tokens);
            ks_free(&kstr);
            return bc_sim_v;
        }

    }
    if (bret < -1){
        *ret = err_msg(-1, 0, "bc_sim_read_file: failed to read from file '%s'", file);
        return bc_sim_v;
    }

    free(tokens);
    free(subtokens);
    ks_free(&kstr);

    bgzf_close(fp);
    mv_resize(bc_simv, &bc_sim_v, mv_size(&bc_sim_v));

    return bc_sim_v;
}

/*******************************************************************************
 * sc_stream_t
 ******************************************************************************/

void sc_streams_init(sc_streams_t *scs){
    if (scs == NULL) return;
    scs->rna_r1_fs = NULL;
    scs->rna_r2_fs = NULL;
    scs->atac_r1_fs = NULL;
    scs->atac_r2_fs = NULL;
    scs->atac_r3_fs = NULL;
}

void sc_streams_close(sc_streams_t *scs) {
    if (scs == NULL) return;

    if (scs->rna_r1_fs) bgzf_close(scs->rna_r1_fs);
    scs->rna_r1_fs = NULL;
    if (scs->rna_r2_fs) bgzf_close(scs->rna_r2_fs);
    scs->rna_r2_fs = NULL;
    if (scs->atac_r1_fs) bgzf_close(scs->atac_r1_fs);
    scs->atac_r1_fs = NULL;
    if (scs->atac_r2_fs) bgzf_close(scs->atac_r2_fs);
    scs->atac_r2_fs = NULL;
    if (scs->atac_r3_fs) bgzf_close(scs->atac_r3_fs);
    scs->atac_r3_fs = NULL;
}

void sc_streams_free(sc_streams_t *scs) {
    if (scs == NULL) return;

    sc_streams_close(scs);
}

int sc_streams_set(sc_streams_t *scs, const char *dir) {
    if (scs == NULL || dir == NULL)
        return err_msg(-1, 0, "sc_streams_init: argument is null");

    // RNA
    char rna_suff[] = "/RNA/";
    char *rna_dir = strcat2(dir, rna_suff);
    if (rna_dir == NULL)
        return -1;
    if (mkpath(rna_dir, 0755) < 0)
        return err_msg(-1, 0, "sc_streams_set: failed to create '%s'", 
                rna_dir);

    char rna_r1_suff[] = "sim_S1_L001_R1_001.fastq.gz";
    char rna_r2_suff[] = "sim_S1_L001_R2_001.fastq.gz";
    char *rna_r1_fn = strcat2(rna_dir, rna_r1_suff);
    char *rna_r2_fn = strcat2(rna_dir, rna_r2_suff);
    if (rna_r1_fn == NULL || rna_r2_fn == NULL)
        return err_msg(-1, 0, "sc_streams_set: failed to allocate strings");

    if ( (scs->rna_r1_fs = bgzf_open(rna_r1_fn, "wg1")) == NULL)
        return err_msg(-1, 0, "sc_streams_set: failed to open '%s'", rna_r1_fn);
    if ( (scs->rna_r2_fs = bgzf_open(rna_r2_fn, "wg1")) == NULL)
        return err_msg(-1, 0, "sc_streams_set: failed to open '%s'", rna_r2_fn);

    free(rna_dir);
    free(rna_r1_fn);
    free(rna_r2_fn);

    // ATAC
    char atac_suff[] = "/ATAC/";
    char *atac_dir = strcat2(dir, atac_suff);
    if (atac_dir == NULL)
        return -1;
    if (mkpath(atac_dir, 0755) < 0)
        return err_msg(-1, 0, "sc_streams_set: failed to create '%s'", 
                atac_dir);

    char atac_r1_suff[] = "sim_S1_L001_R1_001.fastq.gz";
    char atac_r2_suff[] = "sim_S1_L001_R2_001.fastq.gz";
    char atac_r3_suff[] = "sim_S1_L001_R3_001.fastq.gz";
    char *atac_r1_fn = strcat2(atac_dir, atac_r1_suff);
    char *atac_r2_fn = strcat2(atac_dir, atac_r2_suff);
    char *atac_r3_fn = strcat2(atac_dir, atac_r3_suff);
    if (atac_r1_fn == NULL || atac_r2_fn == NULL || atac_r3_fn == NULL)
        return err_msg(-1, 0, "sc_streams_set: failed to allocate strings");

    if ( (scs->atac_r1_fs = bgzf_open(atac_r1_fn, "wg1")) == NULL)
        return err_msg(-1, 0, "sc_streams_set: failed to open '%s'", atac_r1_fn);
    if ( (scs->atac_r2_fs = bgzf_open(atac_r2_fn, "wg1")) == NULL)
        return err_msg(-1, 0, "sc_streams_set: failed to open '%s'", atac_r2_fn);
    if ( (scs->atac_r3_fs = bgzf_open(atac_r3_fn, "wg1")) == NULL)
        return err_msg(-1, 0, "sc_streams_set: failed to open '%s'", atac_r3_fn);

    free(atac_dir);
    free(atac_r1_fn);
    free(atac_r2_fn);
    free(atac_r3_fn);

    return 0;
}

/*******************************************************************************
 * sc_sim_t
 ******************************************************************************/

int sc_sim_init(sc_sim_t *sc_sim) {
    if (sc_sim == NULL)
        return err_msg(-1, 0, "sc_sim_init: argument is null");

    mv_init(&sc_sim->barcodes);

    sc_sim->samples = NULL;
    sc_sim->cell_types = NULL;
    sc_sim->gex_prob = gex_prob_alloc();
    if (sc_sim->gex_prob == NULL)
        return -1;
    sc_sim->atac_prob = atac_prob_alloc();
    if (sc_sim->atac_prob == NULL)
        return -1;

    sc_sim->fa = fa_seq_alloc();
    if (sc_sim->fa == NULL)
        return -1;
    sc_sim->sr = bcf_sr_init();
    if (sc_sim->sr == NULL)
        return -1;
    /* 
    sc_sim->vcf_hdr = bcf_hdr_init("r");
    if (sc_sim->vcf_hdr == NULL)
        return -1;
    */

    sc_sim->gv = g_var_alloc();
    if (sc_sim->gv == NULL)
        return -1;
    sc_sim->chrms = init_str_map();
    if (sc_sim->chrms == NULL)
        return -1;
    il_qname_init(&sc_sim->rna_names);
    il_qname_init(&sc_sim->atac_names);

    sc_sim->rna_rd_len = 0;
    sc_sim->rna_umi_len = 0;
    sc_sim->atac_rd_len = 0;
    sc_sim->seq_error = 0.0;
    binom_ds_init(&sc_sim->rna_err_prob);
    binom_ds_init(&sc_sim->atac_err_prob);
    sc_sim->rna_rd_counter = 0;
    sc_sim->atac_rd_counter = 0;
    sc_streams_init(&sc_sim->scs);
    sc_sim->out = NULL;

    return 0;
}

void sc_sim_free(sc_sim_t *sc_sim){
    if (sc_sim == NULL) return;
    size_t i, n_bc = mv_size(&sc_sim->barcodes);
    for (i = 0; i < n_bc; ++i){
        bc_sim_t *bc = &mv_i(&sc_sim->barcodes, i);
        bc_sim_free(bc);
    }
    mv_free(&sc_sim->barcodes);

    gex_prob_dstry(sc_sim->gex_prob);
    sc_sim->gex_prob = NULL;
    atac_prob_dstry(sc_sim->atac_prob);
    sc_sim->atac_prob = NULL;

    fa_seq_dstry(sc_sim->fa);
    sc_sim->fa = NULL;

    if (sc_sim->sr)
        bcf_sr_destroy(sc_sim->sr);
    sc_sim->sr = NULL;
    sc_sim->vcf_hdr = NULL;
    destroy_str_map(sc_sim->samples);
    sc_sim->samples = NULL;
    destroy_str_map(sc_sim->cell_types);
    sc_sim->cell_types = NULL;
    g_var_dstry(sc_sim->gv);
    sc_sim->gv = NULL;
    destroy_str_map(sc_sim->chrms);
    sc_sim->chrms = NULL;
    il_qname_free(&sc_sim->rna_names);
    il_qname_free(&sc_sim->atac_names);
    free(sc_sim->out);
    sc_sim->out = NULL;
    binom_ds_free(&sc_sim->rna_err_prob);
    binom_ds_free(&sc_sim->atac_err_prob);
    sc_streams_free(&sc_sim->scs);
}

sc_sim_t *sc_sim_alloc(){
    sc_sim_t *sc_sim = malloc(sizeof(sc_sim_t));
    if (sc_sim == NULL){
        err_msg(-1, 0, "sc_sim_alloc: %s", strerror(errno));
        return NULL;
    }
    if (sc_sim_init(sc_sim) < 0)
        return NULL;
    return(sc_sim);
}

void sc_sim_dstry(sc_sim_t *sc_sim){
    if (sc_sim == NULL) return;
    sc_sim_free(sc_sim);
    free(sc_sim);
}

int sc_sim_set_rd_len(sc_sim_t *sc_sim, uint32_t rna_rd_len,
        uint32_t rna_umi_len, uint32_t atac_rd_len) {
    if (sc_sim == NULL)
        return err_msg(-1, 0, "sc_sim_set_rd_len: argument is null");

    sc_sim->rna_rd_len = rna_rd_len;
    sc_sim->rna_umi_len = rna_umi_len;
    sc_sim->atac_rd_len = atac_rd_len;

    return 0;
}

int sc_sim_set_seq_error(sc_sim_t *sc_sim, double prob) {
    if (sc_sim == NULL)
        return err_msg(-1, 0, "sc_sim_set_rd_len: argument is null");

    if (prob < 0 || prob > 1)
        return err_msg(-1, 0, "sc_sim_set_rd_len: prob=%f must be within "
                "[0,1]", prob);

    sc_sim->seq_error = prob;

    if (binom_ds_set(&sc_sim->rna_err_prob, prob, sc_sim->rna_rd_len) < 0)
        return -1;
    if (binom_ds_set(&sc_sim->atac_err_prob, prob, sc_sim->atac_rd_len) < 0)
        return -1;

    return 0;
}

int sc_sim_load_vars(sc_sim_t *sc_sim, const char *vcf_fn, 
        const char *sample_fn) {
    if (sc_sim == NULL || vcf_fn == NULL)
        return err_msg(-1, 0, "sc_sim_load_vars: argument is null");

    // if fields are present destroy to avoid memory leaks.
    if (sc_sim->sr) {
        bcf_sr_destroy(sc_sim->sr);
        sc_sim->sr = NULL;
        sc_sim->vcf_hdr = NULL;
    }

    if (load_vcf(vcf_fn, ".", 0, &sc_sim->sr, &sc_sim->vcf_hdr) < 0)
        return -1;

    // if samples are given
    if (sample_fn && sub_vcf_samples(&sc_sim->vcf_hdr, sample_fn))
        return -1;

    // set sample IDs in sc_sim samples field
    int n_samples = bcf_hdr_nsamples(sc_sim->vcf_hdr);
    if (sc_sim->samples) {
        destroy_str_map(sc_sim->samples);
        sc_sim->samples = NULL;
    }
    sc_sim->samples = init_str_map_array(sc_sim->vcf_hdr->samples, n_samples);
    if (sc_sim->samples == NULL)
        return -1;

    if (sc_sim->gv) {
        g_var_dstry(sc_sim->gv);
        sc_sim->gv = NULL;
    }
    sc_sim->gv = g_var_read_vcf(sc_sim->sr, sc_sim->vcf_hdr, 0, 0);

    return 0;
}

int sc_sim_load_cell_types(sc_sim_t *sc_sim, const char *fn) {
    if (sc_sim == NULL || fn == NULL)
        return err_msg(-1, 0, "sc_sim_load_cell_types: argument is null");

    if (sc_sim->cell_types) {
        destroy_str_map(sc_sim->cell_types);
        sc_sim->cell_types = NULL;
    }
    sc_sim->cell_types = read_str_map(fn);
    if (sc_sim->cell_types == NULL)
        return -1;

    if (sc_sim->cell_types->n < 1)
        return err_msg(-1, 0, "sc_sim_load_cell_types: no cell types");

    return 0;
}

int sc_sim_load_gtf(sc_sim_t *sc_sim, const char *file, int tx_basic) {
    if (sc_sim == NULL || file == NULL)
        return err_msg(-1, 0, "sc_sim_load_gtf: argument is null");

    if (sc_sim->gex_prob == NULL)
        return err_msg(-1, 0, "sc_sim_load_gtf: sc_sim->gex_prob is null");

    if (gex_prob_load_gtf(sc_sim->gex_prob, file, tx_basic) < 0)
        return -1;

    return 0;
}

int sc_sim_load_fa(sc_sim_t *sc_sim, const char *file) {
    if (sc_sim == NULL || file == NULL)
        return err_msg(-1, 0, "sc_sim_load_fa: argument is null");

    if (sc_sim->fa == NULL)
        return err_msg(-1, 0, "sc_sim_load_fa: sc_sim->fa is null");

    if (fa_seq_add_fai(sc_sim->fa, file) < 0)
        return err_msg(-1, 0, "sc_sim_load_fa: failed to add fasta index");
    if (fa_seq_load_seq(sc_sim->fa) < 0)
        return err_msg(-1, 0, "sc_sim_load_fa: failed to load sequence");
    return 0;
}

int sc_sim_load_gex(sc_sim_t *sc_sim, const char *rho_file, 
        const char *mmrna_file, const char *gene_ids_file) {
    if (sc_sim == NULL || rho_file == NULL || mmrna_file == NULL || 
            gene_ids_file == NULL)
        return err_msg(-1, 0, "sc_sim_load_gex: argument is null");

    if (sc_sim->gex_prob == NULL)
        return err_msg(-1, 0, "sc_sim_load_gex: sc_sim->gex_prob is null");

    // load gene IDs
    if (gex_prob_load_gene_ids(sc_sim->gex_prob, gene_ids_file) < 0)
        return -1;

    // load gene expression probabilities
    int cell_k;
    if ((cell_k = gex_prob_load_prob(sc_sim->gex_prob, rho_file, 
                    mmrna_file)) < 0)
        return -1;

    return cell_k;
}

int sc_sim_load_atac(sc_sim_t *sc_sim, const char *prob_file, 
        const char *peak_file) {
    if (sc_sim == NULL || prob_file == NULL || peak_file == NULL)
        return err_msg(-1, 0, "sc_sim_load_atac: argument is null");

    if (sc_sim->atac_prob == NULL)
        return err_msg(-1, 0, "sc_sim_load_atac: atac_prob is null");

    if (atac_prob_load_peaks(sc_sim->atac_prob, peak_file) < 0)
        return -1;

    if (atac_prob_load_file(sc_sim->atac_prob, prob_file) < 0)
        return -1;

    return 0;
}

int sc_sim_read_file(sc_sim_t *sc_sim, const char *file){
    if (sc_sim == NULL || file == NULL)
        return err_msg(-1, 0, "sc_sim_read_file: argument is null");

    if (sc_sim->samples == NULL)
        return err_msg(-1, 0, "sc_sim_read_file: samples needs to be set with "
                "call to 'sc_sim_load_vars' first");

    int ret = 0;
    sc_sim->barcodes = bc_sim_read_file(file, sc_sim->samples,
            sc_sim->cell_types, &ret);
    if (ret < 0) return -1;

    // add total num of reads
    sc_sim->rna_nreads = sc_sim->atac_nreads = 0;
    size_t i, n_bc = mv_size(&sc_sim->barcodes);
    for (i = 0; i < n_bc; ++i) {
        sc_sim->rna_nreads += mv_i(&sc_sim->barcodes, i).rna_nreads_cell;
        sc_sim->rna_nreads += mv_i(&sc_sim->barcodes, i).rna_nreads_ambn;
        sc_sim->atac_nreads += mv_i(&sc_sim->barcodes, i).atac_nreads_cell;
        sc_sim->atac_nreads += mv_i(&sc_sim->barcodes, i).atac_nreads_ambn;
    }

    fprintf(stdout, "num. RNA reads = %u, num. ATAC reads = %u\n",
            sc_sim->rna_nreads, sc_sim->atac_nreads);

    return 0;
}

int sc_sim_check_k(sc_sim_t *sc_sim) {
    if (sc_sim == NULL)
        return err_msg(-1, 0, "sc_sim_check_k: argument is null");

    uint16_t k = (uint16_t)sc_sim->cell_types->n;
    if (sc_sim->gex_prob == NULL || sc_sim->atac_prob == NULL)
        return err_msg(-1, 0, "sc_sim_check_k: prob is null");

    uint16_t gk = 0, ak = 0;
    if (sc_sim->gex_prob->k > 0) gk = sc_sim->gex_prob->k - 1;
    if (sc_sim->atac_prob->k > 0) ak = sc_sim->atac_prob->k - 1;
    if (gk && gk != k)
        return err_msg(-1, 0, "sc_sim_check_k: cell types in expr probs (%u) "
                "does not match those given in cell types file (%u)", gk, k);
    if (ak && ak != k)
        return err_msg(-1, 0, "sc_sim_check_k: cell types in peak probs (%u) "
                "does not match those given in cell types file (%u)", ak, k);

    sc_sim->K = k + 1;

    return 0;
}

int sc_sim_intrs_chrms(sc_sim_t *sc_sim) {
    if (sc_sim == NULL)
        return err_msg(-1, 0, "sc_sim_intrs_chrms: argument is null");

    // check input
    if (sc_sim->gex_prob == NULL ||
        sc_sim->gex_prob->anno == NULL ||
        sc_sim->gex_prob->anno->chrm_ix == NULL)
        return err_msg(-1, 0, "sc_sim_intrs_chrms: genes not added");
    if (sc_sim->atac_prob == NULL ||
        sc_sim->atac_prob->peaks == NULL ||
        sc_sim->atac_prob->peaks->chr_map == NULL)
        return err_msg(-1, 0, "sc_sim_intrs_chrms: peaks not added");
    if (sc_sim->fa == NULL ||
        sc_sim->fa->c_names == NULL)
        return err_msg(-1, 0, "sc_sim_intrs_chrms: fasta not added");
    if (sc_sim->gv == NULL ||
        sc_sim->gv->chrm_ix == NULL)
        return err_msg(-1, 0, "sc_sim_intrs_chrms: variants not added");

    destroy_str_map(sc_sim->chrms);
    sc_sim->chrms = init_str_map();
    if (sc_sim->chrms == NULL)
        return -1;

    str_map *gex_chrms = sc_sim->gex_prob->anno->chrm_ix;
    str_map *atac_chrms = sc_sim->atac_prob->peaks->chr_map;
    str_map *fa_chrms = sc_sim->fa->c_names;
    str_map *var_chrms = sc_sim->gv->chrm_ix;
    int found = 0;
    int i;
    for (i = 0; i < gex_chrms->n; ++i) {
        char *chrm = str_map_str(gex_chrms, i);
        if (str_map_ix(atac_chrms, chrm) < 0)
            continue;
        if (str_map_ix(fa_chrms, chrm) < 0)
            continue;
        if (str_map_ix(var_chrms, chrm) < 0)
            continue;
        if (add2str_map(sc_sim->chrms, chrm, &found) < 0)
            return -1;
    }

    return 0;
}

int bc_sim_gen_reads(sc_sim_t *sc_sim, size_t bc_i) {
    if (sc_sim == NULL)
        return err_msg(-1, 0, "bc_sim_gen_reads: argument is null");

    // get bc_sim struct
    if (bc_i >= mv_size(&sc_sim->barcodes))
        return err_msg(-1, 0, "bc_sim_gen_reads: bc_i=%zu >= num. of "
                "barcodes given (%zu)", bc_i, mv_size(&sc_sim->barcodes));
    bc_sim_t *bc_sim = &mv_i(&sc_sim->barcodes, bc_i);

    // array holder for file streams
    BGZF *atac_fs[3] = {sc_sim->scs.atac_r1_fs, sc_sim->scs.atac_r2_fs,
        sc_sim->scs.atac_r3_fs};
    BGZF *gex_fs[2] = {sc_sim->scs.rna_r1_fs, sc_sim->scs.rna_r2_fs};

    uint32_t r_i;

    uint16_t n_sam = sc_sim->samples->n;
    uint16_t n_ct = sc_sim->K - 1; // num. of cell types not including ambient

    uint32_t rna_cell_rd = bc_sim->rna_nreads_cell;
    uint32_t rna_ambn_rd = bc_sim->rna_nreads_ambn;
    uint32_t atac_cell_rd_peak = bc_sim->atac_npeak_cell;
    uint32_t atac_cell_rd_out = bc_sim->atac_nreads_cell - atac_cell_rd_peak;
    uint32_t atac_ambn_rd_peak = bc_sim->atac_npeak_ambn;
    uint32_t atac_ambn_rd_out = bc_sim->atac_nreads_ambn - atac_ambn_rd_peak;

    // input checks
    if (rna_cell_rd && bc_sim->n_cells < 1)
            return err_msg(-1, 0, "bc_sim_gen_reads: %u rna cell reads given for "
                    "an empty droplet", rna_cell_rd); 
    if (bc_sim->atac_nreads_cell && bc_sim->n_cells < 1)
            return err_msg(-1, 0, "bc_sim_gen_reads: %u atac cell reads given for "
                    "an empty droplet", bc_sim->atac_nreads_cell); 
    if (bc_sim->n_cells > mv_size(&bc_sim->samples))
        return err_msg(-1, 0, "bc_sim_gen_reads: not enough samples (%u) for "
                "droplet with %u nuclei", mv_size(&bc_sim->samples), bc_sim->n_cells); 
    if (bc_sim->n_cells > mv_size(&bc_sim->cell_types))
        return err_msg(-1, 0, "bc_sim_gen_reads: not enough cell types (%u) for "
                "droplet with %u nuclei", mv_size(&bc_sim->samples), bc_sim->n_cells); 

    int n_reads = rna_cell_rd + rna_ambn_rd + atac_cell_rd_peak + atac_cell_rd_out + 
        atac_ambn_rd_peak + atac_ambn_rd_out;

    // for RNA cell reads
    for (r_i = 0; r_i < rna_cell_rd; ++r_i) {
        // sample a cell type and a sample
        int d_ix = rand() % bc_sim->n_cells;
        uint8_t s_ix = mv_i(&bc_sim->samples, d_ix);
        uint8_t ct_ix = mv_i(&bc_sim->cell_types, d_ix);

        rna_read_t rna_read;
        rna_read_init(&rna_read);
        // initialize rna read
        int rret = gex_sample_read(sc_sim, ct_ix, s_ix, &rna_read);
        if (rret < 0) return -1;
        rret = rna_read_set_name(&rna_read, &sc_sim->rna_names);
        if (rret < 0) return -1;
        rret = rna_read_set_seq(&rna_read, bc_sim->rna_barcode,
                sc_sim->rna_umi_len);
        if (rret < 0) return -1;
        if (rna_read_seq_error(sc_sim, &rna_read) < 0) return -1;
        if (rna_read_write(&rna_read, gex_fs) < 0)
            return -1;
        rna_read_free(&rna_read);
        ++sc_sim->rna_rd_counter;
    }

    // for RNA ambient reads
    for (r_i = 0; r_i < rna_ambn_rd; ++r_i) {
        uint8_t s_ix = n_sam;
        uint8_t ct_ix = n_ct;

        rna_read_t rna_read;
        rna_read_init(&rna_read);
        int rret = gex_sample_read(sc_sim, ct_ix, s_ix, &rna_read);
        if (rret < 0) return -1;
        rret = rna_read_set_name(&rna_read, &sc_sim->rna_names);
        if (rret < 0) return -1;
        rret = rna_read_set_seq(&rna_read, bc_sim->rna_barcode,
                sc_sim->rna_umi_len);
        if (rret < 0) return -1;
        if (rna_read_seq_error(sc_sim, &rna_read) < 0) return -1;
        if (rna_read_write(&rna_read, gex_fs) < 0)
            return -1;
        rna_read_free(&rna_read);
        ++sc_sim->rna_rd_counter;
    }

    // for ATAC cell peak reads
    for (r_i = 0; r_i < atac_cell_rd_peak; ++r_i) {
        // sample a cell type and a sample
        int d_ix = rand() % bc_sim->n_cells;
        uint8_t s_ix = mv_i(&bc_sim->samples, d_ix);
        uint8_t ct_ix = mv_i(&bc_sim->cell_types, d_ix);

        atac_pair_t atac_pair;
        atac_pair_init(&atac_pair);
        int rret = atac_prob_sample_read(sc_sim, ct_ix, 1, s_ix, &atac_pair);
        if (rret < 0) return -1;
        rret = atac_pair_set_name(&atac_pair, &sc_sim->atac_names);
        if (rret < 0) return -1;
        rret = atac_pair_set_seq(&atac_pair, bc_sim->atac_barcode);
        if (rret < 0) return -1;
        if (atac_pair_seq_error(sc_sim, &atac_pair) < 0) return -1;
        rret = atac_pair_write(&atac_pair, atac_fs);
        if (rret < 0) return -1;

        atac_pair_free(&atac_pair);
        ++sc_sim->atac_rd_counter;
    }

    // for ATAC cell out reads
    for (r_i = 0; r_i < atac_cell_rd_out; ++r_i) {
        // sample a cell type and a sample
        int d_ix = rand() % bc_sim->n_cells;
        uint8_t s_ix = mv_i(&bc_sim->samples, d_ix);
        uint8_t ct_ix = mv_i(&bc_sim->cell_types, d_ix);

        atac_pair_t atac_pair;
        atac_pair_init(&atac_pair);
        int rret = atac_prob_sample_read(sc_sim, ct_ix, 0, s_ix, &atac_pair);
        if (rret < 0) return -1;
        rret = atac_pair_set_name(&atac_pair, &sc_sim->atac_names);
        if (rret < 0) return -1;
        rret = atac_pair_set_seq(&atac_pair, bc_sim->atac_barcode);
        if (rret < 0) return -1;
        if (atac_pair_seq_error(sc_sim, &atac_pair) < 0) return -1;
        rret = atac_pair_write(&atac_pair, atac_fs);
        if (rret < 0) return -1;

        atac_pair_free(&atac_pair);
        ++sc_sim->atac_rd_counter;
    }

    // for ATAC ambient peak reads
    for (r_i = 0; r_i < atac_ambn_rd_peak; ++r_i) {
        // ambient sample, cell type
        uint8_t s_ix = n_sam;
        uint8_t ct_ix = n_ct;

        atac_pair_t atac_pair;
        atac_pair_init(&atac_pair);
        int rret = atac_prob_sample_read(sc_sim, ct_ix, 1, s_ix, &atac_pair);
        if (rret < 0) return -1;
        rret = atac_pair_set_name(&atac_pair, &sc_sim->atac_names);
        if (rret < 0) return -1;
        rret = atac_pair_set_seq(&atac_pair, bc_sim->atac_barcode);
        if (rret < 0) return -1;
        if (atac_pair_seq_error(sc_sim, &atac_pair) < 0) return -1;
        rret = atac_pair_write(&atac_pair, atac_fs);
        if (rret < 0) return -1;

        atac_pair_free(&atac_pair);
        ++sc_sim->atac_rd_counter;
    }

    // for ATAC ambient out reads
    for (r_i = 0; r_i < atac_ambn_rd_out; ++r_i) {
        // ambient sample, cell type
        uint8_t s_ix = n_sam;
        uint8_t ct_ix = n_ct;

        atac_pair_t atac_pair;
        atac_pair_init(&atac_pair);
        int rret = atac_prob_sample_read(sc_sim, ct_ix, 0, s_ix, &atac_pair);
        if (rret < 0) return -1;
        rret = atac_pair_set_name(&atac_pair, &sc_sim->atac_names);
        if (rret < 0) return -1;
        rret = atac_pair_set_seq(&atac_pair, bc_sim->atac_barcode);
        if (rret < 0) return -1;
        if (atac_pair_seq_error(sc_sim, &atac_pair) < 0) return -1;
        rret = atac_pair_write(&atac_pair, atac_fs);
        if (rret < 0) return -1;

        atac_pair_free(&atac_pair);
        ++sc_sim->atac_rd_counter;
    }

    return n_reads;
}

int sc_sim_gen_reads(sc_sim_t *sc_sim){
    if (sc_sim == NULL)
        return err_msg(-1, 0, "sc_sim_gen_reads: argument is null");

    // set streams
    if (sc_sim->out == NULL)
        return err_msg(-1, 0, "sc_sim_gen_reads: sc_sim->out must be "
                "set");

    if (sc_streams_set(&sc_sim->scs, sc_sim->out) < 0)
        return -1;

    // prevent overflow (variables stored in uint8_t)'
    assert(sc_sim->samples->n >= 0);
    assert(sc_sim->samples->n < 256);
    assert(sc_sim->K < 256);

    // Loop over each barcode
    size_t bc_i, bc_n = mv_size(&sc_sim->barcodes);
    int n_reads = 0, r = 0;
    for (bc_i = 0; bc_i < bc_n; ++bc_i) {
        int n;
        if ( (n = bc_sim_gen_reads(sc_sim, bc_i)) < 0)
            return -1;
        n_reads += n;
        int rm = n_reads / 1000000;
        if (rm > r) {
            log_msg("generated %i million reads", rm);
        }
        r = rm;
    }
    log_msg("generated %i million reads", n_reads);

    sc_streams_close(&sc_sim->scs);

    return 0;
}

/*
int sc_sim_sample_gex(sc_sim_t *sc_sim, int n_reads, int read_len) {
    if (sc_sim == NULL)
        return err_msg(-1, 0, "sc_sim_sample_gex: argument is null");

    if (test_sample_gex(sc_sim->gex_prob, sc_sim->fa, sc_sim->gex_prob->k, 
                n_reads, read_len, sc_sim->gv) < 0)
        return -1;

    return 0;
}

int sc_sim_sample_atac(sc_sim_t *sc_sim, int n_reads, int read_len) {
    if (sc_sim == NULL)
        return err_msg(-1, 0, "sc_sim_sample_atac: argument is null");

    if (test_sample_atac(sc_sim->atac_prob, sc_sim->gv, sc_sim->fa, n_reads, 
                read_len) < 0)
        return -1;

    return 0;
}
*/

/*******************************************************************************
 * print
 ******************************************************************************/

void bc_sim_print(FILE *fs, bc_sim_t *bc_sim){
    if (bc_sim == NULL) return;

    fprintf(fs, "%s\t", bc_sim->rna_barcode);
    size_t i;
    fprintf(fs, "samples: ");
    for (i = 0; i < mv_size(&bc_sim->samples); ++i){
        if (i) fprintf(fs, ",");
        fprintf(fs, "%u", mv_i(&bc_sim->samples, i));
    }
    fprintf(fs, "\tcell types: ");
    for (i = 0; i < mv_size(&bc_sim->cell_types); ++i){
        if (i) fprintf(fs, ",");
        fprintf(fs, "%u", mv_i(&bc_sim->cell_types, i));
    }

    fprintf(fs, "\t%u\t%u\t%u\t%u\t%u\t%u\n", 
            bc_sim->rna_nreads_cell, bc_sim->rna_nreads_ambn, 
            bc_sim->atac_nreads_cell, bc_sim->atac_nreads_ambn, 
            bc_sim->atac_npeak_cell, bc_sim->atac_npeak_ambn);
    
}

void sc_sim_print(FILE *fs, sc_sim_t *sc_sim){
    if (sc_sim == NULL) return;

    size_t i, n_bc = mv_size(&sc_sim->barcodes);
    for (i = 0; i < n_bc; ++i){
        bc_sim_t *bc = &mv_i(&sc_sim->barcodes, i);
        bc_sim_print(fs, bc);
    }
}


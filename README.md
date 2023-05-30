# Ambisim

Simulate a multiplexed single cell multiome RNA+ATAC experiment with
ambient molecule contamination.

## Installation

To download ambisim, you can clone from github
```bash
git clone https://github.com/marcalva/ambisim.git
```

To build ambisim, we first build a local htslib library
```bash
make hts
```

Then we can build ambisim
```bash
make
```

## Input files

`ambisim` requires several input files to run. Most importantly, a
droplet data file is required to specify how the droplets are
generated (see below).

A GTF file is required to obtain gene, transcript, and exon coordinates and
strands. It is recommended to put the `--tx-basic` flag in the arguments to
pull only transcripts and exons tagged as `basic` in the GTF file, if it is
present.  A fasta file is required to obtain the genome sequence. A cell type
ID file specifies the cell type IDs. An expression probability file gives a
matrix of gene expression probabilities per cell type, where the last column
specifies the ambient RNA profile. A similar file is required for open
chromatin peaks. A splice probability matrix is also included to specify the
probability of sampling spliced mature mRNA per cell type, with the last column
also giving the ambient probability. This matrix contains 1 row and `K+1`
columns. Finally, a genotype VCF file is required, containing the given sample
IDs.

### Droplet data file

The droplet file contains the following columns

- **RNA_BC** The RNA barcode sequence identifier (take from 10X whitelist if
  mapping with CellRanger). Any trailing characters such as '-1' should be
  removed.
- **ATAC_BC** The corresponding ATAC barcode sequence identifier 
  (take from 10X whitelist if mapping with CellRanger). Any trailing
  characters such as '-1' should be removed.
- **Type** The droplet type, one of empty, singlet, or doublet.
- **n** Number of nuclei in droplet, one of 0, 1, or 2.
- **sam** The sample ID(s) corresponding the one of the VCF sample IDs.
  This value is ignored if n=0. If the droplet is a doublet (n=2), set the 
  join two sample IDS with a comma delimiter, e.g. 'S1, S2'.
- **ct** The numeric index of the cell type of the droplet (0, 1, ... K-1).
  If empty (n=0), this value is ignored. If doublet, join the cell type IDs
  with a comma, e.g. '0,3'.
- **rna_nr_c** The number of gene expression UMIs from the nucleus.
- **rna_nr_a** The number of gene expression UMIs from the ambient pool.
- **atac_nr_c** The number of ATAC DNA reads from the nucleus.
- **atac_nr_a** The number of ATAC DNA reads from the ambient pool
- **atac_nr_cp** The number of ATAC DNA reads from the nucleus/cell inside
  a peak.
- **atac_nr_ap** The number of ATAC DNA reads from the ambient pool inside
  a peak.

## Utility scripts

Utility R scripts are provided to help produce the droplet data input file
required for ambisim.

If you have a feature barcode matrix from 10X and a set of barcodes to cluster,
you can run `utils/probs_from_data.R` to obtain cell types, their proportions,
gene and peak probabilities. You can run the script as

```bash
./utils/probs_from_data.R \
    path/to/counts \
    flt_barcodes.txt \
    out/
```

where `path/to/counts` contains `features.tsv.gz`, `barcodes.tsv.gz`, and
`matrix.mtx.gz` from a 10X multiome run. The `flt_barcodes.txt` file contains
the barcodes to cluster the expression data on, and `out/` is the directory to
place the output files. The file `ct_freq.txt` contains the cell type
frequencies, `cell_types` contains the cell type IDs, `gene_probs.txt` contains
the gene expression probabilities per cell type, `gene_ids.txt` contains the
gene IDs, `peak_probs.txt` contains the accessibility probabilities per peak,
and `spl_probs.txt` contains splice probabilities per cell type (fixed at 0.4
for cell types and 0.6 for ambient).

If you want to obtain the the `n` expressed barcodes for clustering, you can
run, for example, 

```bash
./utils/get_top_n_bc.R \
    path/to/counts \
    3000 \
    top_3k_bcs.txt
```

where `path/to/counts` contains the raw counts output from a 10X multiome run
(same as above), `3000` specifies to output the top 3,000 barcodes ranked by
number of UMIs, and `top_3k_bcs.txt` is the output file listing one barcode
per line.

Proportions of samples must be provided to the droplet data script. If you want
to even proportions for all samples, you can run

```bash
./utils/mk_sam_props.R \
    sample_ids.txt \
    sample_freq.txt
```

where `sample_ids.txt` contains a sample ID per line, and `sample_freq.txt` contains
the output of uniform proportions per sample.

To generate the droplet file, you can run

```bash
./utils/mk_sam_props.R \
    rna_wl_bcs \
    atac_wl_bcs \
    sample_freq.txt \
    ct_freq.txt \
    1 \
    drop_data_rand.txt
```

where `rna_wl_bcs` is the file containing all RNA whitelist barcodes,
`atac_wl_bcs` contains all ATAC whitelist barcodes (each RNA and ATAC barcode
line correspond to matching barcodes). These barcodes can be found in a
10X cellranger arc download in

```bash
lib/python/cellranger/barcodes/737K-arc-v1.txt.gz
lib/python/atac/barcodes/737K-arc-v1.txt.gz
```

The `sample_freq.txt` file contains the desired sample ID proportions,
`sample_freq.txt` contains the desired cell type proportions, `1` specifies the
random seed number, and `drop_data_rand.txt` specifies the output file.

Note that the ambient contamination is sampled randomly and uniformly for
nuclei clusters. The parameters for the random distributions can be
changed, as well as any parts of the script.

Finally, ambisim can be run as
```bash
ambisim \
    --gtf gencode.gtf \
    --fasta fasta.fa \
    --samples sample_ids.txt \
    --cell-types cell_types.txt \
    --drop-file drop_data_rand.txt \
    --expr-prob gene_probs.txt \
    --gene-ids gene_ids.txt \
    --mmrna-prob spl_probs.txt \
    --peak-prob peak_probs.txt \
    --peaks atac_peaks.bed \
    --vcf genotypes.vcf \
    --seed 1 \
    --tx-basic \
    --out ambisim_out
```

# Options

```
Input files

  -v, --vcf                Indexed VCF file.
  -s, --samples            File containing sample IDs in VCF file to simulate genotypes.
                           If this option is not given, all samples in VCF will be used.
  -c, --cell-types         File containing cell type IDs. One ID per line.
  -g, --gtf                Path to GTF file. Required to obtain gene, transcript, and
                           exons coordinates.
  -f, --fasta              FASTA file for genome, required for pulling the genome sequence.
                           Note that the chromosome names must match with the VCF and GTF files.
  -d, --drop-file          Path to droplet data table (see README for explanation).
  -p, --expr-prob          File path to matrix with gene expression probabilities. Each column gives
                           expression probabilities for each gene. The last column is ambient.
                           No row names or header.
  -G, --gene-ids           Path to file containing gene IDs. Each row contains the gene ID of the
                           corresponding row in the 'expr-prob' file. Must have the same number of lines.
  -m, --mmrna-prob         Path to file with splicing probabilities per cell type. Gives a numeric matrix
                           with one row and K + 1 columns, where columns represent cell types.
                           Last column is ambient. No row names or header.
  -r, --peak-prob          Path to file with numeric matrix giving open chromatin peak accessibility probabilities.
                           Rows represent peaks and columns represent cell types. Last column is ambient.
                           No row names or header.
  -o, --peaks              Path to BED file containing peaks. Number of peaks must match number of rows
                           in peak probabilities file.
  -O, --out                File path to output folder [sim].
                           Subdirectories 'RNA' and 'ATAC' will be created, where fastq files will be output.

Sequencing options

  -e, --seq-error          Probability of a sequencing error [0.01].
  -l, --rna-rd-len         Length of RNA read [91].
  -U, --rna-umi-len        Length of RNA UMI [12].
  -L, --atac-rd-len        Length of each ATAC read pair [50].

  -S, --seed               Seed for random number generation.


GTF options

  -t, --tx-basic           Read only transcripts tagged as 'basic' in the GTF file.

  -V, --verbose            Write status on output.
      --help               Print this help screen.

```

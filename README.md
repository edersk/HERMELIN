# HERMELIN

HERMELIN genotypes single nucleotide variants (SNVs) from Oxford Nanopore
long-read single-cell RNA sequencing data. It streams FASTQ files to extract
cell barcodes and UMIs, calls the base at user-defined regions of interest,
performs UMI error correction, and classifies cells as mutant or wild-type with
optional SNP-based phasing.

## Installation

```r
# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "IRanges", "S4Vectors", "ShortRead",
                       "BiocGenerics", "pwalign"))

# Install HERMELIN from GitHub
# install.packages("remotes")
remotes::install_github("edersk/HERMELIN")
```

HERMELIN also depends on the CRAN packages `data.table`, `stringr`, `stringdist`
and `inflection`, which are installed automatically.

## Pipeline

Three functions are run in sequence:

```
FASTQ --> detect_snv() --> summarize_snv() --> flag_snv()
```

1. **`detect_snv()`** — streams the FASTQ, detects the Read 1 adapter, extracts
   the cell barcode (CBC) and UMI, requires a poly(T) sequence, and calls the
   base at the major region of interest (ROI) and an optional minor ROI (a SNP
   for phasing).
2. **`summarize_snv()`** — corrects UMIs by Levenshtein distance, aggregates
   reads per cell, detects knee/inflection thresholds, and reports SNP allele
   frequencies to gauge phasing feasibility.
3. **`flag_snv()`** — classifies each cell as `MUT`, `WT`, or `MUT_WT`, with an
   optional SNP-phased call that separates definitive from ambiguous wild-type.

## Quick run

The minimal set of arguments needed to go from FASTQ to per-cell genotype
calls. Everything not shown here has a sensible default (see
[Function reference](#function-reference) below).

```r
library(HERMELIN)

# Step 1: detect SNVs from FASTQ
detect_snv(
    session_name       = "experiment_name",
    path_input_file    = "path/to/experiment.fastq.gz",
    path_output_folder = "output/",
    read1_seq          = "CTACACGACGCTCTTCCGATCT",  # Read 1 adapter
    read1_mismatch     = 3L,
    major_roi_seq      = "GAGATTTCNCTGTAGCT",       # 'N' marks the SNV position
    major_roi_mismatch = 3L,
    major_roi_mut      = "T",
    major_roi_wt       = "A",
    minor_roi_seq      = "AGAACAATNCCAAATGC",       # optional SNP for phasing
    minor_roi_mismatch = 3L,
    minor_roi_mut      = "C",
    minor_roi_wt       = "T",
    n_cores            = 4L
)

# Step 2: correct UMIs and summarize per cell
res <- summarize_snv(
    session_name       = "experiment_name",
    path_snv_table     = "output/experiment_name_snv_table.csv",
    path_barcodes      = "path/to/barcodes.csv",   # cell-barcode whitelist
    path_output_folder = "output/",
    n_cores            = 4L
)

res$knee_points        # suggested reads-per-UMI thresholds
res$phasing_message    # whether phasing looks viable for this dataset

# Step 3: classify cells (optionally phased)
flags <- flag_snv(
    path_snv_summary   = "output/experiment_name_snv_summary.csv",
    path_output_folder = "output/",
    session_name       = "experiment_name",
    threshold          = 10,       # pick from res$knee_points
    use_phasing        = TRUE,
    minor_on_major     = "MUT"     # minor SNP status expected on the mutant allele
)

table(flags$mutation_call)
table(flags$phased_call)
```

If you do not have a phasing SNP, omit the four `minor_roi_*` arguments in
step 1 and leave `use_phasing = FALSE` in step 3.

## Input requirements

- **FASTQ** (`.fastq` or `.fastq.gz`) of Oxford Nanopore long reads.
- **ROI sequence** — the region flanking the variant, with a single `N` at the
  SNV position, e.g. `GAGATTTCNCTGTAGCT`. Other IUPAC codes are allowed in
  the flanks.
- **Barcode whitelist** — a `.csv` (or plain-text) list of expected cell
  barcodes, one per line; only the first 16 bp are used for matching, so 10x
  suffixes such as `-1` are tolerated.
- A **minor ROI** is optional; provide one only if you want SNP-based phasing.

## Outputs

| File                                 | Written by      | Contents                                                                                                                                |
| ------------------------------------ | --------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `<session>_snv_table.csv`            | `detect_snv`    | One row per read: `read_id`, `CBC`, `UMI`, `strand`, `major_base`/`major_status`, `minor_base`/`minor_status` (+22 columns if enriched) |
| `<session>_log.csv`                  | `detect_snv`    | Run parameters and read-filtering counts                                                                                                |
| `<session>_snv_summary.csv`          | `summarize_snv` | One row per cell, `;`-separated: `cbc;umis;n_umi;n_umi_mut;n_umi_wt;n_snp;n_snp_mut;n_snp_wt` (count columns are `,`-separated per UMI) |
| `<session>_umi_frequencies.csv`      | `summarize_snv` | Ranked reads-per-UMI, for choosing a threshold                                                                                          |
| `<session>_thresholds.csv`           | `summarize_snv` | All knee-point estimates and the suggested range                                                                                        |
| `<session>_*_roi_frequency.csv`      | `summarize_snv` | Major/minor allele frequency tables (plus `_snp_frequency_combined.csv`)                                                                |
| `<session>_umi_distribution.pdf`     | `summarize_snv` | Diagnostic plots (rank-abundance, histogram, cumulative, violins)                                                                       |
| `<session>_flags.csv`                | `flag_snv`      | Per-cell genotype: `cbc`, `barcode_10x`, `has_np_data`, `mutation_call`, `major_call_simple`, `phased_call`, `phased_call_simple`       |

## Genotype categories

**`mutation_call`** (major ROI): `MUT`, `WT`, `MUT_WT` (heterozygous — both
alleles), `below_threshold` (no UMI cleared the read threshold), or `no_np_data`
(barcode present in the whitelist but had no Nanopore reads). `major_call_simple`
collapses these to `MUT` / `WT`.

**`phased_call`** (only meaningful with `use_phasing = TRUE`): the minor SNP is
used to qualify the major call.

| `phased_call`          | `phased_call_simple` | Meaning                                                  |
| ---------------------- | -------------------- | -------------------------------------------------------- |
| `MUT_mnsnp`            | `MUT`                | Mutant, confirmed by the expected minor SNP              |
| `MUT_non_mnsnp`        | `MUT`                | Mutant, minor SNP not observed                           |
| `MUT_both`             | `MUT`                | Mutant with a mix of observed and non-observed minor SNP |
| `def_WT`               | `defWT`              | Definitive wild-type (phasing-consistent)                |
| `amb_WT`               | `ambWT`              | Ambiguous wild-type (could not be phased)                |
| `MUT_ERROR` / `ERROR1` | —                    | Minor SNP contradicts the major call                     |

## Choosing a threshold

`summarize_snv()` returns `res$knee_points` and writes
`<session>_umi_distribution.pdf`. Pick the `flag_snv()` `threshold` within the
suggested range: a **higher** threshold is more specific (fewer, higher-
confidence calls); a **lower** threshold is more sensitive with possible
false-positive calls.

## Function reference

All arguments of the three pipeline functions. "Required" means there is no
default and the argument must be supplied.

### `detect_snv()`

Streams the FASTQ in chunks. For every read it searches for the Read 1 adapter
on both strands, extracts the CBC and UMI immediately downstream, checks for a
poly(N) tail, and then locates the major (and minor) ROI with a two-pass
strategy: an exact IUPAC match first, then a fuzzy match followed by pairwise
alignment for reads that did not match exactly.

#### Session and I/O

| Argument             | Default  | Description                                                                                                                                 |
| -------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `session_name`       | required | Prefix for all output files. Letters, digits, `_` and `-` only.                                                                             |
| `path_input_file`    | required | Path to the FASTQ file (`.fastq` or `.fastq.gz`).                                                                                           |
| `path_output_folder` | required | Existing directory where `<session>_snv_table.csv` and `<session>_log.csv` are written.                                                     |
| `n_cores`            | required | CPU cores for the pairwise-alignment pass. Automatically capped at `parallel::detectCores() - 1`.                                           |
| `n_reads_per_chunk`  | `1e5`    | Reads per streaming chunk. Larger chunks are faster but use more memory.                                                                    |
| `enriched`           | `FALSE`  | If `TRUE`, appends 22 diagnostic columns (read length, adapter/ROI match widths, substitution/insertion/deletion counts, Q-scores, match method). Use for QC and publication analyses; not needed for genotyping. |

#### Read 1 adapter and barcode layout

| Argument             | Default  | Description                                                                                                                                 |
| -------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `read1_seq`          | required | Adapter sequence that immediately precedes the cell barcode (A/C/G/T only), e.g. `CTACACGACGCTCTTCCGATCT` for 10x 3' libraries.             |
| `read1_mismatch`     | required | Maximum edit distance (substitutions + indels) allowed when matching the adapter. Reads with no match, or with a match on both strands (chimeric), are discarded. |
| `read1_within_range` | `70`     | Number of bases from the read start (on each strand) in which the adapter is searched.                                                      |
| `len_CBC`            | `16`     | Cell barcode length, taken directly after the adapter match.                                                                                |
| `len_UMI`            | `12`     | UMI length, taken directly after the barcode.                                                                                               |

#### Poly(N) tail filter

| Argument         | Default | Description                                                                                                                                           |
| ---------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `polyN_base`     | `"T"`   | Base of the homopolymer expected after the UMI (poly(T) for 3' GEX). Set to `NULL` to disable the filter, e.g. for 5' GEX libraries where the TSO follows the UMI. |
| `polyN_length`   | `10`    | Window length checked for the homopolymer. Ignored when `polyN_base = NULL`.                                                                          |
| `polyN_mismatch` | `2`     | Number of non-`polyN_base` positions tolerated within that window. Ignored when `polyN_base = NULL`.                                                  |

#### Major ROI (the SNV to genotype)

| Argument             | Default  | Description                                                                                                                                 |
| -------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `major_roi_seq`      | required | Flanking sequence with exactly one `N` at the SNV position. IUPAC codes (R, Y, S, W, K, M, B, D, H, V) are allowed elsewhere.                |
| `major_roi_mismatch` | required | Maximum edit distance for the fuzzy ROI search (pass 2). Reads that match more than once are left uncalled.                                 |
| `major_roi_mut`      | `NULL`   | Base that labels a read `MUT` in `major_status`. Must be given together with `major_roi_wt`.                                                |
| `major_roi_wt`       | `NULL`   | Base that labels a read `WT` in `major_status`. If both are `NULL`, only the raw `major_base` is reported and `summarize_snv()` cannot be run. |

#### Minor ROI (optional phasing SNP)

| Argument             | Default | Description                                                                                                                           |
| -------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `minor_roi_seq`      | `NULL`  | Flanking sequence of a heterozygous SNP on the same transcript, with one `N` at the SNP position. `NULL` skips minor ROI detection entirely. |
| `minor_roi_mismatch` | `NULL`  | Maximum edit distance for the minor ROI search. Required when `minor_roi_seq` is given.                                              |
| `minor_roi_mut`      | `NULL`  | Base that labels a read `MUT` in `minor_status`. Must be given together with `minor_roi_wt`.                                          |
| `minor_roi_wt`       | `NULL`  | Base that labels a read `WT` in `minor_status`.                                                                                       |

**Returns** `NULL` invisibly; results are written to disk.

### `summarize_snv()`

Loads the SNV table, keeps reads with a `MUT`/`WT` major status and a
whitelisted barcode, collapses UMIs within each cell that are within
`umi_mismatch` edits of a more abundant UMI, counts `MUT`/`WT` reads per
corrected UMI, and estimates reads-per-UMI thresholds with several knee-point
methods (geometric knee, inflection, UIK, CellRanger-style, Otsu, Gaussian
mixture).

| Argument             | Default  | Description                                                                                                                                 |
| -------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `session_name`       | required | Prefix for output files.                                                                                                                    |
| `path_snv_table`     | required | Path to `<session>_snv_table.csv` written by `detect_snv()`. Only the columns needed for summarization are read, so enriched tables are fine. |
| `path_barcodes`      | required | Cell-barcode whitelist (`.csv` read as the first column, otherwise one barcode per line). Reads whose first 16 bp of CBC are not in this list are dropped. |
| `path_output_folder` | required | Existing directory for the summary, frequency tables, thresholds and diagnostic PDF.                                                        |
| `umi_mismatch`       | `2L`     | Maximum Levenshtein distance for UMI collapsing. A UMI within this distance of a more frequent UMI in the same cell is merged into it.        |
| `n_cores`            | `4L`     | CPU cores. Capped at `parallel::detectCores() - 1`.                                                                                        |
| `skip_plots`         | `FALSE`  | If `TRUE`, do not write `<session>_umi_distribution.pdf`.                                                                                   |

**Returns** (invisibly) a list with:

| Element           | Contents                                                                                                   |
| ----------------- | ---------------------------------------------------------------------------------------------------------- |
| `summary`         | Per-cell `data.table` (same content as `<session>_snv_summary.csv`); can be passed directly to `flag_snv()` |
| `umi_frequencies` | Reads per corrected UMI, ranked                                                                             |
| `knee_points`     | List of threshold estimates: `upper_knee`, `inflection`, `lower_knee`, `uik`, `cellranger`, `otsu`, `mixture`, plus the suggested `lower`/`upper` range |
| `snp_frequencies` | Major/minor allele frequency tables                                                                        |
| `phasing_viable`  | `TRUE` when both minor alleles are each ≥ 10 % of UMIs                                                     |
| `phasing_message` | Human-readable phasing assessment                                                                          |
| `has_minor_roi`   | Whether minor ROI data was found in the input                                                              |
| `statistics`      | Read and cell counts at each filtering stage                                                               |

### `flag_snv()`

Assigns a genotype to every cell. Per cell, each UMI is kept if its total
reads (`MUT` + `WT`) reach `threshold`, called `MUT` or `WT` if at least
`consensus_threshold` of its reads agree, and the cell call is derived from the
set of called UMIs. With phasing enabled, each called UMI is additionally
qualified by its minor SNP status.

| Argument              | Default       | Description                                                                                                                                                  |
| --------------------- | ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `path_snv_summary`    | required      | Path to `<session>_snv_summary.csv`, or the `summary` `data.table` returned by `summarize_snv()`.                                                            |
| `path_barcodes`       | `NULL`        | Optional full barcode whitelist used for **backfilling**: barcodes in this list with no Nanopore reads are added to the output with `mutation_call = "no_np_data"`. |
| `path_barcode_filter` | `NULL`        | Optional barcode list used for **filtering**: only barcodes in this list are analysed and reported.                                                           |
| `path_output_folder`  | `NULL`        | Directory for `<session>_flags.csv`. If `NULL`, nothing is written and the result is only returned.                                                          |
| `session_name`        | `"snv_flags"` | Prefix for the output file.                                                                                                                                  |
| `threshold`           | `10`          | Minimum total reads (`MUT` + `WT`) a UMI needs to be considered. Choose from the range reported by `summarize_snv()`.                                        |
| `consensus_threshold` | `0.8`         | Fraction of a UMI's reads that must support one allele for the UMI to be called `MUT` or `WT`. UMIs below this are ignored.                                  |
| `use_phasing`         | `FALSE`       | Use the minor SNP columns to produce `phased_call`. Falls back to unphased calls (with a warning) if the summary has no SNP columns.                          |
| `minor_percentage`    | `0.2`         | A UMI's minor SNP status is only used if its minor reads are at least this fraction of its major reads; otherwise the minor status is `none`.                 |
| `minor_on_major`      | `"MUT"`       | Which minor SNP allele (`"MUT"` or `"WT"`) is expected on the same transcript copy as the major mutation. Determines whether an observed minor SNP confirms (`MUT_mnsnp`, `def_WT`) or contradicts (`MUT_ERROR`) the major call. |

**Returns** a `data.table` with one row per cell: `cbc`, `barcode_10x`,
`has_np_data`, `mutation_call`, `major_call_simple`, `phased_call`,
`phased_call_simple` (see [Genotype categories](#genotype-categories)).

## License

MIT © 2026 Sebastian Eder

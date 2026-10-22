# HERMELIN

**H**aplotype-resolved **E**nrichment of **R**NA for **M**utation and
**E**xpression **L**inkage **IN** single cells.

### Fast, customizable single-cell mutation genotyping.

HERMELIN is an easy-to-use, user-customizable open-source R package for
per-cell genotyping of low-abundance single-nucleotide variants (SNVs).
It directly links driver mutations to cellular phenotypes through
single-cell gene expression profiles, helping researchers compare mutant
and wild-type cell states.

-   **Fast analysis:** Go from Nanopore FASTQ files to per-cell mutation
    calls in three steps, without prior genome alignment.
-   **Cost-effective workflow:** Use surplus full-length cDNA from
    standard scRNA-seq libraries for targeted Oxford Nanopore
    sequencing.
-   **Flexible settings:** Define your target variants, sequence motifs,
    and filtering thresholds.
-   **Haplotype-aware comparisons:** When linked SNPs allow phasing,
    distinguish definite wild-type (`defWT`) from ambiguous wild-type
    (`ambWT`) calls.

> [!NOTE]  
> **Regular updates and user support**  
>
> We will continue to provide regular updates and support for HERMELIN.
> Questions, feedback, bug reports, and feature requests are welcome via
> <a href="https://github.com/edersk/HERMELIN/issues">
> <img src="https://cdn.simpleicons.org/github" alt="GitHub" width="18" height="18"/>
> GitHub Issues </a>.
>  

## Contents

-   [Installation](#installation)
-   [Inputs and library
    compatibility](#inputs-and-library-compatibility)
-   [Workflow](#workflow)
-   [Quick start without phasing](#quick-start-without-phasing)
-   [Add haplotype phasing](#add-haplotype-phasing)
-   [Interpret genotype calls](#interpret-genotype-calls)
-   [Choose thresholds](#choose-thresholds)
-   [Add calls to a Seurat object](#add-calls-to-a-seurat-object)
-   [Output files](#output-files)
-   [Function reference](#function-reference)
-   [License](#license)

## Installation {#installation}

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

BiocManager::install(c(
    "Biostrings", "IRanges", "S4Vectors", "ShortRead",
    "BiocGenerics", "pwalign"
))

remotes::install_github("edersk/HERMELIN")
```

The CRAN dependencies `data.table`, `stringr`, `stringdist`, and
`inflection` are installed automatically. Use an R version supported by
your selected [Bioconductor release](https://bioconductor.org/install/).

## Inputs and library compatibility {#inputs-and-library-compatibility}

| Input                  | What to provide                                                                                                                        |
|------------------------------------|------------------------------------|
| FASTQ                  | A `.fastq` or `.fastq.gz` file of Oxford Nanopore reads retaining the cell barcode, UMI, and target region.                            |
| Target motif           | The sequence surrounding the SNV, with exactly one `N` at the variant position. Other IUPAC ambiguity codes are allowed in the flanks. |
| Target alleles         | The mutant and wild-type bases in the orientation of the supplied motif.                                                               |
| Barcode whitelist      | A headerless `.csv` or plain-text file with one expected cell barcode per line. For CSV input, the first column is used.               |
| Optional phasing motif | A linked heterozygous SNP on the same transcript, with a known phase relationship to the target mutation.                              |

Use a sample-specific list of retained scRNA-seq cell barcodes when you
want the output to represent the cells in that experiment. A
chemistry-wide whitelist contains possible barcodes, rather than the
cells retained in your transcriptomic analysis.

The default layout is a Read 1 adapter followed by a 16-base cell
barcode and a 12-base UMI, with a downstream poly(T) sequence. Check the
layout for your library chemistry. For 5′ libraries where a
template-switch oligonucleotide follows the UMI, set `polyN_base = NULL`
and verify the adapter and UMI settings.

Current downstream whitelist matching uses the first 16 bases,
tolerating suffixes such as `-1`. Although `detect_snv()` exposes
`len_CBC`, changing that argument alone does not adapt downstream
matching to other barcode lengths. Barcodes with sample prefixes need an
explicit mapping to their original barcode sequence. Process samples
separately before combining their metadata.

## Workflow {#workflow}

Run three functions in sequence:

``` text
FASTQ -> detect_snv() -> summarize_snv() -> flag_snv() -> per-cell calls
```

1.  **`detect_snv()`** streams the FASTQ, detects the Read 1 adapter,
    extracts cell barcodes and UMIs, checks the configured homopolymer,
    and calls bases at the target motif and an optional phasing motif.
2.  **`summarize_snv()`** filters to the barcode whitelist, corrects
    UMIs within cells, summarizes read counts, estimates reads-per-UMI
    thresholds, and reports allele frequencies for an initial assessment
    of phasing feasibility.
3.  **`flag_snv()`** applies read-depth and consensus thresholds to
    corrected UMIs, then assigns cell-level calls. Optional phasing
    qualifies wild-type calls as definite or ambiguous.

The target motif is called the **major ROI** in the function arguments.
The optional phasing motif is called the **minor ROI**. Here, major and
minor describe their roles, not their allele frequencies.

HERMELIN uses sequence motifs to locate a target single-nucleotide
variant (SNV), so the computational workflow does not require prior
alignment to a reference genome. Reads must retain the adapter, cell
barcode, UMI, and target sequence needed for the analysis. Target
enrichment and long-read sequencing are performed before running this R
package.

## Quick start without phasing {#quick-start-without-phasing}

This example shows the required files and the order of operations.
Replace the paths and example motif with those for your experiment. The
motif and allele values below are illustrative and are not a universal
variant assay. Allele definitions must match the orientation of the
motif, which may differ from the genomic reference strand.

``` r
library(HERMELIN)

session_name <- "experiment_name"
fastq_path <- "path/to/experiment.fastq.gz"
barcode_path <- "path/to/barcodes.csv"
output_dir <- "output"
n_cores <- if (.Platform$OS.type == "windows") 1L else 4L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(fastq_path), file.exists(barcode_path))

# 1. Extract barcode/UMI information and call the target base.
detect_snv(
    session_name       = session_name,
    path_input_file    = fastq_path,
    path_output_folder = output_dir,
    read1_seq          = "CTACACGACGCTCTTCCGATCT",
    read1_mismatch     = 3L,
    major_roi_seq      = "GAGATTTCNCTGTAGCT",  # N marks the target position
    major_roi_mismatch = 3L,
    major_roi_mut      = "T",                  # example allele definitions
    major_roi_wt       = "A",
    n_cores            = n_cores
)

# 2. Filter to expected cells, correct UMIs, and inspect coverage.
res <- summarize_snv(
    session_name       = session_name,
    path_snv_table     = file.path(
        output_dir, paste0(session_name, "_snv_table.csv")
    ),
    path_barcodes      = barcode_path,
    path_output_folder = output_dir,
    n_cores            = n_cores
)

res$knee_points
# Inspect output/experiment_name_umi_distribution.pdf before proceeding.
# Here the inflection estimate is used as an explicit starting choice.
umi_threshold <- res$knee_points$inflection
if (length(umi_threshold) != 1L ||
    !is.finite(umi_threshold) || umi_threshold < 1) {
    stop("No usable inflection estimate. Inspect the plots and choose a threshold.")
}

# 3. Assign cell-level calls and retain the complete expected cell list.
flags <- flag_snv(
    path_snv_summary   = res$summary,
    path_barcodes      = barcode_path,
    path_output_folder = output_dir,
    session_name       = session_name,
    threshold          = umi_threshold,
    use_phasing        = FALSE
)

table(flags$mutation_call, useNA = "ifany")
table(flags$major_call_simple, useNA = "ifany")
```

The example uses three allowed edits for adapter and motif matching.
These values are explicit settings to review for your data. The depth
estimate is also a choice to inspect, not a universal cutoff.

Supplying `path_barcodes` to `flag_snv()` includes expected cells that
are absent from the filtered summary. These receive `no_np_data`.
Without that argument, the output covers only barcodes present in the
summary. Keep this distinction in mind when reporting the fraction of
cells genotyped.

## Add haplotype phasing {#add-haplotype-phasing}

For a heterozygous target mutation, phasing requires a linked
heterozygous SNP and enough reads covering both sites on the same
molecules. Establish which SNP allele is linked to the target-mutant
haplotype in your sample before setting `minor_on_major`.

The SNP's `minor_roi_mut` and `minor_roi_wt` labels simply name its two
alleles. They do not imply that the SNP itself is pathogenic. If the
allele labeled `MUT` is linked to the target mutation, use
`minor_on_major = "MUT"`. If the allele labeled `WT` is linked to it,
use `minor_on_major = "WT"`.

Rerun detection and summarization with the minor motif included. A
separate session name keeps the phased and unphased output files
distinct.

``` r
phased_session <- paste0(session_name, "_phased")

detect_snv(
    session_name       = phased_session,
    path_input_file    = fastq_path,
    path_output_folder = output_dir,
    read1_seq          = "CTACACGACGCTCTTCCGATCT",
    read1_mismatch     = 3L,
    major_roi_seq      = "GAGATTTCNCTGTAGCT",
    major_roi_mismatch = 3L,
    major_roi_mut      = "T",
    major_roi_wt       = "A",
    minor_roi_seq      = "AGAACAATNCCAAATGC",  # replace with your linked SNP
    minor_roi_mismatch = 3L,
    minor_roi_mut      = "C",                 # label C as minor MUT
    minor_roi_wt       = "T",                 # label T as minor WT
    n_cores            = n_cores
)

res_phased <- summarize_snv(
    session_name       = phased_session,
    path_snv_table     = file.path(
        output_dir, paste0(phased_session, "_snv_table.csv")
    ),
    path_barcodes      = barcode_path,
    path_output_folder = output_dir,
    n_cores            = n_cores
)

res_phased$phasing_message
if (!isTRUE(res_phased$has_minor_roi)) {
    stop("No usable minor-ROI calls. Check the SNP motif and coverage.")
}
if (!isTRUE(res_phased$phasing_viable)) {
    warning("Review the minor-allele frequencies before using phased calls.")
}

phased_threshold <- res_phased$knee_points$inflection
if (length(phased_threshold) != 1L ||
    !is.finite(phased_threshold) || phased_threshold < 1) {
    stop("Inspect the phased-session coverage plots and choose a threshold.")
}

flags_phased <- flag_snv(
    path_snv_summary   = res_phased$summary,
    path_barcodes      = barcode_path,
    path_output_folder = output_dir,
    session_name       = phased_session,
    threshold          = phased_threshold,
    use_phasing        = TRUE,
    minor_on_major     = "MUT"  # only if C is known to mark the mutant haplotype
)

table(flags_phased$phased_call, useNA = "ifany")
table(flags_phased$phased_call_simple, useNA = "ifany")
```

`phasing_viable` is a dataset-level heuristic based on the
representation of both SNP alleles, each accounting for at least 10% of
the minor-ROI UMI counts. It does not establish the phase relationship
or guarantee that every cell can be phased.

For the example allele labels and `minor_on_major = "MUT"`, the key
logic is:

| Target observation | Linked SNP observation      | Interpretation                                                                                                   |
|------------------------|------------------------|------------------------|
| Target `MUT`       | SNP `MUT`                   | Supports a phase-consistent mutant call.                                                                         |
| Target `WT`        | SNP `MUT`                   | Supports `defWT`: the marker haplotype expected to carry the mutation instead carries the wild-type target base. |
| Target `WT`        | SNP `WT`, or unresolved SNP | Remains `ambWT`: the mutant-linked haplotype was not resolved as wild-type.                                      |

These rows describe molecule-level evidence. Final cell calls combine
the retained molecules and account for mutant evidence and phase
conflicts.

## Interpret genotype calls {#interpret-genotype-calls}

For heterozygous mutations, an unphased `WT` call means that wild-type
transcripts were observed. It does not establish that the cell lacks the
mutation, because the mutant allele may not have been captured. Phasing
helps qualify these observations but cannot recover transcripts missing
from the original cDNA library.

### Unphased target calls

| `mutation_call`   | `major_call_simple` | Meaning                                                                                                                                               |
|------------------------|------------------------|------------------------|
| `MUT`             | `MUT`               | Retained called UMIs support the mutant allele.                                                                                                       |
| `WT`              | `WT`                | Retained called UMIs support the wild-type allele. This does not exclude a heterozygous mutation without phasing.                                     |
| `MUT_WT`          | `MUT`               | Both mutant and wild-type transcript molecules were detected. This is an RNA observation, not independent confirmation of genomic zygosity.           |
| `below_threshold` | `NA`                | No UMI passed both the read-depth and allele-consensus requirements.                                                                                  |
| `no_np_data`      | `NA`                | An expected barcode is absent from the filtered summary. This can reflect absent reads or earlier filtering, not necessarily zero raw Nanopore reads. |

`has_np_data` indicates whether the barcode has a summary entry with
target counts. It does not report raw FASTQ coverage independently of
filtering. Uncalled states stay `NA` in the simplified columns and must
not be treated as wild-type.

### Phased calls

Interpret these columns as phased evidence only when phasing was enabled
and usable SNP data were supplied. The output can contain fallback
labels when phasing is disabled, so the presence of `phased_call` alone
does not establish that phasing was performed.

| `phased_call`                    | `phased_call_simple` | Meaning                                                                                                         |
|------------------------|------------------------|------------------------|
| `MUT_mnsnp`                      | `MUT`                | Mutant evidence supported by the expected SNP allele.                                                           |
| `MUT_non_mnsnp`                  | `MUT`                | Mutant evidence without resolved SNP support.                                                                   |
| `MUT_both`                       | `MUT`                | Both SNP-supported and SNP-unresolved mutant UMIs were observed.                                                |
| `def_WT`                         | `defWT`              | Wild-type target evidence on the marker haplotype expected to carry the mutation, with no retained mutant call. |
| `amb_WT`                         | `ambWT`              | Wild-type target evidence without enough phase information to assign `defWT`.                                   |
| `MUT_ERROR` / `ERROR1`           | `NA`                 | Target and marker evidence conflict with the configured phase relationship. Review these cells.                 |
| `below_threshold` / `no_np_data` | `NA`                 | No qualifying call.                                                                                             |

Use `defWT` as the phase-supported wild-type comparison group when
available. Keep `ambWT` separate rather than treating it as a confirmed
wild-type group.

## Choose thresholds {#choose-thresholds}

These settings apply to different stages of the analysis:

| Setting                                                      | Applies to                    | Meaning                                                                                                                                |
|------------------------|------------------------|------------------------|
| `threshold`                                                  | Each corrected UMI            | Minimum number of target-called reads (`MUT` + `WT`) required, inclusive.                                                              |
| `consensus_threshold`                                        | Allele calling within a UMI   | Minimum fraction of reads supporting one allele, default `0.8`. The same consensus requirement is used for the minor SNP when phasing. |
| `minor_percentage`                                           | Phasing coverage within a UMI | Minimum minor-SNP read count relative to target-called read count, default `0.2`. It is not a cell-level VAF cutoff.                   |
| `read1_mismatch`, `major_roi_mismatch`, `minor_roi_mismatch` | Read matching                 | Allowed edit distances for sequence detection.                                                                                         |
| `umi_mismatch`                                               | UMI correction within a cell  | Maximum edit distance for collapsing a less abundant UMI into a more abundant one, default `2L`.                                       |

`summarize_snv()` returns estimates in `res$knee_points` and writes
`<session>_umi_distribution.pdf`. Inspect the ranked reads-per-UMI curve
and the estimates before selecting a value for `flag_snv()`.

A higher depth threshold retains fewer, more deeply supported molecules.
A lower threshold retains more molecules, including potentially
erroneous ones. Neither choice alone establishes accuracy for a new
dataset. The inflection estimate is used explicitly in the quick start,
while other estimates are available for comparison. If an estimate is
missing or implausible, inspect the plots and choose a justified value
rather than silently substituting the default `10`.

This threshold concerns reads per corrected UMI. It is distinct from any
separate cell-barcode selection based on total reads per cell.

## Add calls to a Seurat object {#add-calls-to-a-seurat-object}

This optional example assumes `seurat_obj` is an existing Seurat object
for the same sample and its cell names match `flags$barcode_10x`.

`flag_snv()` returns both barcode formats:

| Output column | Example              | Format                                            |
|------------------------|------------------------|------------------------|
| `cbc`         | `AAACCCAAGAAACCAT`   | Cell barcode without a suffix.                    |
| `barcode_10x` | `AAACCCAAGAAACCAT-1` | Cell barcode with an automatically appended `-1`. |

HERMELIN constructs `barcode_10x` using `paste0(cbc, "-1")`. It does not
restore original sample prefixes or other suffixes. The join below
therefore uses `barcode_10x` for Seurat cells named in the usual
`BARCODE-1` format. Use `flags_phased` instead of `flags` if you want
the phased-session results.

``` r
metadata <- as.data.frame(flags)
stopifnot(!anyDuplicated(metadata$barcode_10x))
rownames(metadata) <- metadata$barcode_10x

cell_names <- colnames(seurat_obj)
idx <- match(cell_names, rownames(metadata))
if (all(is.na(idx))) {
    stop("No matching barcodes. Check sample prefixes and barcode suffixes.")
}
message(sum(!is.na(idx)), " of ", length(cell_names), " Seurat cells matched")

# Align by barcode, never by row position. Unmatched cells retain NA.
call_columns <- c(
    "mutation_call", "major_call_simple", "phased_call", "phased_call_simple"
)
metadata <- metadata[idx, call_columns, drop = FALSE]
rownames(metadata) <- cell_names
colnames(metadata) <- paste0("hermelin_", colnames(metadata))

seurat_obj <- SeuratObject::AddMetaData(
    object = seurat_obj,
    metadata = metadata
)
```

For renamed cells, map each result to its original sample and barcode
before joining. Do not remove sample prefixes from a merged object
without retaining sample identity, since different samples can reuse the
same barcode. See the [SeuratObject metadata
documentation](https://satijalab.github.io/seurat-object/reference/AddMetaData.html).

## Output files {#output-files}

| File                                   | Written by        | Contents                                                                                                                                                                                   |
|------------------------|------------------------|------------------------|
| `<session>_snv_table.csv`              | `detect_snv()`    | One row per retained read, including `read_id`, `CBC`, `UMI`, `strand`, and target/minor bases and statuses. `enriched = TRUE` adds diagnostic columns.                                    |
| `<session>_log.csv`                    | `detect_snv()`    | Run parameters and read-filtering counts.                                                                                                                                                  |
| `<session>_snv_summary.csv`            | `summarize_snv()` | One row per cell. Semicolon-separated columns contain comma-separated counts per corrected UMI. Read with `data.table::fread(path, sep = ";")`, or pass the file directly to `flag_snv()`. |
| `<session>_umi_frequencies.csv`        | `summarize_snv()` | Ranked reads per corrected UMI.                                                                                                                                                            |
| `<session>_thresholds.csv`             | `summarize_snv()` | Threshold estimates and suggested range.                                                                                                                                                   |
| `<session>_major_roi_frequency.csv`    | `summarize_snv()` | Target allele frequencies.                                                                                                                                                                 |
| `<session>_minor_roi_frequency.csv`    | `summarize_snv()` | Marker allele frequencies, when minor-ROI data are available.                                                                                                                              |
| `<session>_snp_frequency_combined.csv` | `summarize_snv()` | Combined target and marker frequency tables, where available.                                                                                                                              |
| `<session>_umi_distribution.pdf`       | `summarize_snv()` | Diagnostic plots, unless `skip_plots = TRUE`.                                                                                                                                              |
| `<session>_flags.csv`                  | `flag_snv()`      | Per-cell calls, written when `path_output_folder` is supplied.                                                                                                                             |

The final table contains `cbc`, `barcode_10x`, `has_np_data`,
`mutation_call`, `major_call_simple`, `phased_call`, and
`phased_call_simple`.

## Function reference {#function-reference}

All arguments of the three pipeline functions. "Required" means there is
no default and the argument must be supplied. The examples above explain
how to choose the main settings.

### `detect_snv()`

Streams the FASTQ in chunks. For every read it searches for the Read 1
adapter on both strands, extracts the CBC and UMI immediately
downstream, checks for the configured homopolymer, and then locates the
major (and minor) ROI with a two-pass strategy: an exact IUPAC match
first, then a fuzzy match followed by pairwise alignment for reads that
did not match exactly.

#### Session and I/O

| Argument             | Default  | Description                                                                                                                                                                                                       |
|--------------|--------------|--------------------------------------------|
| `session_name`       | required | Prefix for all output files. Letters, digits, `_` and `-` only.                                                                                                                                                   |
| `path_input_file`    | required | Path to the FASTQ file (`.fastq` or `.fastq.gz`).                                                                                                                                                                 |
| `path_output_folder` | required | Existing directory where `<session>_snv_table.csv` and `<session>_log.csv` are written.                                                                                                                           |
| `n_cores`            | required | CPU cores for the pairwise-alignment pass. Automatically capped at `parallel::detectCores() - 1`.                                                                                                                 |
| `n_reads_per_chunk`  | `1e5`    | Reads per streaming chunk. Larger chunks are faster but use more memory.                                                                                                                                          |
| `enriched`           | `FALSE`  | If `TRUE`, appends 22 diagnostic columns (read length, adapter/ROI match widths, substitution/insertion/deletion counts, Q-scores, match method). Use for QC and publication analyses; not needed for genotyping. |

#### Read 1 adapter and barcode layout

| Argument             | Default  | Description                                                                                                                                                       |
|--------------|--------------|--------------------------------------------|
| `read1_seq`          | required | Adapter sequence that immediately precedes the cell barcode (A/C/G/T only), e.g. `CTACACGACGCTCTTCCGATCT` for 10x 3' libraries.                                   |
| `read1_mismatch`     | required | Maximum edit distance (substitutions + indels) allowed when matching the adapter. Reads with no match, or with a match on both strands (chimeric), are discarded. |
| `read1_within_range` | `70`     | Number of bases from the read start (on each strand) in which the adapter is searched.                                                                            |
| `len_CBC`            | `16`     | Cell barcode length extracted directly after the adapter match. Downstream whitelist matching currently assumes 16-base barcodes.                                 |
| `len_UMI`            | `12`     | UMI length, taken directly after the barcode. Match this to the library chemistry.                                                                                |

#### Poly(N) tail filter

| Argument         | Default | Description                                                                                                                                                                    |
|-------------|-------------|----------------------------------------------|
| `polyN_base`     | `"T"`   | Base of the homopolymer searched for downstream of the UMI (poly(T) for 3' GEX). Set to `NULL` to disable the filter, e.g. for 5' GEX libraries where the TSO follows the UMI. |
| `polyN_length`   | `10`    | Length of the homopolymer match sought in the downstream search window. Ignored when `polyN_base = NULL`.                                                                      |
| `polyN_mismatch` | `2`     | Number of non-`polyN_base` positions tolerated within that window. Ignored when `polyN_base = NULL`.                                                                           |

#### Major ROI (the SNV to genotype)

| Argument             | Default  | Description                                                                                                                                    |
|--------------|--------------|--------------------------------------------|
| `major_roi_seq`      | required | Flanking sequence with exactly one `N` at the SNV position. IUPAC codes (R, Y, S, W, K, M, B, D, H, V) are allowed elsewhere.                  |
| `major_roi_mismatch` | required | Maximum edit distance for the fuzzy ROI search (pass 2). Reads that match more than once are left uncalled.                                    |
| `major_roi_mut`      | `NULL`   | Base that labels a read `MUT` in `major_status`. Must be given together with `major_roi_wt`.                                                   |
| `major_roi_wt`       | `NULL`   | Base that labels a read `WT` in `major_status`. If both are `NULL`, only the raw `major_base` is reported and `summarize_snv()` cannot be run. |

#### Minor ROI (optional phasing SNP)

| Argument             | Default | Description                                                                                                                                  |
|--------------|--------------|--------------------------------------------|
| `minor_roi_seq`      | `NULL`  | Flanking sequence of a heterozygous SNP on the same transcript, with one `N` at the SNP position. `NULL` skips minor ROI detection entirely. |
| `minor_roi_mismatch` | `NULL`  | Maximum edit distance for the minor ROI search. Required when `minor_roi_seq` is given.                                                      |
| `minor_roi_mut`      | `NULL`  | Marker allele labeled `MUT` in `minor_status`, not necessarily a disease-associated allele. Must be given together with `minor_roi_wt`.      |
| `minor_roi_wt`       | `NULL`  | Other marker allele, labeled `WT` in `minor_status`. These labels distinguish the two SNP alleles.                                           |

**Returns** `NULL` invisibly; results are written to disk.

### `summarize_snv()`

Loads the SNV table, keeps reads with a `MUT`/`WT` major status and a
whitelisted barcode, collapses UMIs within each cell that are within
`umi_mismatch` edits of a more abundant UMI, counts `MUT`/`WT` reads per
corrected UMI, and estimates reads-per-UMI thresholds with several
knee-point methods (geometric knee, inflection, UIK, CellRanger-style,
Otsu, Gaussian mixture).

| Argument             | Default  | Description                                                                                                                                                                                                           |
|--------------|--------------|--------------------------------------------|
| `session_name`       | required | Prefix for output files.                                                                                                                                                                                              |
| `path_snv_table`     | required | Path to `<session>_snv_table.csv` written by `detect_snv()`. Only the columns needed for summarization are read, so enriched tables are fine.                                                                         |
| `path_barcodes`      | required | Cell-barcode whitelist (`.csv` read as the first column, otherwise one barcode per line). Whitelist entries are normalized to their first 16 bases and matched against extracted CBCs. Nonmatching reads are dropped. |
| `path_output_folder` | required | Existing directory for the summary, frequency tables, thresholds and diagnostic PDF.                                                                                                                                  |
| `umi_mismatch`       | `2L`     | Maximum Levenshtein distance for UMI collapsing. A UMI within this distance of a more frequent UMI in the same cell is merged into it.                                                                                |
| `n_cores`            | `4L`     | CPU cores. Capped at `parallel::detectCores() - 1`.                                                                                                                                                                   |
| `skip_plots`         | `FALSE`  | If `TRUE`, do not write `<session>_umi_distribution.pdf`.                                                                                                                                                             |

**Returns** (invisibly) a list with:

| Element           | Contents                                                                                                                                                |
|-----------------|-------------------------------------------------------|
| `summary`         | Per-cell `data.table` (same content as `<session>_snv_summary.csv`); can be passed directly to `flag_snv()`                                             |
| `umi_frequencies` | Reads per corrected UMI, ranked                                                                                                                         |
| `knee_points`     | List of threshold estimates: `upper_knee`, `inflection`, `lower_knee`, `uik`, `cellranger`, `otsu`, `mixture`, plus the suggested `lower`/`upper` range |
| `snp_frequencies` | Major/minor allele frequency tables                                                                                                                     |
| `phasing_viable`  | `TRUE` when both minor alleles each account for at least 10% of the minor-ROI UMI counts. A feasibility heuristic, not proof of the phase relationship  |
| `phasing_message` | Human-readable phasing assessment                                                                                                                       |
| `has_minor_roi`   | Whether minor ROI data was found in the input                                                                                                           |
| `statistics`      | Read and cell counts at each filtering stage                                                                                                            |

### `flag_snv()`

Assigns a genotype to every cell. Per cell, each UMI is kept if its
total reads (`MUT` + `WT`) reach `threshold`, called `MUT` or `WT` if at
least `consensus_threshold` of its reads agree, and the cell call is
derived from the set of called UMIs. With phasing enabled, each called
UMI is additionally qualified by its minor SNP status.

| Argument              | Default       | Description                                                                                                                                                                                                                             |
|--------------|--------------|--------------------------------------------|
| `path_snv_summary`    | required      | Path to `<session>_snv_summary.csv`, or the `summary` `data.table` returned by `summarize_snv()`.                                                                                                                                       |
| `path_barcodes`       | `NULL`        | Optional full barcode whitelist used for **backfilling**: barcodes in this list that are absent from the filtered summary are added with `mutation_call = "no_np_data"`. This does not establish that no raw reads existed.             |
| `path_barcode_filter` | `NULL`        | Optional barcode list used for **filtering**: only barcodes in this list are analysed and reported.                                                                                                                                     |
| `path_output_folder`  | `NULL`        | Directory for `<session>_flags.csv`. If `NULL`, nothing is written and the result is only returned.                                                                                                                                     |
| `session_name`        | `"snv_flags"` | Prefix for the output file.                                                                                                                                                                                                             |
| `threshold`           | `10`          | Minimum total target-called reads (`MUT` + `WT`) per corrected UMI, inclusive. Inspect the estimates and diagnostic plots from `summarize_snv()` before choosing a value. This is a reads-per-UMI threshold, not a per-cell VAF cutoff. |
| `consensus_threshold` | `0.8`         | Fraction of a UMI's reads that must support one allele for the UMI to be called `MUT` or `WT`. UMIs without sufficient consensus are ignored. A cell with no remaining called UMIs is reported as `below_threshold`.                    |
| `use_phasing`         | `FALSE`       | Use the minor SNP columns to produce `phased_call`. Falls back to unphased calls (with a warning) if the summary has no SNP columns.                                                                                                    |
| `minor_percentage`    | `0.2`         | A UMI's minor SNP status is only used if its minor reads are at least this fraction of its major reads; otherwise the minor status is `none`.                                                                                           |
| `minor_on_major`      | `"MUT"`       | Which minor SNP allele (`"MUT"` or `"WT"`) is expected on the same transcript copy as the major mutation. Determines whether an observed minor SNP confirms (`MUT_mnsnp`, `def_WT`) or contradicts (`MUT_ERROR`) the major call.        |

**Returns** a `data.table` with one row per cell: `cbc`, `barcode_10x`,
`has_np_data`, `mutation_call`, `major_call_simple`, `phased_call`,
`phased_call_simple` (see [Interpret genotype
calls](#interpret-genotype-calls)).

## License {#license}

MIT © 2026 Sebastian Eder

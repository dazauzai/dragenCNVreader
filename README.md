# dragenCNVreader

> Read and parse Illumina **DRAGEN** CNV VCF files into tidy data frames.

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/dazauzai/dragenCNVreader/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dazauzai/dragenCNVreader/actions/workflows/R-CMD-check.yaml)

`dragenCNVreader` provides a single function, `read_dragen_cnv_vcf()`,
that:

- reads one VCF file **or** every VCF in a folder (with optional recursion);
- splits the DRAGEN `FORMAT` field
  (`GT:CN:MCN:CNQ:MCNQ:CNF:MCNF:MF:SM:SD:MAF:BC:AS:PE`) into separate columns;
- maps the DRAGEN ALT codes
  (`<DEL>` / `<DUP>` / `<LOH>` / `<REF>` / `<INV>` / `<INS>` / `<BND>`,
  plus the VCF missing marker `.` which DRAGEN uses for reference
  segments) into a tidy `cnv` category
  (`loss` / `gain` / `cn-LOH` / `LOH` / `ref` / `inv` / `ins` / `bnd`),
  and **errors out** on any ALT outside that set;
- returns both the **PASS-only** call set and the **full** call set, plus
  the `##FORMAT=` lines from the VCF header;
- optionally writes the chosen table(s) to disk via `out_dir` / `out_type`,
  with overridable separator (`out_sep`) and header (`out_col`).

---

## Requirements

- **R ≥ 4.0.0**
- R packages (will be installed automatically by `remotes::install_github`):
  - [`data.table`](https://CRAN.R-project.org/package=data.table)
  - [`dplyr`](https://CRAN.R-project.org/package=dplyr)
  - [`tidyr`](https://CRAN.R-project.org/package=tidyr)
  - [`stringr`](https://CRAN.R-project.org/package=stringr)
  - [`magrittr`](https://CRAN.R-project.org/package=magrittr)

No system dependencies (gzipped VCFs are handled via R's built-in `gzfile`).

---

## Installation

### From GitHub (recommended)

```r
# install remotes if you don't have it
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

remotes::install_github("dazauzai/dragenCNVreader")
```

Equivalent via `devtools`:

```r
# install devtools if you don't have it
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

devtools::install_github("dazauzai/dragenCNVreader")
```

### From a local clone (for development)

```r
# in the project root
remotes::install_local(".")
# or:
devtools::install()
```

### From a release tarball

```r
install.packages("dragenCNVreader_0.1.0.tar.gz", repos = NULL, type = "source")
```

---

## Usage

```r
library(dragenCNVreader)
```

### 1. Read a single VCF

```r
res <- read_dragen_cnv_vcf(
  vcf_file = "sample.cnv.vcf.gz"
)

res$Passed_CNV    # PASS-only CNVs: chr, start, end, cnv, CN, sample_id
res$All_CNV       # all CNV records with FORMAT fields expanded
res$format_lines  # the ##FORMAT lines from the VCF header
```

### 2. Read every VCF in a folder (non-recursive)

```r
res <- read_dragen_cnv_vcf(
  vcf_dir   = "cnv_vcfs/",
  recursive = FALSE
)

# per-sample, keyed by sample id from the VCF
res$HSS2675$Passed_CNV
res$HSS2675$All_CNV
res$HSS2675$format_lines

# concatenated across all samples
res$all$Passed_CNV
res$all$All_CNV
res$all$format_lines   # named list, one entry per sample
```

### 3. Read all VCFs including sub-folders

```r
res <- read_dragen_cnv_vcf(
  vcf_dir   = "cnv_vcfs/",
  recursive = TRUE
)
```

### 4. Custom file suffix

```r
res <- read_dragen_cnv_vcf(
  vcf_dir   = "cnv_vcfs/",
  suffix    = "cnv.vcf.gz",
  recursive = TRUE
)
```

You can pass multiple comma-separated suffixes; the default is
`"vcf,vcf.gz"`.

### 5. Write parsed tables to disk

`read_dragen_cnv_vcf()` can dump the parsed tables to a directory for
downstream tools that prefer flat files. The default writer is

```r
write.table(file = ..., row.names = FALSE, col.names = FALSE,
            sep = "\t", quote = FALSE)
```

You only need to set `out_dir` to enable it.

```r
# Per-sample All_CNV TSVs (one file per VCF)
read_dragen_cnv_vcf(
  vcf_dir  = "cnv_vcfs/",
  out_dir  = "cnv_tsv/",
  out_type = "All_CNV"          # default
)

# Both All_CNV and Passed_CNV, CSV with a header row
read_dragen_cnv_vcf(
  vcf_dir  = "cnv_vcfs/",
  out_dir  = "cnv_csv/",
  out_type = c("All_CNV", "pass"),
  out_sep  = ",",               # switches extension to .csv
  out_col  = TRUE               # write a header line
)
```

`out_type` is case-insensitive and accepts synonyms:

| value                            | table written      |
|----------------------------------|--------------------|
| `"All_CNV"` / `"all"`            | full call set      |
| `"Passed_CNV"` / `"passed"` / `"pass"` | PASS-only set |

Unknown values raise an error. Files are named
`<sample_id>.<table>.tsv` (or `.csv` when `out_sep = ","`). The
directory is created if it does not already exist.

---

## Output schema

`Passed_CNV` (a `data.frame`)

| column     | type    | description                                |
|------------|---------|--------------------------------------------|
| chr        | char    | chromosome (e.g. `chr1`)                   |
| start      | integer | 1-based start (from VCF `POS`)             |
| end        | integer | end coordinate (parsed from `ID`)          |
| cnv        | char    | one of `loss` / `gain` / `cn-LOH` / `LOH` / `ref` / `inv` / `ins` / `bnd` |
| CN         | numeric | total copy number from FORMAT/CN           |
| sample_id  | char    | sample name (the column after FORMAT)      |

`All_CNV` includes all the columns above plus `QUAL`, `FILTER`, `ALT`,
and the full FORMAT fields:
`GT, CN, MCN, CNQ, MCNQ, CNF, MCNF, MF, SM, SD, MAF, BC, AS, PE`.

`format_lines` is a character vector (single-file mode) or a named list
of character vectors (batch mode) holding the raw `##FORMAT=` header
lines.

---

## ALT → cnv mapping

Any ALT not in this table raises an error (with the offending value
and the file name) so silent misclassification is impossible.

| DRAGEN ALT | Condition  | `cnv` value |
|------------|------------|-------------|
| `<DEL>`    | —          | `loss`      |
| `<DUP>`    | —          | `gain`      |
| `<LOH>`    | `CN == 2`  | `cn-LOH`    |
| `<LOH>`    | otherwise  | `LOH`       |
| `<REF>` or `.` ¹ | —    | `ref`       |
| `<INV>`    | —          | `inv`       |
| `<INS>`    | —          | `ins`       |
| `<BND>`    | —          | `bnd`       |
| anything else | —       | **error**   |

¹ DRAGEN's WGS CNV caller emits reference (no-CNV) segments with the
VCF missing-value marker `.` in the `ALT` column rather than `<REF>`;
both encodings are accepted.

---

## Citing

If this package helped your analysis, please cite the repository URL
and the DRAGEN version that produced your VCFs.

---

## Issues / Contributions

PRs and issues welcome at
<https://github.com/dazauzai/dragenCNVreader/issues>.

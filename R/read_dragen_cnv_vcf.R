#' Read and parse DRAGEN CNV VCF files
#'
#' Reads CNV (copy-number variation) VCF files produced by the Illumina
#' DRAGEN pipeline, parses the multi-key FORMAT field into separate
#' columns, normalises the DRAGEN ALT codes (DEL / DUP / LOH) into a tidy
#' `cnv` category (loss / gain / cn-LOH / LOH), and returns both the
#' PASS-only and the full call set as data frames.
#'
#' Single-file and recursive batch modes are supported.
#'
#' @param vcf_file  Character. Path to one VCF or VCF.GZ file. Mutually
#'   exclusive with `vcf_dir`.
#' @param vcf_dir   Character. Path to a directory containing VCF files.
#'   When set, every file matching `suffix` (optionally recursively) is
#'   read and combined.
#' @param suffix    Character. Comma-separated list of acceptable file
#'   suffixes. Default `"vcf,vcf.gz"`. Substrings are matched with
#'   `endsWith()` so a custom value like `"cnv.vcf.gz"` works too.
#' @param recursive Logical. When `vcf_dir` is set, whether to recurse
#'   into sub-directories. Default `FALSE`.
#'
#' @return
#' If `vcf_file` is supplied, a `list` with three named elements:
#' \describe{
#'   \item{`Passed_CNV`}{data.frame: chr, start, end, cnv, CN, sample_id;
#'         filtered to FILTER == "PASS" and ALT in DEL/DUP/LOH.}
#'   \item{`All_CNV`}{data.frame: all CNV records with the FORMAT fields
#'         expanded into separate columns (GT, CN, MCN, CNQ, MCNQ, CNF,
#'         MCNF, MF, SM, SD, MAF, BC, AS, PE) plus QUAL, FILTER, ALT,
#'         sample_id.}
#'   \item{`format_lines`}{The raw ##FORMAT=... lines from the VCF header.}
#' }
#'
#' If `vcf_dir` is supplied, a `list` whose elements are the per-sample
#' results above, plus an additional element `all` containing the
#' concatenated `Passed_CNV`, `All_CNV` and the per-sample `format_lines`.
#'
#' @details
#' The FORMAT/sample column is split on `:` into the keys named in the
#' FORMAT string, and pivoted into separate columns. Numeric-like
#' columns (CN, MCN, CNQ, MCNQ, MF, SM, SD, MAF, BC, AS) are coerced to
#' numeric. ALT values are stripped of `<>` brackets and mapped to the
#' tidy `cnv` category: `"DEL" -> "loss"`, `"DUP" -> "gain"`,
#' `"LOH"` with `CN == 2` -> `"cn-LOH"`, other `"LOH"` -> `"LOH"`.
#'
#' @examples
#' \dontrun{
#' # Single VCF
#' res <- read_dragen_cnv_vcf(vcf_file = "sample.cnv.vcf.gz")
#' head(res$Passed_CNV)
#' head(res$All_CNV)
#'
#' # All VCFs in a folder (non-recursive)
#' res <- read_dragen_cnv_vcf(vcf_dir = "cnv_vcfs/", recursive = FALSE)
#' head(res$all$Passed_CNV)
#'
#' # Custom suffix, recursive
#' res <- read_dragen_cnv_vcf(
#'   vcf_dir   = "cnv_vcfs/",
#'   suffix    = "cnv.vcf.gz",
#'   recursive = TRUE
#' )
#' }
#'
#' @export
#' @importFrom data.table fread
#' @importFrom dplyr "%>%" all_of any_of across bind_rows case_when filter left_join mutate rename row_number select
#' @importFrom magrittr "%>%"
#' @importFrom stats setNames
#' @importFrom stringr str_remove_all str_split str_trim
#' @importFrom tidyr pivot_wider unnest_longer
read_dragen_cnv_vcf <- function(vcf_file  = NULL,
                                vcf_dir   = NULL,
                                suffix    = "vcf,vcf.gz",
                                recursive = FALSE) {

  parse_one_vcf <- function(vcf_file) {

    # silence R CMD check no-visible-binding NOTEs
    `#CHROM` <- POS <- ID <- FORMAT <- sample_value <- NULL
    key      <- value <- row_id <- ALT <- CN <- ALT_clean <- NULL

    # 1. Read VCF header and extract FORMAT definition lines
    header       <- readLines(gzfile(vcf_file), n = 10000)
    format_lines <- grep("^##FORMAT=", header, value = TRUE)

    # 2. Read VCF body
    df <- data.table::fread(
      vcf_file,
      skip       = "#CHROM",
      header     = TRUE,
      sep        = "\t",
      data.table = FALSE
    ) %>%
      dplyr::rename(chr = `#CHROM`)

    # 3. Identify sample column (the one right after FORMAT)
    sample_col <- names(df)[which(names(df) == "FORMAT") + 1]
    sample_id  <- sample_col

    # 4. Add start, end, sample_id, row_id
    df <- df %>%
      dplyr::mutate(
        start     = as.integer(POS),
        end       = as.integer(sub(".*-", "", ID)),
        sample_id = sample_id,
        row_id    = dplyr::row_number()
      )

    # 5. Expand FORMAT/sample column into separate columns
    format_wide <- df %>%
      dplyr::select(row_id, FORMAT, sample_value = dplyr::all_of(sample_col)) %>%
      dplyr::mutate(
        key   = stringr::str_split(FORMAT, ":"),
        value = stringr::str_split(sample_value, ":")
      ) %>%
      tidyr::unnest_longer(c(key, value)) %>%
      tidyr::pivot_wider(names_from = key, values_from = value)

    # 6. Merge expanded FORMAT columns back
    df_parsed <- df %>%
      dplyr::left_join(format_wide, by = "row_id") %>%
      dplyr::select(-row_id)

    # 7. Make sure expected FORMAT columns exist
    expected_format_cols <- c(
      "GT", "CN", "MCN", "CNQ", "MCNQ", "CNF", "MCNF",
      "MF", "SM", "SD", "MAF", "BC", "AS", "PE"
    )
    for (col in expected_format_cols) {
      if (!col %in% names(df_parsed)) df_parsed[[col]] <- NA
    }

    # 8. Convert numeric-like FORMAT columns
    numeric_cols <- c(
      "CN", "MCN", "CNQ", "MCNQ",
      "MF", "SM", "SD", "MAF", "BC", "AS"
    )
    df_parsed <- df_parsed %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::any_of(numeric_cols),
          ~ suppressWarnings(as.numeric(.x))
        )
      )

    # 9. Standardise ALT
    df_parsed <- df_parsed %>%
      dplyr::mutate(ALT_clean = stringr::str_remove_all(ALT, "[<>]"))

    # 10. Convert DRAGEN ALT into CNV category
    df_parsed <- df_parsed %>%
      dplyr::mutate(
        cnv = dplyr::case_when(
          ALT_clean == "DEL"               ~ "loss",
          ALT_clean == "DUP"               ~ "gain",
          ALT_clean == "LOH" & CN == 2     ~ "cn-LOH",
          ALT_clean == "LOH"               ~ "LOH",
          TRUE                             ~ ALT_clean
        )
      )

    # 11. All_CNV table
    All_CNV <- df_parsed %>%
      dplyr::select(
        chr, start, end, cnv, CN, QUAL, FILTER, ALT,
        dplyr::all_of(expected_format_cols), sample_id
      )

    # 12. Passed_CNV table
    Passed_CNV <- df_parsed %>%
      dplyr::filter(FILTER == "PASS",
                    ALT_clean %in% c("DEL", "DUP", "LOH")) %>%
      dplyr::select(chr, start, end, cnv, CN, sample_id)

    list(
      Passed_CNV   = Passed_CNV,
      format_lines = format_lines,
      All_CNV      = All_CNV
    )
  }

  # ----------------------------------------------------------
  # Batch mode: scan a directory and process all VCF / VCF.GZ files
  # ----------------------------------------------------------
  if (!is.null(vcf_dir)) {

    suffix_vec <- stringr::str_split(suffix, ",")[[1]] %>% stringr::str_trim()

    all_files <- list.files(path = vcf_dir, recursive = recursive,
                            full.names = TRUE)
    if (length(all_files) == 0) {
      stop("No VCF files found in the input directory with suffix: ", suffix)
    }
    keep <- vapply(all_files, function(x)
      any(vapply(suffix_vec, function(suf) endsWith(x, suf), logical(1))),
      logical(1))
    vcf_files <- all_files[keep]

    if (length(vcf_files) == 0)
      stop("No VCF files found in the input directory with suffix: ", suffix)

    result_list  <- lapply(vcf_files, parse_one_vcf)
    sample_names <- sapply(result_list,
                           function(x) unique(x$All_CNV$sample_id)[1])
    names(result_list) <- make.unique(sample_names)

    all_unit <- list(
      Passed_CNV   = dplyr::bind_rows(lapply(result_list, function(x) x$Passed_CNV)),
      All_CNV      = dplyr::bind_rows(lapply(result_list, function(x) x$All_CNV)),
      format_lines = stats::setNames(
        lapply(result_list, function(x) x$format_lines),
        names(result_list)
      )
    )
    result_list$all <- all_unit
    return(result_list)
  }

  # ----------------------------------------------------------
  # Single-file mode
  # ----------------------------------------------------------
  if (!is.null(vcf_file)) return(parse_one_vcf(vcf_file))

  stop("Please provide either `vcf_file` or `vcf_dir`.")
}

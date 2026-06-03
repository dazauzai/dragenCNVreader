#' Read and parse DRAGEN CNV VCF files
#'
#' Reads CNV (copy-number variation) VCF files produced by the Illumina
#' DRAGEN pipeline, parses the multi-key FORMAT field into separate
#' columns, normalises the DRAGEN ALT codes into a tidy `cnv`
#' category, and returns both the PASS-only and the full call set as
#' data frames. Optionally writes the chosen table to disk.
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
#' @param out_dir   Character or `NULL`. When set, the chosen table(s)
#'   (see `out_type`) are written to this directory using
#'   `utils::write.table()`. The directory is created if it does not
#'   exist. Default `NULL` (no files written).
#' @param out_type  Character vector. Which table(s) to write when
#'   `out_dir` is set. Case-insensitive; accepts `"All_CNV"` /
#'   `"all"`, and `"Passed_CNV"` / `"passed"` / `"pass"`. Multiple
#'   values are allowed. Default `"All_CNV"`. Unknown values raise an
#'   error.
#' @param out_sep   Character. Field separator passed to `write.table()`
#'   when writing output. Default `"\t"`. If set to `","` the file
#'   extension switches to `.csv`, otherwise `.tsv`.
#' @param out_col   Logical. Value of `col.names` passed to
#'   `write.table()` when writing output. Default `FALSE`.
#'
#' @return
#' If `vcf_file` is supplied, a `list` with three named elements:
#' \describe{
#'   \item{`Passed_CNV`}{data.frame: chr, start, end, cnv, CN, sample_id;
#'         filtered to FILTER == "PASS" and ALT in DEL/DUP/LOH (i.e.
#'         actual CNV events; REF segments and INV/INS/BND structural
#'         variants are excluded — they remain in All_CNV).}
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
#' numeric.
#'
#' ALT values are stripped of `<>` brackets and validated against the
#' supported set `.` (DRAGEN's REF segment marker), `DEL`, `DUP`,
#' `LOH`, `REF`, `INV`, `INS`, `BND`. Any other value raises an error.
#' Supported ALTs are mapped to the `cnv` category:
#'
#' * `DEL` -> `loss`
#' * `DUP` -> `gain`
#' * `LOH` with `CN == 2` -> `cn-LOH`
#' * `LOH` otherwise -> `LOH`
#' * `REF` or `.` -> `ref`
#' * `INV` -> `inv`
#' * `INS` -> `ins`
#' * `BND` -> `bnd`
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
#'
#' # Write per-sample All_CNV tables to disk as TSV
#' read_dragen_cnv_vcf(
#'   vcf_dir  = "cnv_vcfs/",
#'   out_dir  = "cnv_tsv/",
#'   out_type = "All_CNV"
#' )
#'
#' # Write both All_CNV and Passed_CNV as CSV, with a header row
#' read_dragen_cnv_vcf(
#'   vcf_dir  = "cnv_vcfs/",
#'   out_dir  = "cnv_csv/",
#'   out_type = c("All_CNV", "pass"),
#'   out_sep  = ",",
#'   out_col  = TRUE
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
#' @importFrom utils write.table
read_dragen_cnv_vcf <- function(vcf_file  = NULL,
                                vcf_dir   = NULL,
                                suffix    = "vcf,vcf.gz",
                                recursive = FALSE,
                                out_dir   = NULL,
                                out_type  = "All_CNV",
                                out_sep   = "\t",
                                out_col   = FALSE) {

  # Canonical ALT -> cnv mapping. DEL/DUP/LOH keep their original
  # semantic names (loss/gain/cn-LOH/LOH); newer ALT codes are lower-cased.
  # DRAGEN encodes reference (no-CNV) segments with the VCF missing
  # marker `.` rather than `<REF>`, so both are accepted as REF.
  # Extend `allowed_alts` here to support new ALT codes in future.
  allowed_alts <- c(".", "DEL", "DUP", "LOH", "REF", "INV", "INS", "BND")

  # out_type lookup. Keys are lower-case; values are slot names in the
  # per-sample result list. Add new (key, slot) pairs here to expose
  # additional output tables.
  out_type_map <- list(
    "all_cnv"    = "All_CNV",
    "all"        = "All_CNV",
    "passed_cnv" = "Passed_CNV",
    "passed"     = "Passed_CNV",
    "pass"       = "Passed_CNV"
  )

  resolve_out_types <- function(out_type) {
    keys <- tolower(out_type)
    bad  <- keys[!keys %in% names(out_type_map)]
    if (length(bad) > 0) {
      stop("Unknown out_type value(s): ",
           paste(unique(bad), collapse = ", "),
           ". Allowed: ",
           paste(unique(names(out_type_map)), collapse = ", "))
    }
    unique(unlist(out_type_map[keys], use.names = FALSE))
  }

  write_one_result <- function(result, sample_name, out_dir, slots,
                               out_sep, out_col) {
    ext <- if (identical(out_sep, ",")) ".csv" else ".tsv"
    for (slot in slots) {
      tbl <- result[[slot]]
      if (is.null(tbl)) next
      out_file <- file.path(out_dir,
                            paste0(sample_name, ".", slot, ext))
      utils::write.table(tbl, file = out_file,
                         row.names = FALSE, col.names = out_col,
                         sep = out_sep, quote = FALSE)
    }
  }

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

    # 9. Standardise ALT and validate against the supported set
    df_parsed <- df_parsed %>%
      dplyr::mutate(ALT_clean = stringr::str_remove_all(ALT, "[<>]"))

    bad_alts <- setdiff(unique(df_parsed$ALT_clean), allowed_alts)
    if (length(bad_alts) > 0) {
      stop("Unsupported ALT value(s) in '", vcf_file, "': ",
           paste(bad_alts, collapse = ", "),
           ". Allowed: ", paste(allowed_alts, collapse = ", "))
    }

    # 10. Convert DRAGEN ALT into the cnv category
    df_parsed <- df_parsed %>%
      dplyr::mutate(
        cnv = dplyr::case_when(
          ALT_clean == "DEL"             ~ "loss",
          ALT_clean == "DUP"             ~ "gain",
          ALT_clean == "LOH" & CN == 2   ~ "cn-LOH",
          ALT_clean == "LOH"             ~ "LOH",
          ALT_clean %in% c("REF", ".")   ~ "ref",
          ALT_clean == "INV"             ~ "inv",
          ALT_clean == "INS"             ~ "ins",
          ALT_clean == "BND"             ~ "bnd"
        )
      )

    # 11. All_CNV table
    All_CNV <- df_parsed %>%
      dplyr::select(
        chr, start, end, cnv, CN, QUAL, FILTER, ALT,
        dplyr::all_of(expected_format_cols), sample_id
      )

    # 12. Passed_CNV table: PASS rows that are real CNV events.
    # REF segments (.) and the structural-variant ALTs (INV/INS/BND)
    # are still parsed into All_CNV but excluded here.
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
  # Validate out_* up front so we fail fast on bad config
  # ----------------------------------------------------------
  slots <- NULL
  if (!is.null(out_dir)) {
    slots <- resolve_out_types(out_type)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
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

    if (!is.null(out_dir)) {
      for (s in names(result_list)) {
        write_one_result(result_list[[s]], s, out_dir, slots,
                         out_sep, out_col)
      }
    }

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
  if (!is.null(vcf_file)) {
    res <- parse_one_vcf(vcf_file)
    if (!is.null(out_dir)) {
      sample_name <- unique(res$All_CNV$sample_id)[1]
      write_one_result(res, sample_name, out_dir, slots,
                       out_sep, out_col)
    }
    return(res)
  }

  stop("Please provide either `vcf_file` or `vcf_dir`.")
}

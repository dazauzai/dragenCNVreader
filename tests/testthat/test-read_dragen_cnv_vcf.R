# Minimal end-to-end test using a tiny synthetic DRAGEN-style CNV VCF
# that lives at inst/extdata/mini.cnv.vcf.

test_that("single-file mode parses a tiny DRAGEN CNV VCF", {
  vcf <- system.file("extdata", "mini.cnv.vcf",
                     package = "dragenCNVreader")
  skip_if(vcf == "", "Example VCF not installed; skipping.")

  res <- read_dragen_cnv_vcf(vcf_file = vcf)

  expect_named(res, c("Passed_CNV", "format_lines", "All_CNV"))

  # All_CNV should have FORMAT fields expanded
  expect_true(all(c("chr", "start", "end", "cnv", "CN", "QUAL",
                    "FILTER", "ALT", "GT", "MCN", "BC", "sample_id")
                  %in% names(res$All_CNV)))

  # Passed_CNV should be a subset with only the tidy columns
  expect_named(res$Passed_CNV,
               c("chr", "start", "end", "cnv", "CN", "sample_id"))

  # All Passed rows must have FILTER == PASS in All_CNV
  pass_keys <- paste(res$Passed_CNV$chr,
                     res$Passed_CNV$start,
                     res$Passed_CNV$end)
  all_pass <- paste(
    res$All_CNV$chr[res$All_CNV$FILTER == "PASS"],
    res$All_CNV$start[res$All_CNV$FILTER == "PASS"],
    res$All_CNV$end[res$All_CNV$FILTER == "PASS"]
  )
  expect_true(all(pass_keys %in% all_pass))

  # cnv categorisation: DEL/DUP/LOH keep semantic names; the new
  # codes added in 0.2.0 are lower-case
  expect_true(all(res$Passed_CNV$cnv %in%
                    c("loss", "gain", "cn-LOH", "LOH",
                      "ref", "inv", "ins", "bnd")))

  # format_lines must be character of ##FORMAT lines
  expect_type(res$format_lines, "character")
  expect_true(all(startsWith(res$format_lines, "##FORMAT=")))
})

test_that("input validation works", {
  expect_error(read_dragen_cnv_vcf(),
               "Please provide either")
  expect_error(read_dragen_cnv_vcf(vcf_dir = tempdir()),
               "No VCF files found")
})

test_that("ALT='.' (DRAGEN REF segments) is accepted and mapped to 'ref'", {
  vcf <- system.file("extdata", "mini.cnv.vcf",
                     package = "dragenCNVreader")
  skip_if(vcf == "", "Example VCF not installed; skipping.")

  with_ref <- tempfile(fileext = ".vcf")
  lines    <- readLines(vcf)
  # Append a REF-segment row in DRAGEN's actual encoding (ALT = ".")
  ref_row  <- paste(
    "chr5", "5000000", "DRAGEN:REF:chr5:5000000-5100000", "N", ".",
    "1000", "PASS", "END=5100000",
    "GT:CN:MCN:CNQ:MCNQ:CNF:MCNF:MF:SM:SD:MAF:BC:AS:PE",
    "0/0:2:1:1000:1000:2.0:1.0:.:1.0:0.05:0.5:100:50:0,0",
    sep = "\t"
  )
  writeLines(c(lines, ref_row), with_ref)

  res <- read_dragen_cnv_vcf(vcf_file = with_ref)
  expect_true("ref" %in% res$All_CNV$cnv)
  expect_true("ref" %in% res$Passed_CNV$cnv)
})

test_that("unsupported ALT codes raise an error", {
  vcf <- system.file("extdata", "mini.cnv.vcf",
                     package = "dragenCNVreader")
  skip_if(vcf == "", "Example VCF not installed; skipping.")

  bad <- tempfile(fileext = ".vcf")
  lines <- readLines(vcf)
  # Replace the first data row's ALT (<DEL>) with an unsupported code
  data_idx <- grep("^chr", lines)[1]
  lines[data_idx] <- sub("<DEL>", "<FOO>", lines[data_idx], fixed = TRUE)
  writeLines(lines, bad)

  expect_error(read_dragen_cnv_vcf(vcf_file = bad),
               "Unsupported ALT")
})

test_that("out_dir writes the requested tables", {
  vcf <- system.file("extdata", "mini.cnv.vcf",
                     package = "dragenCNVreader")
  skip_if(vcf == "", "Example VCF not installed; skipping.")

  out <- file.path(tempfile("cnvout_"))
  res <- read_dragen_cnv_vcf(
    vcf_file = vcf,
    out_dir  = out,
    out_type = c("All_CNV", "pass")
  )
  sample_name <- unique(res$All_CNV$sample_id)[1]

  expect_true(file.exists(file.path(out,
                                    paste0(sample_name, ".All_CNV.tsv"))))
  expect_true(file.exists(file.path(out,
                                    paste0(sample_name, ".Passed_CNV.tsv"))))
})

test_that("unknown out_type values raise an error", {
  vcf <- system.file("extdata", "mini.cnv.vcf",
                     package = "dragenCNVreader")
  skip_if(vcf == "", "Example VCF not installed; skipping.")

  expect_error(
    read_dragen_cnv_vcf(vcf_file = vcf,
                        out_dir  = tempfile("cnvout_"),
                        out_type = "nope"),
    "Unknown out_type"
  )
})

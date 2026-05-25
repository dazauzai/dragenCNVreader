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

  # cnv categorisation
  expect_true(all(res$Passed_CNV$cnv %in%
                    c("loss", "gain", "cn-LOH", "LOH")))

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

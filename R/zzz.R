# Silence R CMD check NOTEs for tidyverse NSE column references that aren't
# visible at static analysis time.
utils::globalVariables(c(
  "chr", "start", "end", "cnv",
  "QUAL", "FILTER", "ALT", "ALT_clean", "CN",
  "POS", "ID", "FORMAT", "sample_value",
  "key", "value", "row_id"
))

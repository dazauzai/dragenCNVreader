# dragenCNVreader 0.2.0

- `read_dragen_cnv_vcf()`:
  - **Breaking:** the `cnv` category is now lower-case throughout
    (`cn-loh` / `loh` instead of `cn-LOH` / `LOH`).
  - **Breaking:** unrecognised ALT codes now raise an error instead of
    falling through to `NA`. Allowed ALTs: `DEL`, `DUP`, `LOH`, `REF`,
    `INV`, `INS`, `BND`.
  - Added ALT support for `REF`, `INV`, `INS`, `BND` (mapped to `ref`,
    `inv`, `ins`, `bnd`).
  - New optional output parameters:
    - `out_dir` — directory to write parsed tables to.
    - `out_type` — which table(s) to write; accepts `"All_CNV"` /
      `"all"` and `"Passed_CNV"` / `"passed"` / `"pass"`. Vector input
      is allowed.
    - `out_sep` — override the column separator (default `"\t"`;
      `","` switches the extension to `.csv`).
    - `out_col` — override `col.names` (default `FALSE`).
  - Files are written as `<sample_id>.<table>.<ext>` using
    `utils::write.table(..., row.names = FALSE, col.names = out_col, sep = out_sep, quote = FALSE)`.

# dragenCNVreader 0.1.0

Initial release.

- `read_dragen_cnv_vcf()`:
  - single-file mode (`vcf_file`) and batch mode (`vcf_dir`, optional
    `recursive`);
  - parses `GT:CN:MCN:CNQ:MCNQ:CNF:MCNF:MF:SM:SD:MAF:BC:AS:PE`;
  - tidy `cnv` category (`loss` / `gain` / `cn-LOH` / `LOH`);
  - returns `Passed_CNV` (PASS-only), `All_CNV` (full), and the
    raw `##FORMAT=` header lines.
- Minimal end-to-end testthat suite with a synthetic 4-record VCF.
- MIT license, GitHub Actions R-CMD-check workflow.

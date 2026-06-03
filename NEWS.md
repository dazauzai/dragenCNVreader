# dragenCNVreader 0.2.2

- `read_dragen_cnv_vcf()`:
  - Restore the original `Passed_CNV` ALT filter (DEL / DUP / LOH).
    0.2.0 had dropped this filter, which let REF segments (and the
    new INV / INS / BND structural variants) leak into the PASS-only
    CNV table. They still appear in `All_CNV` — only the CNV-focused
    `Passed_CNV` view is restricted.

# dragenCNVreader 0.2.1

- `read_dragen_cnv_vcf()`:
  - Accept ALT `.` (VCF missing marker) as DRAGEN's reference-segment
    encoding and map it to `cnv = "ref"`. The WGS CNV caller emits
    REF segments this way rather than as `<REF>`, so 0.2.0 erroneously
    rejected real WGS VCFs.
  - Restore the original mixed-case `cn-LOH` / `LOH` for `<LOH>` ALTs.
    The 0.2.0 lower-casing was an over-reach of the
    "convert REF/INV/INS/BND to lower-case" request.

# dragenCNVreader 0.2.0

- `read_dragen_cnv_vcf()`:
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

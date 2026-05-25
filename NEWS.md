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

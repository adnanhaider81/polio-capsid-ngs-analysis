# Changelog

## 1.0.1 - 2026-04-30
- Fixed antigenic-site TSV writing and script compilation.
- Declared BAM/index outputs in the Snakemake DAG so downstream rules can be resolved reproducibly.
- Added documented coverage outputs to the default workflow target.
- Improved lightweight validation commands for clean external checkouts.
- Added repository hygiene rules for private data, generated outputs, and local environment files.

## 1.0.0 - 2025-09-29
- Initial public pipeline release for post-culture whole-capsid poliovirus NGS analysis.

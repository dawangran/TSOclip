# Changelog

## 0.2.0-dev

- Replaced shift-only Hamming matching with bounded edit-distance scoring for full and partial TSO hits.
- Added tail-anchored hit-chain selection to avoid trimming internal TSO-like matches that are not connected to the 3' end.
- Added explicit cut coordinates and edit/mismatch/insertion/deletion fields to TSV output.
- Added portable OpenMP auto-detection in the Makefile.
- Added shell-based regression tests, Docker build support, and GitHub Actions CI.

## 0.1.0

- Initial single-file C implementation for tail-window TSO detection and trimming.

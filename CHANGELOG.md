# Changelog

## Unreleased

### Added

- Added `code/scripts/gpFlatfile2SQLite.py`, a Python utility that parses `genomeProperties.txt` and writes a SQLite database using a schema aligned with the existing Genome Properties relational model.
- Added GO term enrichment during SQLite export via the EBI QuickGO API, with local caching and a `--skip-go-lookup` option for offline or faster runs.
- Added `gp_sqlite.yml` to provide a conda environment for the SQLite export workflow.
- Added a generated SQLite database artifact at `flatfiles/genome_properties.db`.

### Notes

- This branch introduces a SQLite-based distribution path for Genome Properties data alongside the existing flatfile assets.

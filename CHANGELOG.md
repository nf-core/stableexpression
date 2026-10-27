# nf-core/stableexpression: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - [2026-09-28]

### `Added`

- [#58](https://github.com/nf-core/stableexpression/issues/58) Tries to fetch genome annotation from the [Ensembl organisms](https://ftp.ebi.ac.uk/pub/ensemblorganisms/) FTP, in addition to Ensembl and NCBI.

### `Fixed`

- [#38](https://github.com/nf-core/stableexpression/issues/38) Fix bug with species having multiple entries in `g:Profiler` (like `canis lupus`). Now querying `g:Profiler` beforehand to get the exact list of available species
- [#37](https://github.com/nf-core/stableexpression/issues/37) Fix too short time limit for `IMPUTE_MISSING_VALUES` module, which resulting in module failture for species comprising a large number of datasets, like `homo sapiens`, `mus musculus`, or `arabidopsis thaliana`
- [#36](https://github.com/nf-core/stableexpression/issues/36) Check correspondance between Ensembl gene IDs obtained from `g:Profiler` and gene IDs contained in the downloaded genome annotation. Now multiple annotations are downloaded and the one displaying the greatest number of matching gene IDs is chosen. If no annotation displayed matching gene IDs, the module tries with NCBI genome annotation.
- [#61](https://github.com/nf-core/stableexpression/issues/61) Fixed an issue with the R Bioconductor ExpressionAtlas package resulting in a bug when trying to download some experiments (for example `Musa acuminata`).

### `Changed`

- [#PR46](https://github.com/nf-core/stableexpression/pull/46) Updated to nf-core pipeline template v4.1.0
- [#PR46](https://github.com/nf-core/stableexpression/pull/46) Changed default branch to `main`

## v1.0.0 - [2026-09-14]

Initial release of nf-core/stableexpression, created with the [nf-core](https://nf-co.re/) template.

## v1.0dev - [2025-01-26]

Initial pre-release of nf-core/stableexpression, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

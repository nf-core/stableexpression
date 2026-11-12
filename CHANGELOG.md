# nf-core/stableexpression: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0 - [2026-11-10]

Important refactoring of the normalisation process.

### `Added`

- A new `gene_max_zero_ratio` parameter allows setting a threshold to exclude genes showing a too high proportion of zeros across datasets in a specific platform
- The new `quantile_normalisation` parameter, allows users to select quantile normalisation instead of the default `reComBat` normalisation pathway (not recommended, unless scalability issues are encountered)
- The `multiqc_nb_genes_plotted` parameter to let users control how many genes have their expression plotted in each section of MultiQC report. Set default to 0 to avoid loading issues in the web browser when the report is too large.

### `Changed`

- [#59] All raw RNA-seq counts now undergo [GeTMM](https://link.springer.com/article/10.1186/s12859-018-2246-7) normalisation. RNA(seq data already normalised are converted to TPM.
- Instead of quantile normalisation, correction of batch effect is now performed using [reComBat](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac071/6728137?login=false) by default
- `nb_candidates_per_section` now defaults to `null`, meaning that all genes are kept as candidates for subsequent analyses. This improves the precision of stability scores.
- Replaced the `genorm` subworkflow by a single module, wrapper a highly optimised implementation of Genorm.
- Updated the metromap.

### `Removed`

- the `normalisation_method` parameter. All raw RNA-seq data now undergo GeTMM normalisation, while normalised data are converted to TPM.
- the `quantile_norm_target_distrib` parameter. Now quantile normalisation must use a uniform target distribution, to avoid negative values in the normalised dataframe.
- the `missing_value_imputer` parameter. Now [KNN imputation](https://scikit-learn.org/stable/modules/generated/sklearn.impute.KNNImputer.html) is mandatory (iterative imputation cannot hold scalability as well as knn imputation)
- the `skip_genorm` parameter. `GeNorm` is now mandatory.

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

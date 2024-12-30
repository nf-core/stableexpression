<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-stableexpression_logo_dark.png">
    <img alt="nf-core/stableexpression" src="docs/images/nf-core-stableexpression_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/stableexpression/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/stableexpression/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/stableexpression/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/stableexpression/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/stableexpression/results)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/stableexpression)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23stableexpression-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/stableexpression)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/stableexpression** is a bioinformatics pipeline that aims at finding the most stable genes among a single or multiple public / local count datasets. It takes as input a species name (mandatory), keywords for expression atlas search (optional) and / or a CSV input file listing local raw / normalized count datasets (optional). **A typical usage is to find the most suitable qPCR housekeeping genes for a specific species (and optionally specific conditions)**.

<p align="center">
    <img title="Stableexpression Workflow" src="docs/images/nf-core-stableexpression_metro_map.png" width=100%>
</p>

## Pipeline summary

1. Get Expression Atlas accessions corresponding to the provided species (and optionally keywords) ([Expression Atlas](https://www.ebi.ac.uk/gxa/home); optional)
2. Download Expression Atlas data ([Expression Atlas](https://www.ebi.ac.uk/gxa/home); optional)
3. Normalize raw data (using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or [EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html))
4. Map gene IDS to Ensembl IDS for standardisation among datasets ([g:Profiler](https://biit.cs.ut.ee/gprofiler/gost))
5. Compute gene variation coefficients and get the most stable genes

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet listing the different count datasets:

`datasets.csv`:

```csv
counts,design,normalized
path/to/normalized.counts.csv,path/to/normalized.design.csv,true
path/to/raw.counts.csv,path/to/raw.design.csv,false
```

Make sure to format your datasets properly:

`counts.csv`:

```csv
,sample_A,sample_B,sample_C
gene_1,1,2,3
gene_2,1,2,3
...
```

`design.csv`:

```csv
sample,condition
sample_A,condition_1
sample_B,condition_2
...
```

Now you can tun the pipeline as follows:

> ```bash
> nextflow run nf-core/stableexpression \
>   -profile docker \
>   --species <SPECIES> \
>   --eatlas_accessions <ACCESSIONS> \
>   --eatlas_keywords <KEYWORDS> \
>   --datasets ./datasets.csv \
>   --outdir ./results
> ```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/stableexpression/usage) and the [parameter documentation](https://nf-co.re/stableexpression/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/stableexpression/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/stableexpression/output).

## Credits

nf-core/stableexpression was originally written by Olivier Coen.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#stableexpression` channel](https://nfcore.slack.com/channels/stableexpression) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

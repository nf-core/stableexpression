<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-stableexpression_logo_dark.png">
    <img alt="nf-core/stableexpression" src="docs/images/nf-core-stableexpression_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new/nf-core/stableexpression)
[![GitHub Actions CI Status](https://github.com/nf-core/stableexpression/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/stableexpression/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/stableexpression/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/stableexpression/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/stableexpression/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.4.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.4.1)
[![run with apptainer](https://custom-icon-badges.demolab.com/badge/run%20with-apptainer-4545?logo=apptainer&color=teal&labelColor=000000)](https://apptainer.org/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?logo=singularity&labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/stableexpression)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23stableexpression-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/stableexpression)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/stableexpression** is a bioinformatics pipeline aiming to aggregate multiple count datasets (public / provided by the user) for a specific species and find the most stable genes.

It takes as main inputs :
  * a species name (mandatory)
  * keywords for Expression Atlas / GEO search (optional)
  * a CSV input file listing your own raw / normalised count datasets (optional).

**Use cases**:
  * **find the most suitable genes as RT-qPCR reference genes for a specific species (and optionally specific conditions)**
  * download all Expression Atlas and NCBI GEO datasets (microarray only!) for a species



<p align="center">
    <img title="Stableexpression Workflow" src="docs/images/nf-core-stableexpression_metro_map.png" width=100%>
</p>

## Basic usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

To search the most stable genes in a species considering all public datasets, simply run:

```bash
nextflow run nf-core/stableexpression \
   -r dev \
   -profile <PROFILE (examples: docker / apptainer / conda / micromamba)> \
   --species <SPECIES (examples: arabidopsis_thaliana / "drosophila melanogaster")> \
   --outdir <OUTDIR (example: ./results)>
 ```

> [!IMPORTANT]
 > For more specific scenarios, __like fetching only specific conditions or using your own expression dataset(s)__, please refer to the [usage documentation](https://nf-co.re/stableexpression/usage).

> [!NOTE]
> See [here](https://nf-co.re/stableexpression/usage#profiles) for more information about profiles.

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/stableexpression/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/stableexpression/output).

## Support us

If you like nf-core/stableexpression, please make sure you give it a star on GitHub.

[![stars - stableexpression](https://img.shields.io/github/stars/nf-core/stableexpression?style=social)](https://github.com/nf-core/stableexpression)

## Credits

nf-core/stableexpression was originally written by Olivier Coen.

<!-- We thank the following people for their extensive assistance in the development of this pipeline: -->

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#stableexpression` channel](https://nfcore.slack.com/channels/stableexpression) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/stableexpression for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

# nf-core/stableexpression: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/stableexpression/usage](https://nf-co.re/stableexpression/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Basic run

This pipeline fetches Expression Atlas and GEO accessions for the provided species and downloads the corresponding data.

```bash
nextflow run nf-core/stableexpression \
   -profile <PROFILE (examples: docker / apptainer / conda / micromamba)> \
   --species <SPECIES (examples: arabidopsis_thaliana / "drosophila melanogaster")> \
   --outdir <OUTDIR (example: ./results)> \
   -resume
```

## Specific public datasets

You can provide keywords to restrict downloaded datasets to specific conditions.

```bash
nextflow run nf-core/stableexpression \
   -profile <PROFILE> \
   --species <SPECIES> \
   --keywords <KEYWORDS (examples: "leaf" / "flower,stress")>
   --outdir <OUTDIR>
```

> [!NOTE]
>
> - Multiple keywords must be separated by commas.
> - Please note that keywords are additive, so you will receive datasets that match with **either of the provided keywords**.
> - A dataset will be downloaded if a keyword is found in its summary or in the name of a sample.
> - The natural language processing [`nltk`](https://www.nltk.org/) python package is used to find keywords as well as derived words. For example, the `leaf` keyword should match 'leaf', 'leaves', 'leafy', and so on.

Users can also choose to download data from a specific platform (`RNA-seq` / `Microarray`) only:

```bash
nextflow run nf-core/stableexpression \
   -profile <PROFILE> \
   --species <SPECIES> \
   --platform <rnaseq/microarray> \
   --outdir <OUTDIR>
```

## Provide your own accessions

You may already have specific Expression Atlas/GEO accession numbers in mind for use in the analysis.
If so, you can provide these directly to the pipeline.

```bash
nextflow run nf-core/stableexpression \
   -profile <PROFILE> \
   --species <SPECIES> \
   --skip_fetch_eatlas_accessions \
   [--eatlas_accessions <ACCESSION(S) (example: "E-MTAB-7711,E-GEOD-51720")>] \
   [--eatlas_accessions_file <FILE>] \
   [--geo_accessions <ACCESSION(S) (example: "GSE262492,GSE305365")>] \
   [--geo_accessions_file <FILE>] \
   --outdir <OUTDIR>
```

> [!NOTE]
>
> - To download only the datasets corresponding to the supplied accessions, set the `--skip_fetch_eatlas_accessions` parameter.
> - The files supplied with either the `--eatlas_accessions_file` or `--geo_accessions_file` parameters should contain one accession per line. The file extension does not matter.
> - For `homo sapiens`, the experiment `E-GTEX-8` comprises 17350 samples :exploding_head: :exploding_head: :exploding_head:... It it therefore excluded by default, but you can include it anyway by setting `--accessions E-GTEX-8`.

## Fetch accessions only

If you are unsure which accessions you want, but would like to control precisely which datasets are included in your analysis, you can first run:

```bash
nextflow run nf-core/stableexpression \
   -profile <PROFILE> \
   --species <SPECIES> \
   --accessions_only \
   --outdir <OUTDIR>
```

The list of fetched accessions with their respective metadata will be displayed in the MultiQC report, at `<OUTDIR>/reporting/multiqc_report.html`, in the sections `Expression Atlas datasets` and `GEO datasets`.

## Use your own expression datasets

Of course, you can provide your own count datasets and experimental designs.

> [!IMPORTANT]
>
> - To ensure that all RNA-seq datasets are processed in the same way, users should provide **raw counts**.
> - If normalised counts are provided, the same normalisation process must be applied to all of them.
> - Microarray data must already be normalised. When combining your own datasets with public ones in a single run, you should use the RMA method to ensure compliance with Expression Atlas datasets.

First, prepare a CSV samplesheet listing the datasets you want to use. Each row represents a specific dataset and must contain the following information:

| Column       | Description                                                                               |
| ------------ | ----------------------------------------------------------------------------------------- |
| `counts`     | Path to the count dataset (a CSV / TSV file)                                              |
| `design`     | Path to the experimental design file associated to this dataset (a CSV / TSV file)        |
| `platform`   | Platform used to generate the counts (`rnaseq` or `microarray`)                           |
| `normalised` | Boolean (`true` / `false`) representing whether the counts are already normalised or not. |

It should look as follows:

```csv title=datasets.csv
counts,design,platform,normalised
path/to/normalised.counts.csv,path/to/normalised.design.csv,rnaseq,true
path/to/raw.counts.csv,path/to/raw.design.csv,rnaseq,false
path/to/microarray.counts.csv,path/to/microarray.design.csv,microarray,true
```

It can also be a YAML file:

```yaml title=datasets.yaml
- counts: path/to/normalised.counts.csv
  design: path/to/normalised.design.csv
  platform: rnaseq
  normalised: true
- counts: path/to/raw.counts.csv
  design: path/to/raw.design.csv
  platform: rnaseq
  normalised: false
- counts: path/to/microarray.counts.csv
  design: path/to/microarray.design.csv
  platform: microarray
  normalised: true
```

The counts should have the following structure:

```csv title=counts.csv
gene_id,sample_A,sample_B,sample_C
gene_1,1,2,3
gene_2,1,2,3
```

The design should look like:

```csv title=design.csv
sample,condition
sample_A,condition_1
sample_B,condition_2
sample_C,condition_1
```

> [!NOTE]
>
> - The first column in the counts file must contain gene IDs and can be given any name (**as long as a comma is placed before the name of the first count column**). The count file should only contain the first column (gene IDs) and the sample columns.
> - Both counts and design files can also be supplied as TSV files.

Now run the pipeline with:

```bash
nextflow run nf-core/stableexpression \
   -profile <PROFILE> \
   --species <SPECIES> \
   --datasets <CSV / YAML FILE> \
   --skip_fetch_eatlas_accessions \
   --outdir <OUTDIR>
```

> [!NOTE]
>
> - In the example above, the `--skip_fetch_eatlas_accessions` parameter is provided to demonstrate how to analyse your own dataset only. If you want to mix your dataset(s) with public ones, you can remove this parameter.

> [!IMPORTANT]
> By default, the pipeline attempts to map gene IDs to Ensembl gene IDs. **All genes that cannot be mapped are discarded from the analysis**. This ensures that all genes are named consistently across datasets, enabling comparison of multiple datasets. However, if you are confident that the names of your genes are consistent across your different datasets, or if you think that your gene IDs simply won't map properly, you can disable this mapping by adding the '--skip_id_mapping' parameter. In this case, we recommend that users supply their own gene ID mapping and gene metadata files using the `--gene_id_mapping` and `--gene_metadata` parameters, respectively.
>
> Both files are totally optional, however:
> - a custom gene ID mapping might help to merge datasets properly.
> - custom gene metadata (the association between gene ID, gene name and gene description) will provide the MultiQC report with relevant metadata.
>
> See [next section](custom-gene-id-mapping--metadata) for further details.

> [!TIP]
> You can check if your gene IDs can be mapped using the [g:Profiler server](https://biit.cs.ut.ee/gprofiler/convert).

### Custom gene ID mapping / metadata

Users can provide their own gene ID mapping and/or gene metadata using the '--gene_id_mapping' and '--gene_metadata' parameters, respectively. The gene ID mapping file maps gene IDs in the count table(s) (local or downloaded) to more generic IDs, which are used as a basis for subsequent steps. The gene metadata file provides supplementary details about the genes, such as their common name and description.

Structure of the gene id mapping file:

| Column             | Description                                   |
| ------------------ | --------------------------------------------- |
| `original_gene_id` | Gene ID used in the provided count dataset(s) |
| `gene_id`          | Mapped gene ID                                |

Example:

```csv title=gene_id_mapping.csv
original_gene_id,gene_id
gene_A,ENSG1234567890
geneB,OTHERmappedgeneID
```

Structure of the gene metadata file:

| Column        | Description      |
| ------------- | ---------------- |
| `gene_id`     | Mapped gene ID   |
| `name`        | Gene common name |
| `description` | Gene description |

Example:

```csv title=gene_metadata.csv
gene_id,name,description
ENSG1234567890,Gene A,Description of gene A
OTHERmappedgeneID,My OTHER Gene,Another description
```

### RNA-seq count normalisation

By default, raw RNA-seq counts are normalised using [GeTMM](https://link.springer.com/article/10.1186/s12859-018-2246-7), which requires the knowledge of gene length. To achieve this, the pipeline first downloads the public genome annotation that best aligns with the gene IDs present in the datasets. It then computes the length of the longest transcript for each gene. If the species of interest does not have a public annotation or if you encounter network issues, you can provide your own genome annotation or a file associating gene IDs with gene lengths directly using the `--gff` (local), `--gff_url` (remote URL) and `--gene_length` parameters, respectively.

If users do not have access to genome annotations or gene lengths, they can use the `--skip_gene_length_normalisation` flag. The pipeline will then run edgeR without computing RPK normalisation. As the counts will not be normalised by length, this will result in a bias towards long genes.

The genome annotation must be in `GFF` format and have the `.gff` or `.gff3` extension. The [`AGAT`](https://github.com/NBISweden/AGAT) package can be used to convert other genome annotation formats to GFF.

The gene length file must be in `CSV` or `TSV` format and have the following structure:

| Column    | Description                      |
| --------- | -------------------------------- |
| `gene_id` | Mapped gene ID                   |
| `length`  | Gene length (longest transcript) |

Example:

```csv title=gene_length.csv
gene_id,length
ENSG1234567890,1000
OTHERmappedgeneID,2000
```

### More advanced scenarios

For advanced scenarios, you can see the list of available parameters in the [parameter documentation](https://nf-co.re/stableexpression/parameters).

## Pipeline output

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

For a detailed description of the output files, please consult the [nf-core stableexpression output directory structure](https://nf-co.re/stableexpression/output).

## Parameters

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run -r dev nf-core/stableexpression -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
species: 'Homo sapiens'
datasets: './datasets.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/stableexpression
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/stableexpression releases page](https://github.com/nf-core/stableexpression/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### [`-profile`](#profiles)

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Apptainer (Singularity) or Docker containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

> [!TIP]

> When running the pipeline of multi-user server or on a cluster, the best practice is to use Apptainer (formerly Singularity). You can install Apptainer by following these [instructions](https://apptainer.org/docs/admin/main/installation.html#).
> In case you encounter the following error when running Apptainer:
>
> ```
> ERROR  : Could not write info to setgroups: Permission denied
> ERROR  : Error while waiting event for user namespace mappings: no event received
> ```
>
> you may need to install the `apptainer-suid` package instead of `apptainer`:
>
> ```
> # Debian / Ubuntu
> sudo apt install apptainer-suid
> # RHEL / CentOS
> sudo yum install apptainer-suid
> # Fedora
> sudo dnf install apptainer-suid
> ```

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.
- `micromamba`
  - A faster, more lightweight alternative to Conda. As for Conda, use Micromamba as a last resort.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/main/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/main/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

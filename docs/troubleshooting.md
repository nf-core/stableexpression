# nf-core/stableexpression: Troubleshooting

## Ǹo dataset found

>[!IMPORTANT]
> For the time being, only Microarray count datasets are fetched from NCBI GEO.

For species that are not on Expression Atlas and that do not have microarray data on NCBI data, the pipeline will not be able to find suitable datasets and will log the following message:

```
WARN: No dataset found. Please note that for the moment only Microarray count datasets are fetched from NCBI GEO.
You can check at https://www.ncbi.nlm.nih.gov/gds if there are raw RNA-seq count datasets for this species.
```

You may want to check if there are any raw RNA-seq count datasets available for this species on [NCBI GEO](https://www.ncbi.nlm.nih.gov/gds). You can then relaunch the pipeline by providing your own count datasets.

## Java heap space

In some cases, in particular when running the pipeline on a very large number of datasets (such as for `Homo sapiens`), the Nextflow Java virtual machines can start to request a large amount of memory. You may happen to see the following error:

```
java.lang.OutOfMemoryError: Java heap space
```

We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):
```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

or running the pipeline with:
```bash
NXF_OPTS='-Xms1g -Xmx4g' nextflow run nf-core/stableexpression ...

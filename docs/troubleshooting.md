# nf-core/stableexpression: Troubleshooting

## Error 139 on macOS

If you are running the pipeline on macOS with containers (`docker`, `apptainer`, `singularity`, ...), you may encounter issues like:

```
NOTE: Process `NFCORE_STABLEEXPRESSION:STABLEEXPRESSION:ID_MAPPING:CLEAN_GENE_IDS (<accession>)` terminated with an error exit status (139) -- Execution is retried (1)
```

eventually leading to pipeline failure.

This is likely due to the python polars library not being compatible with macOS when run inside a container.

You should run the pipeline with `-profile micromamba` or `-profile conda`.

## Ǹo dataset found

For species that are not on Expression Atlas, the pipeline will not be able to find suitable datasets and will log the following message:

```
ERROR: Could not find any readily usable public dataset
...
```

> [!TIP]
> You can first try to have the pipeline fetch suitable datasets from NCBI GEO by providing the `--fetch_geo_accessions` flag.

In case no datasets are found, you'll have to find a way to get count datasets and to prepare them for the pipeline.
A good start is to check in the folder `<outdir>/public_data/geo/datasets/` if there are `rejected` subfolders. Such subfolders contain datasets that were downloaded (together with their experimental design) but failed to pass checks. Quite often, some of them be manually reprocessed to be suitable for the pipeline.

Finally, you may want to check by yourself on [NCBI GEO](https://www.ncbi.nlm.nih.gov/gds).

Alternatively, some public websites contain expression datasets that may be suitable for the pipeline, such as:

- [Bgee](https://www.bgee.org/)

## Not enough memory

The pipeline limits the number of downloaded datasets to a certain number in order to limit RAM usage, especially for `homo sapiens`.

However, on small computers, the limit may be too permissive and lead to RAM overhead. You can reduce the number of datasets downloaded by setting the `--random_sampling_size` to a lower value.

## Why do I get only a fraction of the public datasets available on Expression Atlas or NCBI GEO? Give them back!

To reduce the RAM overhead, the pipeline selects randomly a certain number of datasets, based on the number of samples they contain. To increase the number of collected datasets, you can increase the `--random_sampling_size` parameter.

[!TIP]

> A seed is also set in order to make the runs reproducible. You can change the subset of chosen datasets by changing the `--random_sampling_seed`.

## The pipeline failed to find a genome annotation for the specified species

If you know the length of the longest cDNA for each gene, you can provide gene lengths yourself with the `--gene_length` flag (see [Custom gene ID mapping / metadata / length](docs/usage.md#5-custom-gene-id-mapping-and-metadata)). In case you do not have access to gene length, TPM normalisation cannot be formed. A fallback is to use CPM normalisation by setting `--normalisation_method cpm`. It will introduce a small bias towards long genes, but this should not result in big changes.

## Java heap space

In some cases, in particular when running the pipeline on a very large number of datasets (such as for `Homo sapiens`), the Nextflow Java virtual machines can start to request a large amount of memory. You may happen to see the following error:

```
java.lang.OutOfMemoryError: Java heap space
```

We recommend to increase the memory available to Java:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

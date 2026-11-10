# nf-core/stableexpression: Guidelines for the selection of RT-qPCR reference gene candidates

Here are some guidelines that will help you find the best reference genes for your experiment.

## Specify a target gene

Using reference genes with different expression levels to those of your target genes can lead to important issues that are often overlooked. For example:
- erronate expression calculation (in case of non-linear efficiency of the primers across the full range of expression)
- signal saturation in digital PCR.

**nf-core/stableexpression** can help users identify stable reference gene candidates whose expression levels fall within a similar range to those of the target genes. The parameters `--target_genes` and `--target_gene_file` allow users to provide a list of target genes of interest. This will not affect the pipeline itself, but the MultiQC report will highlight the sections corresponding to these target genes, facilitating the subsequent selection of candidate reference genes.

Example (with `--target_genes PRUPE_4G144300`):

<p align="center">  
    <img title="MultiQC section with target gene" src="images/multiqc/multiqc_sections_with_target_gene.png">
</p>

> [!NOTE]
> To avoid any unrecognised target gene IDs, we recommend that you provide Ensembl gene IDs instead.

## Always test a set of candidate reference genes beforehand

Although the results of the `nf-core/stableexpression` pipeline should be highly reliable, users should consider the genes displayed as the most stable in their respective section as candidates, especially since the conditions of the samples used in the pipeline may differ slightly from those of the users. We recommend testing at least three different candidates in each section.

## Number of reference genes to use

The number of reference genes to use in RT-qPCR is a matter of debate. However, we recommend using at least three reference genes, as this will provide a more robust normalisation of your target genes.

Reference papers:
- [Bustin, S. A. et al. MIQE 2.0: Revision of the Minimum Information for Publication of Quantitative Real-Time PCR Experiments Guidelines. Clin Chem 71, 634–651 (2025).
](https://academic.oup.com/clinchem/article/71/6/634/8119148?login=false)
- [Vandesompele, J. et al. Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome Biol 3, research0034.1 (2002).](https://link.springer.com/article/10.1186/gb-2002-3-7-research0034)
- [Dowling, C. M., Walsh, D., Coffey, J. C. & Kiely, P. A. The importance of selecting the appropriate reference genes for quantitative real time PCR as illustrated using colon cancer cells and tissue. F1000Res 5, 99 (2016).
](https://f1000research.com/articles/5-99/v2)

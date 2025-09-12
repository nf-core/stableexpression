import polars as pl
import sys
from tqdm import tqdm
from statistics import mean

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"


count_file = sys.argv[1]
design_file = sys.argv[2]

count_lf = pl.scan_parquet(count_file)
design_df = pl.read_csv(design_file)

design_df = (
    design_df
    .with_columns(pl.concat_str([pl.col('batch'), pl.col('condition')], separator="_").alias("group"))
    .select("sample", "group")
)

group_to_sample_df = (
    design_df
    .group_by("group", maintain_order=True) # maintain order is better for repeatability and testing
    .agg("sample")
)

group_to_samples_dict = {
    d["group"]: d["sample"]
    for d in group_to_sample_df.to_dicts()
}
del group_to_sample_df

groups = list(group_to_samples_dict.keys())
n_groups = len(groups)

genes = count_lf.select(ENSEMBL_GENE_ID_COLNAME).collect().to_series().to_list()
k = len(genes)
renaming_dict = { f"column_{i}": gene for i, gene in enumerate(genes) }

if k <= 2:
    raise ValueError("Too few genes")

def get_overall_mean_for_group(df_with_means_over_samples):
    return df_with_means_over_samples.mean().item()


def get_means_over_samples(df):
    return (
        df
        .with_columns(
            mean_over_samples_for_gene=pl.concat_list(pl.all()).list.drop_nulls().list.mean()
        )
        .select("mean_over_samples_for_gene")
    )


def compute_minvar(df, target_gene) -> float:
    return (
        df
        .select((pl.col(target_gene) - pl.col(col)).alias(col) for col in genes if col != target_gene) # makes all pairwise differences with other genes
        .var(ddof=1) # computes variance
        .select(
            ( pl.concat_list(pl.all()).list.drop_nulls().list.min() / 4 ).alias("min") # get min of variances and divides by 4
        )
        .item()
    )


def get_unbiased_intragroup_variance(df, means_over_samples_df, group_overall_mean, samples):

    # TODO: see if it's correct
    # if only one sample in the group, there's no variance
    if len(samples) == 1:
        data = { gene: [0] for gene in genes }
        return pl.DataFrame(data)

    # lf is a lazyframe with a column being the gene ids (ensembl_gene_id)
    # and other columns being the samples
    # the current chunk corresponds to only one group
    # means_over_samples_df is a single column dataframe containing the means across each row (ie for each gene across samples)
    ng = len(samples)

    means_over_samples = means_over_samples_df.to_series().rename("mean_over_samples_for_gene")

    mean_over_genes = df.mean().transpose().to_series().rename("mean_over_genes_for_sample")

    sample_variance_df = (
        df
        .hstack([means_over_samples]) # adding column containing means over all samples in this group (for each gene)
        .select([
            (pl.col(c) - pl.col("mean_over_samples_for_gene")).alias(c) # y_igj - mean(y_ig*)
            for c in samples
        ])
        .transpose(include_header=True, column_names=genes) # columns are now genes
        .hstack([mean_over_genes]) # adding column containing means over all genes (for each sample)
        .select([
            ((pl.col(c) - pl.col("mean_over_genes_for_sample") + group_overall_mean) ** 2).alias(c) # r_igj ^2 = (y_igj - mean(y_ig*) -mean(y_*gj) + mean(y_*g*) ) ^ 2
            for c in genes
        ])
        .transpose(include_header=True, column_names=samples)
        .with_columns(
            sample_variance=pl.concat_list(samples).list.drop_nulls().list.sum() / ((ng - 1) * (1 - 2 / k)) # sum over j (samples) of r_igj ^2 terms
        )
        .select("sample_variance")
        .transpose()
        .rename(renaming_dict)
    )

    # sum of all sample variances for all genes
    sample_variance_sum_over_genes = sample_variance_df.select(pl.sum_horizontal(pl.all())).item() # sum of all s_ij² over all genes

    intra_var_df = (
        sample_variance_df
        .select([
            ( pl.col(c) - sample_variance_sum_over_genes / (k * (k -1)) ).alias(c)
            for c in genes
        ])
    )

    # if some values are negative, we need a special process
    genes_with_negative_values = (
        intra_var_df
        .select(col for col in genes if (intra_var_df[col] < 0).all()) # intra_var_df has only one row but it is a dataframe
        .columns
    )

    minvar = {}
    if genes_with_negative_values:
        transposed_df = df.transpose(include_header=True, column_names=genes).select(genes)
        for gene in genes_with_negative_values:
            minvar[gene] = compute_minvar(transposed_df, gene)

    return (
        intra_var_df
        .with_columns([
            pl.lit(val).alias(col) for col, val in minvar.items()
        ])
    )


def adjust_for_nb_of_samples_in_groups(unbiased_intragroup_variance_df, n_samples_list):
    return (
        unbiased_intragroup_variance_df
        .with_columns(
            n_samples=pl.Series(n_samples_list)
        )
        .select([
            (pl.col(c) / pl.col("n_samples")).alias(c)
            for c in genes
        ])
    )


def get_unbiased_intergroup_variance(gene_means_in_groups_df, dataset_overall_mean):
    mean_over_genes = gene_means_in_groups_df.mean().transpose().to_series().rename("mean_over_genes_for_group")

    return (
        gene_means_in_groups_df
        .with_columns(
            mean_over_groups_for_gene=pl.concat_list(pl.all()).list.drop_nulls().list.mean()
        )
        .select([
            (pl.col(c) - pl.col("mean_over_groups_for_gene")).alias(c)
            for c in gene_means_in_groups_df.columns
        ])
        .transpose(column_names=genes)
        .hstack([mean_over_genes])
        .select([
            (pl.col(c) - pl.col("mean_over_genes_for_group") + dataset_overall_mean).alias(c)
            for c in genes
        ])
        .select([ (pl.col(c) ** 2).alias(c) for c in genes ]) # square to get variance
    )


def compute_gamma_factor(diff_df, vardiff_df):
    first_term = (
        diff_df
        .with_columns(
            sum_of_squares=pl.concat_list(pl.all()).list.drop_nulls().list.sum() # sum over columns
        )
        .select("sum_of_squares")
        .sum() # sum over rows
        .select(
            ( pl.col("sum_of_squares") / ((n_groups - 1) * (k - 1)) ).alias("normalised_sum_of_squares")

        )
        .item()
    )

    second_term = (
        vardiff_df
        .with_columns(
            sum=pl.concat_list(pl.all()).list.drop_nulls().list.sum() # sum over columns
        )
        .select("sum")
        .sum() # sum over rows
        .select(
            (pl.col("sum") / (n_groups* k)).alias("normalised_sum")

        )
        .item()
    )

    return max(first_term - second_term, 0)


def apply_gamma_factor(gamma, diff_df, vardiff_df):
    difnew = diff_df * gamma / (gamma + vardiff_df)
    varnew = vardiff_df + gamma * vardiff_df / (gamma + vardiff_df)
    return difnew, varnew


def get_stability_values(unbiased_intragroup_variance_dfs, means_over_samples_dfs, group_overall_means):

    dataset_overall_mean = mean(group_overall_means)

    # putting together all intragroup variances
    unbiased_intragroup_variance_df = pl.concat(unbiased_intragroup_variance_dfs)

    group_mean_variance_df = adjust_for_nb_of_samples_in_groups(unbiased_intragroup_variance_df, n_samples_list)

    # putting together all means over samples for each group
    gene_means_in_groups_df = pl.concat(means_over_samples_dfs, how="horizontal")
    # adding mean over genes for each group (no need to compute it again)
    # gene_means_in_groups_df = pl.concat([gene_means_in_groups_df, group_overall_mean_df])

    intergroup_variance_df = get_unbiased_intergroup_variance(gene_means_in_groups_df, dataset_overall_mean)

    gamma = compute_gamma_factor(intergroup_variance_df, group_mean_variance_df)

    shrunk_intergroup_variance_df, shrunk_group_mean_variance_df = apply_gamma_factor(gamma, intergroup_variance_df, group_mean_variance_df)

    return (
        (
            shrunk_intergroup_variance_df.select([pl.col(c).abs() for c in genes])
            + shrunk_group_mean_variance_df.select([pl.col(c).sqrt() for c in genes])
        )
        .mean()
    )


unbiased_intragroup_variance_dfs = []
means_over_samples_dfs = []
group_overall_means = []
n_samples_list = []
cpt=0
for group, samples in tqdm(group_to_samples_dict.items()):
    cpt+=1
    if cpt > 30:
        break
    chunk_df = count_lf.select(samples).collect()
    means_over_samples_df = get_means_over_samples(chunk_df)
    group_overall_mean = get_overall_mean_for_group(means_over_samples_df)
    unbiased_intragroup_variance_df = get_unbiased_intragroup_variance(chunk_df, means_over_samples_df, group_overall_mean, samples)
    # storing intragroup values for each gene in this group
    unbiased_intragroup_variance_dfs.append(unbiased_intragroup_variance_df)
    # storing means over samples in this group for each gene
    means_over_samples_df = means_over_samples_df.rename({"mean_over_samples_for_gene": group})
    means_over_samples_dfs.append(means_over_samples_df)
    # storing overall mean of expression in this group, for all genes and samples
    group_overall_means.append(group_overall_mean)
    # storing nb of samples in this group
    n_samples_list.append(len(samples))






stab = get_stability_values(unbiased_intragroup_variance_dfs, means_over_samples_dfs, group_overall_means)
#print(stab)





"""
import time
s_squared = {}
overall_s_squared_factor = {}
for group, samples in tqdm(batch_condition_to_sample_dict.items()):
    overall_group_mean = batch_condition_mean_dict[group]
    nb_samples = len(samples)
    s_squared[group] = {}

    for gene in tqdm(genes):
        sample_ys = {
            k: v
            for k, v in lf.filter(pl.col(ENSEMBL_GENE_ID_COLNAME) == gene).select(samples).collect().to_dicts()[0].items()
            if v is not None
        }
        if not sample_ys:
            s_squared[group][gene] = None
            continue
        mean_over_samples_for_gene = mean(sample_ys.values())
        sample_r_squares = []
        for sample in sample_ys:
            r_square = sample_ys[sample] - mean_over_samples_for_gene - sample_mean_dict[sample] + overall_group_mean
            sample_r_squares.append(r_square)
        factor = nb_samples - 1 if nb_samples > 1 else 1
        s_squared[group][gene] = sum(sample_r_squares) / ( factor * (1 - 2 / nb_genes))

    # sum of s squared for all genes
    overall_s_squared_factor[group] = sum(s_squared[group].values())

sigma_squared = {}
for group, samples in batch_condition_to_sample_dict.items():
    sigma_squared[group] = {}
    s_squared_factor = overall_s_squared_factor[group]
    for gene in genes:
        if gene not in s_squared[group]:
            sigma_squared[group][gene] = None
            continue
        sigma_squared[group][gene] = s_squared[group][gene] - s_squared_factor / nb_genes * (nb_genes - 1)
"""



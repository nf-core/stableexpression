#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE = "gene_transcript_lengths.csv"

GFF_COLUMNS = [
    "chromosome",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]

DTYPES = {
    "chromosome": str,
    "source": str,
    "feature": str,
    "start": int,
    "end": int,
    "score": str,
    "strand": str,
    "phase": str,
    "attributes": str,
}


##################################################################
##################################################################
# FUNCTIONS
##################################################################
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Get CDNA lengths from GFF3 annotation file")
    parser.add_argument(
        "--annotation",
        type=Path,
        dest="annotation_file",
        required=True,
        help="Annotation file in GFF3 format",
    )
    return parser.parse_args()


def parse_gff3_file(annotation_file: Path):
    return pd.read_csv(
        annotation_file,
        sep="\t",
        names=GFF_COLUMNS,
        dtype=DTYPES,
        comment="#",
        on_bad_lines="warn",
    )


def compute_transcript_lengths(df: pd.DataFrame):
    exon_df = df.loc[df["feature"] == "exon"].copy()
    # extract transcript ID from attributes column for each exon
    exon_df["transcript_id"] = exon_df["attributes"].str.extract(
        r"Parent=transcript:([^;]+)"
    )
    # compute transcript length
    exon_df[config.CDNA_LENGTH_COLNAME] = exon_df["end"] - exon_df["start"] + 1
    exon_df = exon_df[["transcript_id", config.CDNA_LENGTH_COLNAME]]
    return exon_df.groupby("transcript_id", as_index=False).agg(
        {config.CDNA_LENGTH_COLNAME: "sum"}
    )


def compute_max_transcript_lengths_per_gene(
    df: pd.DataFrame, transcript_lengths_df: pd.DataFrame
):
    rna_cols = [
        feature
        for feature in df["feature"].unique()
        if "RNA" in feature and "gene" not in feature
    ]
    rna_df = df.loc[df["feature"].isin(rna_cols)].copy()

    # extract gene ID from attributes column for each transcript
    rna_df[config.GENE_ID_COLNAME] = rna_df["attributes"].str.extract(
        r"Parent=gene:([^;]+)"
    )
    # extract transcript ID from attributes column
    rna_df["transcript_id"] = rna_df["attributes"].str.extract(r"ID=transcript:([^;]+)")

    # merge with transcript lengths dataframe to get length
    merged_df = rna_df.merge(transcript_lengths_df, how="left", on="transcript_id")
    logger.info(
        f"Got length for {len(merged_df) / len(rna_df) * 100:.2f}% of transcripts"
    )
    # compute max transcript length per gene
    merged_df = merged_df[[config.GENE_ID_COLNAME, config.CDNA_LENGTH_COLNAME]]
    return merged_df.groupby(config.GENE_ID_COLNAME, as_index=False).agg(
        {config.CDNA_LENGTH_COLNAME: "max"}
    )


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################


def main():
    args = parse_args()

    logger.info("Parsing annotation file")
    df = parse_gff3_file(args.annotation_file)

    logger.info("Computing transcript lengths")
    transcript_lengths_df = compute_transcript_lengths(df)

    # keep only mRNA and exon features
    logger.info("Getting max transcript length per gene")
    gene_length_df = compute_max_transcript_lengths_per_gene(df, transcript_lengths_df)

    logger.info(f"Writing to {OUTFILE}")
    gene_length_df.to_csv(OUTFILE, index=False, header=True)


if __name__ == "__main__":
    main()

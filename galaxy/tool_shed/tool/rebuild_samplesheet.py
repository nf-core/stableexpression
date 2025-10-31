#!/usr/env/bin python
"""
Script dedicated to renaming files in the samplesheet provided.
In Galaxy, data files provided by users are given a new file name.
However, original file names can be retrieved from the name attribute of the file object (inside the tool XML file).
In this script, we replace the original name with the actual Galaxy path.

"""

import argparse
import logging
from pathlib import Path
import csv

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", dest="samplesheet", type=Path, required=True)
    parser.add_argument("--count-files", dest="count_files", type=str, required=True)
    parser.add_argument(
        "--count-filenames", dest="count_filenames", type=str, nargs="+", required=True
    )
    parser.add_argument("--design-files", dest="design_files", type=str, required=True)
    parser.add_argument(
        "--design-filenames",
        dest="design_filenames",
        type=str,
        nargs="+",
        required=True,
    )
    parser.add_argument("--out", dest="outfile", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # files and names arrive in the same order
    count_files = args.count_files.split(",")
    design_files = args.design_files.split(",")

    count_names_to_files = {
        name: file for file, name in zip(count_files, args.count_filenames)
    }
    design_names_to_files = {
        name: file for file, name in zip(design_files, args.design_filenames)
    }

    renamed_rows = []
    with open(args.samplesheet, "r", newline="") as fin:
        reader = csv.DictReader(fin)
        header = reader.fieldnames
        for row in reader:
            # getting original names (file names as written in the samplesheet)
            original_count_filename = Path(row["counts"]).name
            original_design_filename = Path(row["design"]).name
            # turning original names into new names (Galaxy file names)
            row["counts"] = count_names_to_files[original_count_filename]
            row["design"] = design_names_to_files[original_design_filename]
            renamed_rows.append(row)

    with open(args.outfile, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=header)

        writer.writeheader()
        for row in renamed_rows:
            writer.writerow(row)

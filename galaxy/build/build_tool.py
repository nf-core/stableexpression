#!/usr/bin/env python3

import logging
from pathlib import Path

from formatters import ConfigFormatter, SchemaFormatter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

tool_boilerplate_file = Path(__file__).parent / "static/template.xml"
tool_file = Path(__file__).parents[1] / "tool_shed/tool/nf_core_{}.xml"


def main():
    logger.info("Formatting config")
    # package_versions = ConfigFormatter.get_package_versions()
    pipeline_metadata = ConfigFormatter.get_pipeline_metadata()

    logger.info("Formatting schema")
    schema_formatter = SchemaFormatter()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # REPLACING ACTUAL PARAMS IN STATIC TOOL
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    with open(tool_boilerplate_file, "r") as fin:
        static_string = fin.read()

    # checking if package versions were filled by the user
    for package_version in ["OPENJDK_VERSION"]:
        if package_version in static_string:
            raise ValueError(
                f"You must fill the package version in place of {package_version} before building"
            )

    logger.info("Building tool XML file")
    tool_string = (
        static_string
        # .replace("NEXTFLOW_VERSION", package_versions["nextflow"])
        # .replace("APPTAINER_VERSION", package_versions["apptainer"])
        # .replace("OPENJDK_VERSION", package_versions["openjdk"])
        .replace("PIPELINE_VERSION", pipeline_metadata["version"])
        .replace("DESCRIPTION", schema_formatter.pipeline_description)
        .replace("PARAMETERS", schema_formatter.params_cli)
        .replace("INPUTS", schema_formatter.inputs)
    )

    pipeline_name = pipeline_metadata["name"].replace("nf-core/", "")
    outfile = Path(str(tool_file).format(pipeline_name))
    with open(outfile, "w") as fout:
        fout.write(tool_string)

    logger.info("Done")


if __name__ == "__main__":
    main()

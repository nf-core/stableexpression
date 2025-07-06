import logging

from formatters import SchemaFormatter, ConfigFormatter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

STATIC_TOOL_FILENAME = "static/nf_core_stableexpression.boilerplate.xml"
OUTPUT_TOOL_FILENAME = "../tool/nf_core_stableexpression.xml"


def main():
    logger.info("Formatting config")
    config_formatter = ConfigFormatter()

    logger.info("Formatting schema")
    schema_formatter = SchemaFormatter()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # REPLACING ACTUAL PARAMS IN STATIC TOOL
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    with open(STATIC_TOOL_FILENAME, "r") as fin:
        static_string = fin.read()

    logger.info("Building tool XML file")
    tool_string = (
        static_string.replace(
            "NEXTFLOW_VERSION", config_formatter.package_version["nextflow"]
        )
        .replace("APPTAINER_VERSION", config_formatter.package_version["apptainer"])
        .replace("OPENJDK_VERSION", config_formatter.package_version["openjdk"])
        .replace("PIPELINE_VERSION", config_formatter.pipeline_version)
        .replace("DESCRIPTION", schema_formatter.pipeline_description)
        .replace("PARAMETERS", schema_formatter.params_cli)
        .replace("INPUTS", schema_formatter.inputs)
    )

    with open(OUTPUT_TOOL_FILENAME, "w") as fout:
        fout.write(tool_string)

    logger.info("Done")


if __name__ == "__main__":
    main()

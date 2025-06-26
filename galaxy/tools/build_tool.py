import logging

from schema_formatter import SchemaFormatter
from config_formatter import ConfigFormatter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

STATIC_TOOL_FILENAME = "static_tool.xml"
OUTPUT_TOOL_FILENAME = "tool.xml"


def main():
    cformatter = ConfigFormatter()
    sformatter = SchemaFormatter()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # REPLACING ACTUAL PARAMS IN STATIC TOOL
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    with open(STATIC_TOOL_FILENAME, "r") as fin:
        static_string = fin.read()

    tool_string = (
        static_string.replace(
            "NEXTFLOW_VERSION", cformatter.package_version["nextflow"]
        )
        .replace("SINGULARITY_VERSION", cformatter.package_version["singularity"])
        .replace("PIPELINE_VERSION", cformatter.pipeline_version)
        .replace("DESCRIPTION", sformatter.pipeline_description)
        .replace("PARAMETERS", sformatter.params_cli)
        .replace("INPUTS", sformatter.inputs)
        .replace("USAGE_OPTIONS", sformatter.usage_options)
    )

    with open(OUTPUT_TOOL_FILENAME, "w") as fout:
        fout.write(tool_string)

    logger.info("Done")


if __name__ == "__main__":
    main()

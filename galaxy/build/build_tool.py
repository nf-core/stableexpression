import logging
from pathlib import Path

from formatters import SchemaFormatter, ConfigFormatter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

tool_boilerplate_file = Path(__file__).parent / "static/boilerplate.xml"
tool_file = Path(__file__).parents[1] / "tool/nf_core_{}.xml"


def main():
    logger.info("Formatting config")
    package_versions = ConfigFormatter.get_package_versions()
    pipeline_metadata = ConfigFormatter.get_pipeline_metadata()

    logger.info("Formatting schema")
    schema_formatter = SchemaFormatter()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # REPLACING ACTUAL PARAMS IN STATIC TOOL
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    with open(tool_boilerplate_file, "r") as fin:
        static_string = fin.read()

    logger.info("Building tool XML file")
    tool_string = (
        static_string.replace("NEXTFLOW_VERSION", package_versions["nextflow"])
        .replace("APPTAINER_VERSION", package_versions["apptainer"])
        .replace("OPENJDK_VERSION", package_versions["openjdk"])
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

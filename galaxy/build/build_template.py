#!/usr/bin/env python3

import logging
from pathlib import Path

from formatters import ConfigFormatter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

STATIC_TOOL_FILENAME = Path(__file__).parent / "static/template.xml"
BOILERPLATE_FILENAME = Path(__file__).parent / "static/template.boilerplate.xml"


def main():
    logger.info("Parsing config")
    pipeline_metadata = ConfigFormatter.get_pipeline_metadata()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # REPLACING ACTUAL PARAMS IN STATIC TOOL
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    with open(BOILERPLATE_FILENAME, "r") as fin:
        boilerplate_string = fin.read()

    pipeline_name = pipeline_metadata["name"].replace("nf-core/", "")

    logger.info("Building template XML file")
    template_string = boilerplate_string.replace("PIPELINE_NAME", pipeline_name)

    with open(STATIC_TOOL_FILENAME, "w") as fout:
        fout.write(template_string)

    logger.info("Done")


if __name__ == "__main__":
    main()

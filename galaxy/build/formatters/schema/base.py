from pathlib import Path
import json
from dataclasses import dataclass, field
from typing import ClassVar
from . import parameter


@dataclass
class SchemaFormatter:
    SCHEMA_FILE: ClassVar[Path] = Path(__file__).parents[4] / "nextflow_schema.json"
    PARAMS_TO_IGNORE: ClassVar[list] = ["outdir", "email", "multiqc_title"]
    SECTIONS_TO_IGNORE: ClassVar[list] = [
        "institutional_config_options",
        "generic_options",
    ]
    SECTIONS_TO_EXPAND: ClassVar[list] = ["input_output_options"]

    pipeline_description: str = field(init=False)
    inputs: str = field(init=False)
    params_cli: str = field(init=False)
    _pipeline_params: dict = field(init=False)

    _inputs: list = field(init=False, default_factory=list)
    _params_cli: list = field(init=False, default_factory=list)

    def __post_init__(self):
        self.parse_schema_file()

    def parse_schema_file(self):
        with open(self.SCHEMA_FILE, "r") as f:
            pipeline_schema = json.load(f)

        self.pipeline_description = pipeline_schema["description"].strip("\n")
        self._pipeline_params = pipeline_schema["$defs"]

        # PARSING PARAMETERS AND BUILDING STRINGS
        for section, section_dict in self._pipeline_params.items():
            if section in self.SECTIONS_TO_IGNORE:
                continue

            section_inputs, section_params_cli, section_usage_options = (
                self.format_input_section(section, section_dict)
            )

            self._inputs += section_inputs
            self._params_cli += section_params_cli

        self.inputs = "\n".join(self._inputs)
        self.params_cli = "\n".join(self._params_cli)

    def format_input_section(
        self, section: str, section_dict: dict
    ) -> tuple[list, list, list]:
        section_inputs = []
        section_params_cli = []
        section_usage_options = []

        section_title = ""
        section_help = ""

        if title := section_dict.get("title"):
            section_title = f' title="{title}"'
        if description := section_dict.get("description"):
            section_help = f' help="{description}"'

        section_expanded = (
            ' expanded="true"'
            if section in self.SECTIONS_TO_EXPAND
            else ' expanded="false"'
        )

        section_inputs.append(
            f'\t\t<section name="{section}"{section_title}{section_help}{section_expanded}>'
        )
        section_usage_options.append("\n\t" + section.capitalize().replace("_", " "))

        required_params = section_dict.get("required", [])

        for param, param_dict in section_dict["properties"].items():
            if param not in self.PARAMS_TO_IGNORE:
                optional = param not in required_params

                # checking if param must be parsed in a generic or in a custom way
                if param in parameter.PARAMETER_TO_CUSTOM_CLASS:
                    class_ = parameter.PARAMETER_TO_CUSTOM_CLASS[param]
                else:
                    class_ = parameter.BaseParameterFormatter

                param_formatter = class_(param, section, param_dict, optional)

                # input arguments
                section_inputs.append(param_formatter.get_input())
                # cli
                section_params_cli.append(param_formatter.get_cli())

        section_inputs.append("\t\t</section>")

        return section_inputs, section_params_cli, section_usage_options

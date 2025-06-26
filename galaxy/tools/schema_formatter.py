from pathlib import Path
import json
from dataclasses import dataclass, field
from typing import ClassVar


@dataclass
class BaseSchemaFormatter:
    SCHEMA_FILE: ClassVar[Path] = Path(__file__).parents[2] / "nextflow_schema.json"
    PARAMS_TO_IGNORE: ClassVar[list] = ["outdir", "email", "multiqc_title"]
    SECTIONS_TO_IGNORE: ClassVar[list] = [
        "institutional_config_options",
        "generic_options",
    ]
    SECTIONS_TO_EXPAND: ClassVar[list] = ["input_output_options"]
    NF_TYPES_TO_GALAXY: ClassVar[dict] = {
        "string": "text",
        "boolean": "boolean",
        "integer": "integer",
        "number": "float",
    }

    pipeline_description: str = field(init=False)
    inputs: str = field(init=False)
    params_cli: str = field(init=False)
    usage_options: str = field(init=False)
    _pipeline_params: dict = field(init=False)

    _inputs: list = field(init=False, default_factory=list)
    _params_cli: list = field(init=False, default_factory=list)
    _usage_options: list = field(init=False, default_factory=list)

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
            self._usage_options += section_usage_options

        self.inputs = "\n".join(self._inputs)
        self.params_cli = "\n".join(self._params_cli)
        self.usage_options = "\n".join(self._usage_options)

    def format_input_param(self, param: str, param_dict: dict, optional: bool) -> str:
        """
        building input param
        """

        input_param_str = '\t\t\t<param name="{param}" type="{type}" {label}{format}{value}{min}{max}{true_false}{help}{optional} />'
        param_format = ""
        param_label = ""
        param_help = ""
        param_true_false = ""
        param_value = ""
        param_min = ""
        param_max = ""
        param_optional = ' optional="true"' if optional else ' optional="false"'

        param_type = param_dict["type"]
        default_value = param_dict.get("default")

        if param_type == "string" and param_dict.get("format") == "file-path":
            input_type = "data"
            if pattern := param_dict.get("pattern"):
                # TODO: handle multiple extensions
                extension = pattern.split(".")[-1].strip("$")
                param_format = f' format="{extension}"'

            # if param is an optional file with multiple possible values, it requires special handling
            # see https://docs.galaxyproject.org/en/latest/dev/schema.html#id51

        else:
            input_type = self.NF_TYPES_TO_GALAXY[param_type]

            if param_type == "boolean":
                param_true_false = f' truevalue="--{param}" falsevalue=""'

            elif param_type in ["integer", "number"]:
                if minimum := param_dict.get("minimum"):
                    param_min = f' min="{minimum}"'
                if maximum := param_dict.get("maximum"):
                    param_max = f' max="{maximum}"'

        # handle parameter with enum (options)
        if options := param_dict.get("enum"):
            input_type = "select"
            input_param_str = input_param_str.replace(" />", ">\n")
            base_option = (
                '\t\t\t<option value="{option}"{selected_arg}>{label}</option>\n'
            )

            for option in options:
                selected_arg = ' selected="true"' if option == default_value else ""
                option_param = base_option.format(
                    option=option,
                    label=option.capitalize(),
                    selected_arg=selected_arg,
                )
                input_param_str += "\t" + option_param

            input_param_str += "\t\t\t</param>"

        else:
            if default_value:
                param_value = f' value="{default_value}"'

        if description := param_dict.get("description"):
            param_label = f' label="{description}"'
        if help_text := param_dict.get("help_text"):
            param_help = f' help="{help_text}"'

        return input_param_str.format(
            param=param,
            type=input_type,
            label=param_label,
            format=param_format,
            value=param_value,
            min=param_min,
            max=param_max,
            true_false=param_true_false,
            help=param_help,
            optional=param_optional,
        )

    @staticmethod
    def format_input_param_cli(param: str, section: str, optional: bool) -> str:
        if optional:
            return f"\t\t\t#if {section}.{param}\n\t\t\t  --{param} {section}.{param}\n\t\t\t#end if"
        else:
            return f"--{param} {section}.{param}"

    @staticmethod
    def format_input_param_usage(param: str, param_dict: dict, optional: bool) -> str:
        required_param = "" if optional else "[REQUIRED] "
        return f'\t\t--{param} <{param_dict["type"]}> {required_param}: {param_dict["description"]}'

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
                # input arguments
                input_param = self.format_input_param(param, param_dict, optional)
                section_inputs.append(input_param)
                # cli
                param_cli = self.format_input_param_cli(param, section, optional)
                section_params_cli.append(param_cli)
                # usage (help)
                input_param_usage = self.format_input_param_usage(
                    param, param_dict, optional
                )
                section_usage_options.append(input_param_usage)

        section_inputs.append("\t\t</section>")

        return section_inputs, section_params_cli, section_usage_options


class SchemaFormatter(BaseSchemaFormatter):
    pass

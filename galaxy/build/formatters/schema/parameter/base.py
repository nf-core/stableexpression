from dataclasses import dataclass
from typing import ClassVar, override


@dataclass
class Validator:
    """ """

    PATTERN: ClassVar[str] = (
        '\t\t\t<validator type="{type}" message="{message}">{expression}</validator>\n'
    )

    type: str
    message: str
    expression: str

    @override
    def __str__(self):
        return self.PATTERN.format(
            type=self.type, message=self.message, expression=self.expression
        )


@dataclass
class Option:
    """
    Represents an option for a parameter.

    Attributes:
        value (str): The value of the option.
        default_value (str): The default value of the option.
        optional (bool): Whether the option is optional.
    """

    PATTERN: ClassVar[str] = (
        '\t\t\t<option value="{option}"{selected_arg}>{label}</option>\n'
    )

    value: str
    default_value: str | None
    optional: bool

    @override
    def __str__(self):
        selected_arg = ' selected="true"' if self.value == self.default_value else ""
        return self.PATTERN.format(
            option=self.value, label=self.value.capitalize(), selected_arg=selected_arg
        )


@dataclass
class BaseParameterFormatter:
    NF_TYPES_TO_GALAXY: ClassVar[dict] = {
        "string": "text",
        "boolean": "boolean",
        "integer": "integer",
        "number": "float",
    }
    BASE_INPUT_PARAM: ClassVar[str] = (
        '\t\t\t<param name="{param}" type="{type}" {label}{format}{value}{min}{max}{true_false}{help}{optional} />'
    )

    param: str
    section: str
    param_dict: dict
    optional: bool

    @staticmethod
    def enrich_input_param(input_param_str: str, args: list[str]) -> str:
        # opening param for enrichment
        input_param_str = input_param_str.replace(" />", ">\n")
        # adding each arg in a separate line
        for arg in args:
            input_param_str += "\t" + arg
        # closing
        input_param_str += "\t\t\t</param>"
        return input_param_str

    @staticmethod
    def extract_extensions(extension_str: str):
        def clean_extension(ext: str) -> str:
            ext = ext.strip().lower()
            if ext == "yml":
                return "yaml"
            return ext

        # removing the .dat extension, that is only used in the pipeline
        # in order to allow files from the Galaxy file system (all renamed in .dat)
        base_extensions = [ext for ext in extension_str.split("|") if ext != "dat"]
        # Galaxy does not allow 'yml', only 'yaml'
        return list(set([clean_extension(ext) for ext in base_extensions]))

    def process_file_param(self):
        input_type = "data"
        # removing extension check as files are renamed in <hash>.dat files by Galaxy
        if pattern := self.param_dict.get(
            "pattern"
        ):  # going from something like "^\\S+\\.(csv|yaml)$" to "csv,ya
            # getting the extensions part
            extension_str = pattern.split(".")[-1]
            # removes recursively all leading and traling "(", ")" and "$"
            extension_str = extension_str.strip("$()")
            # getting list of extensions; removing dat because this extension is specifically made to handle Galaxy filename
            formated_extensions_str = ",".join(self.extract_extensions(extension_str))
            param_format = f' format="{formated_extensions_str}"'
        else:
            # there is no specific pattern provided in the schema, this means that the format does not matter much
            # however, the planemo linter needs a format, so we specify format="data"
            param_format = ' format="data"'
        return input_type, param_format

    def get_input(self) -> str:
        """
        building input param
        """

        # making copy of base input param string
        input_param_str = self.BASE_INPUT_PARAM

        param_format = ""
        param_label = ""
        param_help = ""
        param_true_false = ""
        param_value = ""
        param_min = ""
        param_max = ""
        param_optional = ' optional="true"' if self.optional else ' optional="false"'

        param_type = self.param_dict["type"]
        default_value = self.param_dict.get("default")

        # special case when parameter is a file
        if param_type == "string" and self.param_dict.get("format") == "file-path":
            input_type, param_format = self.process_file_param()

        # all other types
        else:
            input_type = self.NF_TYPES_TO_GALAXY[param_type]

            if param_type == "boolean":
                param_true_false = f' truevalue="--{self.param}" falsevalue=""'

            elif param_type in ["integer", "number"]:
                if minimum := self.param_dict.get("minimum"):
                    param_min = f' min="{minimum}"'
                if maximum := self.param_dict.get("maximum"):
                    param_max = f' max="{maximum}"'

            elif param_type == "string":
                # if there is a pattern for this string, we need to enrich this XML section with a validator
                # TODO: handle (rare) case where bot enum and pattern are given
                if pattern := self.param_dict.get("pattern"):  # regex
                    msg = f"must match regular expression {pattern}"
                    validator = Validator(type="regex", message=msg, expression=pattern)
                    input_param_str = self.enrich_input_param(
                        input_param_str, args=[str(validator)]
                    )

        # handle parameter with enum (options)
        if option_values := self.param_dict.get("enum"):
            input_type = "select"
            options = [
                Option(value, default_value, self.optional) for value in option_values
            ]
            input_param_str = self.enrich_input_param(
                input_param_str, args=[str(option) for option in options]
            )

        else:
            if default_value is not None:
                param_value = f' value="{default_value}"'

        if description := self.param_dict.get("description"):
            param_label = f'label="{description}"'
        if help_text := self.param_dict.get("help_text"):
            param_help = f' help="{help_text}"'

        return input_param_str.format(
            param=self.param,
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

    def get_cli(self) -> str:
        # extra quotes if string parameter
        value = (
            f'"${self.section}.{self.param}"'
            if self.param_dict["type"] == "string"
            else f"${self.section}.{self.param}"
        )
        if self.optional:
            return f"\t\t\t#if ${self.section}.{self.param}\n\t\t\t  --{self.param} {value}\n\t\t\t#end if"
        else:
            return f"\t\t\t--{self.param} {value}"

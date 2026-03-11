import re
from dataclasses import dataclass
from typing import override

from .base import BaseParameterFormatter


@dataclass
class DatasetsParameterFormatter(BaseParameterFormatter):
    # if param is an optional file with multiple possible values, it requires special handling
    # see https://docs.galaxyproject.org/en/latest/dev/schema.html#id51

    @override
    def get_input(self) -> str:
        input_param_str = super().get_input()
        # setting to required
        # changing param name
        input_param_str = input_param_str.replace(
            'optional="true"', 'optional="false"'
        ).replace(self.param, "samplesheet")
        # changing label
        input_param_str = re.sub(
            r'label="[\s\w]*"', 'label="Samplesheet"', input_param_str
        )

        # adding conditional statement
        return f""" \t\t\t<conditional name="datasets">
                <param name="provide_datasets" type="select" label="Provide custom count datasets?" >
                    <option value="true">Yes</option>
                    <option selected="true" value="false">No</option>
                </param>
                <when value="true">
        {input_param_str}
                    <param name="count_datasets" label="Count datasets" type="data" format="csv" multiple="true" optional="false" help="User count datasets in CSV format" />
                    <param name="experimental_designs" label="Experimental designs" type="data" format="csv" multiple="true" optional="true" help="Experimental designs relative to the provided count datasets" />
                </when>
                <when value="false">
                </when>
            </conditional>"""

    @override
    def get_cli(self) -> str:
        # see https://planemo.readthedocs.io/en/latest/writing_advanced.html#consuming-collections
        return f"""
        \t#if ${self.section}.datasets.provide_datasets == "true":
        \t\t--datasets renamed_samplesheet.csv
        \t#end if"""

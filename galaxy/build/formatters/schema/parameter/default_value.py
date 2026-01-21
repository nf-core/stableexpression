from dataclasses import dataclass
from .base import BaseParameterFormatter


@dataclass
class DefaultValueParameterFormatter(BaseParameterFormatter):
    def __post_init__(self):
        self.param_dict["default"] = "Solanum tuberosum"

from dataclasses import dataclass
from .base import BaseParameterFormatter


@dataclass
class RequiredParameterFormatter(BaseParameterFormatter):
    def __post_init__(self):
        self.optional = False

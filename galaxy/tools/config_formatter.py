from pathlib import Path
import subprocess
import re
from dataclasses import dataclass, field
from typing import ClassVar


@dataclass
class BaseConfigFormatter:
    CONFIG_FILE: ClassVar[Path] = Path(__file__).parents[2] / "nextflow.config"
    MAIN_FILE: ClassVar[Path] = Path(__file__).parents[2] / "main.nf"
    NXF_VERSION_COMMAND_TEMPLATE: ClassVar[str] = (
        "micromamba search --override-channels --channel bioconda 'pkg' | grep pkg | awk '{print $2}' | sort | tail -1"
    )
    PACKAGES: ClassVar[list] = ["nextflow", "singularity"]

    pipeline_version: str = field(init=False)
    package_version: dict = field(init=False, default_factory=dict)

    def __post_init__(self):
        # CONDA PACKAGE VERSIONS
        for package in self.PACKAGES:
            self.package_version[package] = self.get_package_version(package)

        # PARSING CONFIG
        with open(self.CONFIG_FILE, "r") as f:
            pipeline_config = f.read()

        self.pipeline_version = self.get_pipeline_version(pipeline_config)

    @classmethod
    def get_package_version(cls, package_name: str) -> str:
        """
        Get latest conda version of package
        """
        nxf_version_command = cls.NXF_VERSION_COMMAND_TEMPLATE.replace(
            "pkg", package_name
        )
        result = subprocess.run(
            nxf_version_command, shell=True, capture_output=True, text=True, check=True
        )

        if result.stderr:
            raise RuntimeError(f"Command error: {result.stderr}")

        return result.stdout.strip("\n")

    @staticmethod
    def get_pipeline_version(pipeline_config: str):
        # Regular expression to find the manifest block and extract the version
        manifest_pattern = re.compile(r"manifest\s*{\s*(.*?)\s*}", re.DOTALL)
        manifest_match = manifest_pattern.search(pipeline_config)
        version = None

        if manifest_match:
            manifest_content = manifest_match.group(1)
            # Regular expression to find the version field
            version_pattern = re.compile(r'version\s*=\s*[\'"](.*?)[\'"]')
            version_match = version_pattern.search(manifest_content)

            if version_match:
                version = version_match.group(1)

        if version is None:
            raise ValueError("No version found in pipeline config")

        return version


@dataclass
class ConfigFormatter(BaseConfigFormatter):
    pass

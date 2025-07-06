from pathlib import Path
import requests
import re
from dataclasses import dataclass, field
from typing import ClassVar
from packaging.version import parse as vparse
import logging

logger = logging.getLogger(__name__)


@dataclass
class BaseConfigFormatter:
    CONFIG_FILE: ClassVar[Path] = Path(__file__).parents[4] / "nextflow.config"
    MAIN_FILE: ClassVar[Path] = Path(__file__).parents[4] / "main.nf"
    PACKAGES_REPOS: ClassVar[dict] = {
        "nextflow": "bioconda",
        "apptainer": "conda-forge",
        "openjdk": "conda-forge",
    }

    pipeline_version: str = field(init=False)
    package_version: dict = field(init=False, default_factory=dict)
    executable: str = field(init=False)

    def __post_init__(self):
        # CONDA PACKAGE VERSIONS
        for package, repo in self.PACKAGES_REPOS.items():
            self.package_version[package] = self.get_package_version(package, repo)

        # PARSING CONFIG
        with open(self.CONFIG_FILE, "r") as f:
            pipeline_config = f.read()

        self.pipeline_version = self.get_pipeline_version(pipeline_config)

    @classmethod
    def get_package_version(cls, package: str, repo: str) -> str:
        """
        Get latest pip version of package
        """
        logger.info(f"Getting latest version of package {package}")
        url = f" https://api.anaconda.org/package/{repo}/{package}"
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            versions = sorted(
                data["versions"], reverse=True, key=vparse
            )  # from latest to oldest
            return versions[0]  # most recent
        except requests.RequestException as e:
            raise RuntimeError(f"Error fetching version info: {e}")

    @staticmethod
    def get_pipeline_version(pipeline_config: str):
        # regular expression to find the manifest block and extract the version
        manifest_pattern = re.compile(r"manifest\s*{\s*(.*?)\s*}", re.DOTALL)
        manifest_match = manifest_pattern.search(pipeline_config)
        version = None

        if manifest_match:
            manifest_content = manifest_match.group(1)
            # regular expression to find the version field
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

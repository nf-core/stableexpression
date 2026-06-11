import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

import requests
from packaging.version import parse as vparse

logger = logging.getLogger(__name__)


@dataclass
class BaseConfigFormatter:
    """
    Base class for extracting metadata from the pipeline's config files.
    """

    CONFIG_FILE: ClassVar[Path] = Path(__file__).parents[4] / "nextflow.config"
    MAIN_FILE: ClassVar[Path] = Path(__file__).parents[4] / "main.nf"
    PACKAGES_REPOS: ClassVar[dict] = {
        "nextflow": "bioconda",
        "micromamba": "conda-forge",
        "openjdk": "conda-forge",
    }

    @classmethod
    def get_package_versions(cls) -> dict:
        # CONDA PACKAGE VERSIONS
        package_version = {}
        for package, repo in cls.PACKAGES_REPOS.items():
            package_version[package] = cls.get_package_version(package, repo)
        return package_version

    @staticmethod
    def get_package_version(package: str, repo: str) -> str:
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

    @classmethod
    def get_pipeline_metadata(cls) -> dict:
        #  PARSING CONFIG
        with open(cls.CONFIG_FILE, "r") as f:
            pipeline_config = f.read()

        # regular expression to find the manifest block and extract the version
        manifest_pattern = re.compile(r"manifest\s*{\s*(.*?)\s*}", re.DOTALL)
        manifest_match = manifest_pattern.search(pipeline_config)
        version = None
        name = None

        if manifest_match:
            manifest_content = manifest_match.group(1)

            # regular expression to find the version field
            name_pattern = re.compile(r'name\s*=\s*[\'"](.*?)[\'"]')
            name_match = name_pattern.search(manifest_content)
            if name_match:
                name = name_match.group(1)
            else:
                raise ValueError("No name found in pipeline config")

            # regular expression to find the version field
            version_pattern = re.compile(r'version\s*=\s*[\'"](.*?)[\'"]')
            version_match = version_pattern.search(manifest_content)
            if version_match:
                version = version_match.group(1)
            else:
                raise ValueError("No version found in pipeline config")

        return dict(name=name, version=version)


@dataclass
class ConfigFormatter(BaseConfigFormatter):
    pass

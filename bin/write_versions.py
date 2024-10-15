#!/usr/bin/env python3
import platform
import argparse
from pathlib import Path

VERSION_FILENAME = 'versions.yml'
CONDA_ENVIRONMENT_FILENAME = 'environment.yml'
PIP_ENVIRONMENT_FILENAME = 'requirements.txt'


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('module_dir', type=str)
    parser.add_argument('--no-python', action='store_true')
    return parser.parse_args()



def get_installed_conda_packages(module_dir: str):
    """
    Return a dictionary of packages installed with conda in the environment.

    The dictionary returned is {package_name: version}.
    """

    packages = {}
    in_dependencies = False

    env_file = Path(module_dir) / CONDA_ENVIRONMENT_FILENAME

    with open(env_file, 'r') as fin:
        for line in fin.readlines():

            line = line.strip()

            if line.startswith('dependencies:'):
                in_dependencies = True
                continue

            elif in_dependencies:

                if not line.startswith('-'):
                    break
                if line.startswith('- pip:'): # assuming that pip packages are at the end
                    break

                lst = line.split('::')[1].split('==') # puts package and version in a 2-elements list

                try:
                    package, version = lst[0], lst[1]
                except IndexError:
                    raise IndexError(f'No version was found at line {line} in {str(env_file)}')

                packages[package] = version

    return packages


def get_installed_pip_packages(module_dir: str):
    """
    Return a dictionary of packages installed with pip in the environment.

    The dictionary returned is {package_name: version}.
    """

    packages = {}

    env_file = Path(module_dir) / PIP_ENVIRONMENT_FILENAME

    with open(env_file, 'r') as fin:
        for line in fin.readlines():

            line = line.strip()
            lst = line.split('==') # puts package and version in a 2-elements list

            try:
                package, version = lst[0], lst[1]
            except IndexError:
                raise IndexError(f'No version was found at line {line} in {str(env_file)}')

            packages[package] = version

    return packages


def pip_requirements_in_env_file(module_dir: str):

    env_file = Path(module_dir) / PIP_ENVIRONMENT_FILENAME

    pip_section = False

    with open(env_file, 'r') as file:
        for line in file:
            # Strip leading/trailing spaces and handle indentation
            stripped_line = line.strip()

            # Check if we're inside the pip section
            if stripped_line == "- pip:":
                pip_section = True

            elif pip_section:
                # If we're in the pip section, look for the requirement
                if stripped_line == "- -r requirements.txt":
                    return True
                # If the next dependency section starts, exit the pip section
                elif stripped_line.startswith("-"):
                    pip_section = False

    return False


def write_dict_to_yaml(file_path: str, data: dict):
    with open(file_path, 'w') as fout:
        for key, value in data.items():
            if isinstance(value, dict):
                # Write the key for the dictionary and call recursively
                fout.write(f"{key}:\n")
                for subkey, subvalue in value.items():
                    fout.write(f"  {subkey}: {subvalue}\n")
            else:
                # Write simple key-value pairs
                fout.write(f"{key}: {value}\n")


def main():

    args = parse_args()

    packages = {'python': platform.python_version()}

    installed_packages = get_installed_conda_packages(args.module_dir)
    if pip_requirements_in_env_file(args.module_dir):
        pip_packages = get_installed_pip_packages(args.module_dir)
        installed_packages.update(pip_packages)

    packages.update(installed_packages)
    data = {'"${task.process}"': packages}
    write_dict_to_yaml(VERSION_FILENAME, data)


if __name__ == '__main__':
    main()

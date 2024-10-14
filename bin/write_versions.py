#!/usr/bin/env python3
import sys
import argparse
from pathlib import Path

VERSION_FILENAME = 'versions.yml'
ENVIRONMENT_FILENAME = 'environment.yml'


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('module_dir', type=str)
    parser.add_argument('--no-python', action='store_true')
    return parser.parse_args()



def get_installed_packages(module_dir: str):
    """
    Return a dictionary of packages installed in the environment.

    The dictionary returned is {package_name: version}.
    """

    packages = {}
    in_dependencies = False

    env_file = Path(module_dir) / ENVIRONMENT_FILENAME

    with open(env_file, 'r') as fin:
        for line in fin.readlines():

            line = line.strip()

            if line.startswith('dependencies:'):
                in_dependencies = True
                continue

            elif in_dependencies:
                if not line.startswith('-'):
                    break
                lst = line.split('::')[1].split('==') # puts package and version in a 2-elements list

                try:
                    package, version = lst[0], lst[1]
                except IndexError:
                    raise IndexError(f'No version was found at line {line} in {str(env_file)}')

                packages[package] = version

    return packages


def main():

    args = parse_args()

    lines = ['"${task.process}"']
    if not args.no_python:
        lines.append(f'    python: {sys.version}')

    packages = get_installed_packages(args.module_dir)
    for package, version in packages.items():
        lines.append(f'    {package}: {version}')

    with open(VERSION_FILENAME, 'w') as fout:
        fout.write('\\n'.join(lines))


if __name__ == '__main__':
    main()

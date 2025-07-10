#!/usr/bin/env bash

galaxy_dir="$(dirname $(dirname $(readlink -f "$0")))"
tool_dir="${galaxy_dir}/tool"

planemo serve $tool_dir


#!/bin/bash

galaxy_dir="$(dirname $(dirname $(readlink -f "$0")))"
tool_file="${galaxy_dir}/tool/nf_core_stableexpression.xml"

planemo lint $tool_file --fail_level error

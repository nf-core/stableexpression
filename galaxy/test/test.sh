#!/usr/bin/env bash

galaxy_dir="$(dirname $(dirname $(readlink -f "$0")))"
tool_file="${galaxy_dir}/tool/nf_core_stableexpression.xml"

# add --update_test_data to create output file
planemo test $tool_file


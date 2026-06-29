#!/usr/bin/env bash

set -euo pipefail

PREP_DIR="$(dirname $(realpath "$0"))"

echo "Checking consistency of outputs"
${PREP_DIR}/check_consistency.sh

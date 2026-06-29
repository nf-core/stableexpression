#!/usr/bin/env bash

set -euo pipefail

PREP_DIR="$(dirname $(realpath "$0"))"
REPO_ROOT="$(dirname $(dirname $PREP_DIR))"

MAIN_SCRIPT="${REPO_ROOT}/main.nf"
PREP_CONFIG="${PREP_DIR}/.config"
BASE_OUTDIR="${REPO_ROOT}/results/consistency_test"


NB_RUNS=4

rm -rf ${BASE_OUTDIR}
for i in $(seq 1 $NB_RUNS)
do
    outdir=${BASE_OUTDIR}/run_${i}
    echo "Running test ${i} with output directory ${outdir}"
    nextflow run ${MAIN_SCRIPT} -profile apptainer,test -c ${PREP_CONFIG} --outdir "$outdir"
    if [ $i -eq 1 ]; then
        reference_outdir="$outdir"
    else
        res=$(diff --recursive --brief --exclude='pipeline_info' --exclude='multiqc' $outdir $reference_outdir | awk -F' ' '{print $3}')
        if [ $res ]; then
            echo "The following files differ: $res"
            exit 1
        else
            echo "Run $i is consistent with reference run 1"
        fi
    fi
done

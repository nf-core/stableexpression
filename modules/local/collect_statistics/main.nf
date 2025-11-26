process COLLECT_STATISTICS {

    tag "${file.baseName}"
    label "process_high_memory"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/60/604657081a64b39e17bb6ad307e545aa6aebf4133b64d6766515c9789bb2d304/data':
        'community.wave.seqera.io/library/pandas_tqdm:2ca37c1047243549' }"

    input:
    path file

    output:
    path '*.transposed.csv',                                                                                          emit: csv
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions

    script:
    """
    collect_statistics.py $file
    """

}

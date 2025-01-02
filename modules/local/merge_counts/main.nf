process MERGE_COUNTS {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/734bb282c77cae829345845f1280e8298e1d839c6c0411b9826c771c41d45d99/data':
        'community.wave.seqera.io/library/pandas_python:0ae9367e11611863' }"

    input:
    path(f1, stageAs: "?/*")
    path(f2, stageAs: "?/*")


    output:
    path("merged.csv"),                                                                                               emit: merged
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions


    script:
    """
    merge_count_files.py --f1 $f1 --f2 $f2
    """


    stub:
    """
    touch merged.csv
    """

}

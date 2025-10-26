process MERGE_DESIGNS {

    label 'process_high'

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f1/f1c30725ef181337de8749d5b54eacb1a8e1f97ac5e43fe15ec34a61789a7320/data':
        'community.wave.seqera.io/library/pandas:2.3.2--baef3004955c4a32' }"

    input:
    path design_files, stageAs: "?/*"

    output:
    path 'whole_design.csv',                                                                                          emit: design
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions

    script:
    """
    merge_designs.py \\
        --designs "$design_files"
    """

}

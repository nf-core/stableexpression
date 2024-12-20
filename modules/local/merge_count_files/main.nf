process MERGE_COUNT_FILES {

    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e87c143b4e7b31e1d5db5518d5c3d0e82fe20c4a9607e668e3fc8b390257d4f7/data':
        'community.wave.seqera.io/library/pandas:2.2.3--9b034ee33172d809' }"

    input:
    path csv_files

    output:
    path 'all_counts.csv', emit: csv
    tuple val("${task.process}"), val('python'),   eval('python3 --version'),                                         topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    merge_count_files.py $csv_files --outfile all_counts.csv
    """

    stub:
    """
    touch all_counts.csv
    """

}

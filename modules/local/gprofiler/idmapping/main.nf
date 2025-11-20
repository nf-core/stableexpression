process GPROFILER_IDMAPPING {

    label 'process_medium'

    tag "${species} IDs to ${gprofiler_target_db}"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c28c8e613c062828aaee4b950029bc90a1a1aa94d5f61016a588c8ec7be8b65/data':
        'community.wave.seqera.io/library/pandas_requests_tenacity:5ba56df089a9d718' }"

    input:
    path gene_id_file
    val species
    val gprofiler_target_db

    output:
    path('mapped_gene_ids.csv'),                         optional: true,                                              emit: mapping
    path('gene_metadata.csv'),                           optional: true,                                              emit: metadata
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions

    script:
    """
    gprofiler_map_ids.py \\
        --gene-ids $gene_id_file \\
        --species "$species" \\
        --target-db "$gprofiler_target_db"
    """


    stub:
    """
    touch mapped_gene_ids.csv
    touch gene_metadata.csv
    """

}

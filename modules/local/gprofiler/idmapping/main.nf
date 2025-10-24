process GPROFILER_IDMAPPING {

    label 'process_single'

    tag "${meta.dataset}"

    // limiting to 8 threads at a time to avoid 429 errors with the G Profiler API server
    maxForks 8

    errorStrategy {
        if (task.exitStatus == 100) {
            // ignoring cases when the count dataframe is empty
            log.warn("Count file is empty for dataset ${meta.dataset}.")
            return 'ignore'
        } else if (task.exitStatus == 101) {
            // likewise, when no mapping could be found, we do not want to continue with the subsequent steps for this specific dataset
            log.warn("Could not map gene IDs to Ensembl for dataset ${meta.dataset}.")
            return 'ignore'
        } else if (task.exitStatus == 102) {
            // if the server appears to be down, we stop immediately
            log.error("gProfiler server appears to be down, stopping pipeline")
            return 'terminate'
        } else {
            return 'terminate'
        }
    }

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c28c8e613c062828aaee4b950029bc90a1a1aa94d5f61016a588c8ec7be8b65/data':
        'community.wave.seqera.io/library/pandas_requests_tenacity:5ba56df089a9d718' }"

    input:
    tuple val(meta), path(count_file)
    val species
    val gene_id_mapping_file
    val gene_metadata_file

    output:
    tuple val(meta), path('*.renamed.csv'),                                                                           emit: counts
    path('*.metadata.csv'), optional: true,                                                                           emit: metadata
    path('*.mapping.csv'),  optional: true,                                                                           emit: mapping
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions

    script:
    def custom_mapping_arg  = gene_id_mapping_file ? "--custom-mappings $gene_id_mapping_file" : ""
    def custom_metadata_arg = gene_metadata_file   ? "--custom-metadata $gene_metadata_file" : ""
    """
    map_ids_to_ensembl.py \\
        --count-file "$count_file" \\
        --species "$species" \\
        $custom_mapping_arg \\
        $custom_metadata_arg
    """


    stub:
    """
    touch fake_renamed.csv
    touch fake_metadata.csv
    touch fake_mapping.json
    """

}

process GPROFILER_IDMAPPING {

    label 'process_low'
    debug true

    publishDir "${params.outdir}/idmapping"

    tag "${meta.dataset}"

    // limiting to 8 threads at a time to avoid 429 errors with the G Profiler API server
    maxForks 8

    errorStrategy = {
        if (task.exitStatus == 100) {
            // ignoring cases when the count dataframe is empty
            return 'ignore'
        } else if (task.exitStatus == 101) {
            // likewise, when no mapping could be found, we do not want to continue with the subsequent steps for this specific dataset
            return 'ignore'
        } else if (task.exitStatus == 102) {
            // if the server appears to be down, we stop immediately
            return 'terminate'
        } else {
            return 'terminate'
        }
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas_python_requests:e2bc861ffb0fab18':
        'community.wave.seqera.io/library/pandas_python_requests:8c6da05a2935a952' }"

    input:
    tuple val(meta), path(count_file), val(species)
    val gene_id_mapping_file

    output:
    tuple val(meta), path('*_renamed.csv'),                                                                           emit: renamed
    path('*_metadata.csv'), optional: true,                                                                           emit: metadata
    path('*_mapping.csv'),  optional: true,                                                                           emit: mapping
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions


    script:
    """
    map_ids_to_ensembl.py \
        --count-file "$count_file" \
        --species "$species" \
        --custom-mappings $gene_id_mapping_file
    """


    stub:
    """
    touch fake_renamed.csv
    touch fake_metadata.csv
    touch fake_mapping.json
    """

}

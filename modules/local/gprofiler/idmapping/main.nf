process GPROFILER_IDMAPPING {

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
        } else {
            return 'terminate'
        }
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe3f927f5032b9f0749fd5a2d3431b483f1c8cb1613d0290e2326fec10bf8268/data':
        'community.wave.seqera.io/library/pandas_requests:c7451d98ba573475' }"

    input:
    tuple val(meta), path(count_file), val(species)

    output:
    path('*_renamed.csv'),                                                                                            emit: renamed
    path('*_metadata.csv'),                                                                                           emit: metadata
    path('*_mapping.csv'),                                                                                            emit: mapping
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions


    script:
    """
    map_ids_to_ensembl.py --count-file "$count_file" --species "$species"
    """


    stub:
    """
    touch fake_renamed.csv
    touch fake_metadata.csv
    touch fake_mapping.json
    """

}

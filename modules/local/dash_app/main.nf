process DASH_APP {

    label 'process_single'

    conda "${moduleDir}/app/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path all_counts
    path whole_design
    path genes_stats

    output:
    path("*"), emit: app
    //tuple val("${task.process}"), val('python'),                    eval("python3 --version | sed 's/Python //'"),                                                   topic: versions
    //tuple val("${task.process}"), val('dash'),                      eval('python3 -c "import dash; print(dash.__version__)"'),                                       topic: versions
    //tuple val("${task.process}"), val('dash-ag-grid'),              eval('python3 -c "import dash_ag_grid; print(dash_ag_grid.__version__)"'),                       topic: versions
    //tuple val("${task.process}"), val('dash-extensions'),           eval('python3 -c "import dash_extensions; print(dash_extensions.__version__)"'),                 topic: versions
    //tuple val("${task.process}"), val('dash-mantine-components'),   eval('python3 -c "import dash_mantine_components; print(dash_mantine_components.__version__)"'), topic: versions
    //tuple val("${task.process}"), val('polars'),                    eval('python3 -c "import polars; print(polars.__version__)"'),                                   topic: versions
    //tuple val("${task.process}"), val('pandas'),                    eval('python3 -c "import pandas; print(pandas.__version__)"'),                                   topic: versions
    //tuple val("${task.process}"), val('pyarrow'),                   eval('python3 -c "import pyarrow; print(pyarrow.__version__)"'),                                 topic: versions
    //tuple val("${task.process}"), val('scipy'),                     eval('python3 -c "import scipy; print(scipy.__version__)"'),                                     topic: versions

    script:
    """
    mkdir -p data
    mv ${all_counts} ${whole_design} ${genes_stats} data/
    cp ${moduleDir}/app/* .
    """

}

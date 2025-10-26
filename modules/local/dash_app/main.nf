process DASH_APP {

    label 'process_single'

    conda "${moduleDir}/app/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b3/b39ecd56e298b0ba94bed41bb36d67b0a2bc24634bc53baff9773dcc3d422c01/data':
        'community.wave.seqera.io/library/dash-ag-grid_dash-extensions_dash-iconify_dash-mantine-components_pruned:138d9ff01702db68' }"

    errorStrategy {
        if (task.exitStatus == 100) {
            log.warn("Could not start the Dash application.")
            return 'finish' // finishes started processes but reports error
        }
    }

    input:
    path all_counts
    path whole_design
    path top_stable_genes_summary
    path all_genes_stats

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
    mv ${all_counts} ${whole_design} ${top_stable_genes_summary} ${all_genes_stats} data/
    cp -r ${moduleDir}/app/* .

    # trying to launch the app
    # if the resulting exit code is not 124 (exit code of timeout) then there is an error
    timeout 10 python app.py || exit_code=\$?; [ "\$exit_code" -eq 124 ] && exit 0 || exit 100
    """

}

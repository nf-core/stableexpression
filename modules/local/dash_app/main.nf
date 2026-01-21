process DASH_APP {

    label 'process_high'

    conda "${moduleDir}/app/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4eec747f2063edcc2d1b64e3b84a6b154fde1b9cd9d698446321b4a535432272/data':
        'community.wave.seqera.io/library/dash-ag-grid_dash-extensions_dash-iconify_dash-mantine-components_pruned:7cf6396dd8cd850e' }"

    errorStrategy {
        if (task.exitStatus == 100) {
            log.warn("Could not start the Dash application.")
            return 'ignore' // only report errors but ignores it
        } else {
            log.warn("Could not start the Dash application due to unhandled error.")
            return 'ignore' // ignore anyway
        }
    }

    input:
    path all_counts
    path whole_design
    path all_genes_summary

    output:
    path("*"), emit: app
    path "versions.yml", emit: versions

    script:
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    mkdir -p data
    mv ${all_counts} ${whole_design} ${all_genes_summary} data/
    cp -r ${moduleDir}/app/* .

    # as of Nextflow version 25.04.8, having these versions sent to the versions topic channel
    # results in ERROR ~ No such file or directory: <task workdir>/.command.env
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python3 --version | sed "s/Python //" )
        dash: \$( python3 -c "import dash; print(dash.__version__)" )
        dash-extensions: \$( python3 -c "import dash_extensions; print(dash_extensions.__version__)" )
        dash-mantine-components: \$( python3 -c "import dash_mantine_components; print(dash_mantine_components.__version__)" )
        dash-ag-grid: \$( python3 -c "import dash_ag_grid; print(dash_ag_grid.__version__)" )
        polars: \$( python3 -c "import polars; print(polars.__version__)" )
        pandas: \$( python3 -c "import pandas; print(pandas.__version__)" )
        pyarrow: \$( python3 -c "import pyarrow; print(pyarrow.__version__)" )
        scipy: \$( python3 -c "import scipy; print(scipy.__version__)" )
    END_VERSIONS

    # trying to launch the app
    # if the resulting exit code is not 124 (exit code of timeout) then there is an error
    timeout 10 python -B app.py || exit_code=\$?; [ "\$exit_code" -eq 124 ] && exit 0 || exit 100
    """

}

process IDMAPPING {

    debug true

    // limiting to 1 thread at a time to avoid crashing the G Profiler API server
    maxForks 1

    conda "${moduleDir}/environment.yml"

    input:
    tuple path(count_file), val(species)

    output:
    path '*_renamed.csv',                                                                                             emit: csv
    tuple val("${task.process}"), val('python'),   eval('python3 --version'),                                         topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    map_ids_to_ensembl.py --count-file "$count_file" --species "$species"
    """


    stub:
    """
    touch fake_renamed.csv
    """

}

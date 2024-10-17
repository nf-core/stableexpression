process MERGE_DESIGNS {

    debug true

    input:
    path csv_files

    output:
    path 'all_designs.csv',                                                                                           emit: csv
    tuple val("${task.process}"), val('python'),   eval('python3 --version'),                                         topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    concatenate_csv_files.py $csv_files --outfile all_designs.csv
    """

}

process EXPRESSIONATLAS_GETACCESSIONS {

    debug true

    conda "${moduleDir}/environment.yml"


    input:
    tuple val(species), val(keywords)

    output:
    path 'accessions.csv',                                                                                            emit: csv
    tuple val("${task.process}"), val('python'),   eval('python3 --version'),                                         topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions
    tuple val("${task.process}"), val('nltk'),     eval('python3 -c "import nltk; print(nltk.__version__)"'),         topic: versions



    when:
    task.ext.when == null || task.ext.when

    script:
    if (keywords == null) {
        """
        get_expression_atlas_accessions.py --species $species
        """
    } else {
        """
        get_expression_atlas_accessions.py --species "$species" --keywords $keywords
        """
    }


    stub:
    """
    touch accessions.csv
    touch found.json
    """

}

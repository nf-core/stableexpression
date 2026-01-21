process GEO_GETACCESSIONS {

    label 'process_high'

    tag "${species}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/caae35ec5dc72367102a616a47b6f1a7b3de9ff272422f2c08895b8bb5f0566c/data':
        'community.wave.seqera.io/library/biopython_nltk_pandas_parallelbar_pruned:5fc501b07f8e0428' }"

    input:
    val species
    val keywords
    val platform
    path excluded_accessions_file
    val random_sampling_size
    val random_sampling_seed

    output:
    path "accessions.txt",                     optional: true,          emit: accessions
    path "geo_selected_datasets.metadata.tsv", optional: true,          topic: geo_selected_datasets
    path "geo_all_datasets.metadata.tsv",      optional: true,          topic: geo_all_datasets
    path "geo_rejected_datasets.metadata.tsv", optional: true,          topic: geo_rejected_datasets

    tuple val("${task.process}"), val('python'),      eval("python3 --version | sed 's/Python //'"),                          topic: versions
    tuple val("${task.process}"), val('requests'),    eval('python3 -c "import requests; print(requests.__version__)"'),      topic: versions
    tuple val("${task.process}"), val('nltk'),        eval('python3 -c "import nltk; print(nltk.__version__)"'),              topic: versions
    tuple val("${task.process}"), val('pandas'),      eval('python3 -c "import pandas; print(pandas.__version__)"'),          topic: versions
    tuple val("${task.process}"), val('biopython'),   eval('python3 -c "import Bio; print(Bio.__version__)"'),                topic: versions
    tuple val("${task.process}"), val('tqdm'),        eval('python3 -c "import tqdm; print(tqdm.__version__)"'),              topic: versions

    script:
    def keywords_string = keywords.split(',').collect { it.trim() }.join(' ')
    def args = " --species $species"
    if ( keywords_string != "" ) {
        args += " --keywords $keywords_string"
    }
    if ( platform ) {
        args += " --platform $platform"
    }
    if ( excluded_accessions_file ) {
        args += " --exclude-accessions-in $excluded_accessions_file"
    }
    if ( random_sampling_size ) {
        args += " --random-sampling-size $random_sampling_size"
    }
    if ( random_sampling_seed ) {
        args += " --random-sampling-seed $random_sampling_seed"
    }
    // the folder where nltk will download data needs to be writable (necessary for singularity)
    """
    # the Entrez module from biopython automatically stores temp results in <home dir>/.config
    # if this directory is not writable, the script fails
    export HOME=/tmp/biopython
    mkdir -p /tmp/biopython

    export NLTK_DATA=\${PWD}

    get_geo_dataset_accessions.py \\
        $args \\
        --cpus ${task.cpus}
    """

    stub:
    """
    touch accessions.txt \\
        all_experiments.metadata.tsv \\
        filtered_experiments.metadata.tsv \\
        filtered_experiments.keywords.yaml
    """

}

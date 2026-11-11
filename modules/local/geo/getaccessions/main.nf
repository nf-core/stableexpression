process GEO_GETACCESSIONS {

    label 'process_high'

    tag "${species}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e8be45bdbe57d56f7d452513c4799a878fdfeb2f8ff8351f1c02ee99627dc50e/data':
        'community.wave.seqera.io/library/biopython_httpx_nltk_pandas_pruned:f692df8e1f55b14b' }"

    input:
    val species
    path excluded_accessions_file

    output:
    path "accessions.txt",                     optional: true,                                                     emit: accessions
    path "geo_selected_datasets.metadata.tsv", optional: true,                                                     topic: geo_selected_datasets
    path "geo_all_datasets.metadata.tsv",      optional: true,                                                     topic: geo_all_datasets
    path "geo_rejected_datasets.metadata.tsv", optional: true,                                                     topic: geo_rejected_datasets
    tuple val("${task.process}"), val('python'),    eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('httpx'),     eval('python3 -c "import httpx; print(httpx.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('nltk'),      eval('python3 -c "import nltk; print(nltk.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('pandas'),    eval('python3 -c "import pandas; print(pandas.__version__)"'), topic: versions
    tuple val("${task.process}"), val('biopython'), eval('python3 -c "import Bio; print(Bio.__version__)"'),       topic: versions
    tuple val("${task.process}"), val('tqdm'),      eval('python3 -c "import tqdm; print(tqdm.__version__)"'),     topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    # limiting CPU usage
    export OMP_NUM_THREADS=${task.cpus}

    # the Entrez module from biopython automatically stores temp results in <home dir>/.config
    # if this directory is not writable, the script fails
    export HOME=/tmp/biopython
    mkdir -p /tmp/biopython

    # the folder where nltk will download data needs to be writable (necessary for singularity)
    export NLTK_DATA=\${PWD}

    get_geo_dataset_accessions.py \\
        --species $species \\
        --cpus ${task.cpus} \\
        $args
    """

    stub:
    """
    touch accessions.txt
    touch geo_selected_datasets.metadata.tsv
    touch geo_all_datasets.metadata.tsv
    touch geo_rejected_datasets.metadata.tsv
    """
}

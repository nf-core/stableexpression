process GEO_GETACCESSIONS {

    label 'process_medium_requirements'

    tag "${species}"

    conda "${moduleDir}/environment.yml"
    // NOTE: images were built using the wave cli, which does not allow getting the https URL for singularity containers
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_httpx_tenacity_nltk_pruned:5abbf76aa66ff963':
        'community.wave.seqera.io/library/python_httpx_tenacity_nltk_pruned:8fae761ed7b56b59' }"

    input:
    val species
    path excluded_accessions_file

    output:
    path "accessions.txt",                     optional: true,                                                     emit: accessions
    path "geo_selected_datasets.metadata.tsv", optional: true,                                                     topic: geo_selected_datasets
    path "geo_all_datasets.metadata.tsv",      optional: true,                                                     topic: geo_all_datasets
    path "geo_rejected_datasets.metadata.tsv", optional: true,                                                     topic: geo_rejected_datasets
    tuple val("${task.process}"), val('python'),                eval("python3 --version | sed 's/Python //'"),                                          topic: versions
    tuple val("${task.process}"), val('httpx'),                 eval('python3 -c "import httpx; print(httpx.__version__)"'),                            topic: versions
    tuple val("${task.process}"), val('nltk'),                  eval('python3 -c "import nltk; print(nltk.__version__)"'),                              topic: versions
    tuple val("${task.process}"), val('biopython'),             eval('python3 -c "import Bio; print(Bio.__version__)"'),                                topic: versions
    tuple val("${task.process}"), val('tqdm'),                  eval('python3 -c "import tqdm; print(tqdm.__version__)"'),                              topic: versions
    tuple val("${task.process}"), val('sentence_transformers'), eval('python3 -c "import sentence_transformers as st_tfm; print(st_tfm.__version__)"'), topic: versions
    tuple val("${task.process}"), val('scispacy'),              eval('python3 -c "import scispacy; print(scispacy.__version__)"'),                      topic: versions

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

    # the folders where to store files & models
    export MESH_FOLDER=\${PWD}/mesh
    export SENTENCE_TRANSFORMERS_HOME=\${PWD}/sentence_transformers_cache

    # necessary to avoid errors with huggingface
    export HF_HUB_DISABLE_XET=1

    get_geo_dataset_accessions.py \\
        --species $species \\
        --cpus ${task.cpus} \\
        $args

    # removing nltk folders
    rm -rf \${PWD}/corpora \${PWD}/taggers \${PWD}/tokenizers
    # removing Mesh data
    rm -rf \$MESH_FOLDER
    # removing sentence transformers cache
    rm -rf \$SENTENCE_TRANSFORMERS_HOME
    """

    stub:
    """
    touch accessions.txt
    touch geo_selected_datasets.metadata.tsv
    touch geo_all_datasets.metadata.tsv
    touch geo_rejected_datasets.metadata.tsv
    """
}

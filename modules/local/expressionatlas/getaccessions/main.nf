process EXPRESSIONATLAS_GETACCESSIONS {
    debug true
    label 'process_low_requirements'

    tag "${species}"

    conda "${moduleDir}/environment.yml"
    // NOTE: images were built using the wave cli, which does not allow getting the https URL for singularity containers
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_httpx_tenacity_nltk_pruned:1fa403d74b162c3b':
        'community.wave.seqera.io/library/python_httpx_tenacity_nltk_pruned:4a1db8e45d12c0ea' }"

    input:
    val species

    output:
    path "accessions.txt",                    optional: true, emit: accessions
    env("SAMPLING_QUOTA"),                                    emit: sampling_quota
    path "selected_experiments.metadata.tsv", optional: true, topic: eatlas_selected_datasets
    path "rejected_experiments.metadata.tsv", optional: true, topic: eatlas_rejected_datasets
    tuple val("${task.process}"), val('python'),                eval("python3 --version | sed 's/Python //'"),                                          topic: versions
    tuple val("${task.process}"), val('httpx'),                 eval('python3 -c "import httpx; print(httpx.__version__)"'),                            topic: versions
    tuple val("${task.process}"), val('nltk'),                  eval('python3 -c "import nltk; print(nltk.__version__)"'),                              topic: versions
    tuple val("${task.process}"), val('sentence_transformers'), eval('python3 -c "import sentence_transformers as st_tfm; print(st_tfm.__version__)"'), topic: versions
    tuple val("${task.process}"), val('scispacy'),              eval('python3 -c "import scispacy; print(scispacy.__version__)"'),                      topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    # limiting CPU usage
    export OMP_NUM_THREADS=${task.cpus}

    # the folder where nltk will download data needs to be writable (necessary for singularity)
    export NLTK_DATA=\${PWD}

    # the folders where to store files & models
    export MESH_FOLDER=\${PWD}/mesh
    export SENTENCE_TRANSFORMERS_HOME=\${PWD}/sentence_transformers_cache

    # necessary to avoid errors with huggingface
    export HF_HUB_DISABLE_XET=1

    get_eatlas_accessions.py \\
        --species $species \\
        $args

    SAMPLING_QUOTA=\$(cat sampling_quota.txt)
    """

    stub:
    """
    echo "E-stub" > accessions.txt
    SAMPLING_QUOTA=ok
    touch selected_experiments.metadata.tsv
    touch rejected_experiments.metadata.tsv
    """
}

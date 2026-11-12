process EXPRESSIONATLAS_GETACCESSIONS {

    label 'process_high'

    tag "${species}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f943893a85d82f720432e83fd8d4755e5b42a92deca9d49c06930eaa7fc0c968/data':
        'community.wave.seqera.io/library/httpx_nltk_pandas_python_pruned:ab2f10d1d67a7603' }"

    input:
    val species

    output:
    path "accessions.txt",                    optional: true,                                                   emit: accessions
    env("SAMPLING_QUOTA"),                                                                                      emit: sampling_quota
    path "selected_experiments.metadata.tsv", optional: true,                                                   topic: eatlas_selected_datasets
    path "species_experiments.metadata.tsv",  optional: true,                                                   topic: eatlas_all_datasets
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('httpx'),  eval('python3 -c "import httpx; print(httpx.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('nltk'),   eval('python3 -c "import nltk; print(nltk.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('pyyaml'), eval('python3 -c "import yaml; print(yaml.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'), topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    # limiting CPU usage
    export OMP_NUM_THREADS=${task.cpus}

    # the folder where nltk will download data needs to be writable (necessary for singularity)
    export NLTK_DATA=\${PWD}

    get_eatlas_accessions.py \\
        --species $species \\
        --cpus ${task.cpus} \\
        $args

    SAMPLING_QUOTA=\$(cat sampling_quota.txt)
    """

    stub:
    """
    echo "E-stub" > accessions.txt
    SAMPLING_QUOTA=ok
    touch selected_experiments.metadata.tsv
    touch species_experiments.metadata.tsv
    """
}

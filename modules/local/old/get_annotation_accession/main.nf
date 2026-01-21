process GET_ANNOTATION_ACCESSION {

    label 'process_single'

    tag "$species"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b4d686ef63e22bc4d461178fc241cefddd2aa3436e189d3787c8e019448f056e/data':
        'community.wave.seqera.io/library/requests_tenacity_tqdm:126dbed8ef3ff96f' }"

    input:
    val(species)

    output:
    env("ACCESSION"), emit: accession
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                                               topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'),                           topic: versions
    tuple val("${task.process}"), val('tenacity'), eval('python3 -c "from importlib.metadata import version; print(version(\'tenacity\'))"'),   topic: versions

    script:
    """
    get_annotation_accession.py --species $species
    ACCESSION=\$(cat accession.txt)
    """

    stub:
    """
    touch accession.txt
    """

}

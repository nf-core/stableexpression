process COMPUTE_GENE_TRANSCRIPT_LENGTHS {

    label 'process_single'

    tag "${gff3.baseName}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3d7126100b0eb7cb53dfb50291707ea8dda3b9738b76551ab73605d0acbe114b/data':
        'community.wave.seqera.io/library/pandas:2.3.3--5a902bf824a79745' }"

    input:
    path gff3

    output:
    path('gene_transcript_lengths.csv'),                                                                              emit: csv
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions

    script:
    def is_compressed = gff3.getExtension() == "gz" ? true : false
    def gff3_name = is_compressed ? gff3.getBaseName() : gff3
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${gff3} > ${gff3_name}
    fi

    compute_gene_transcript_lengths.py \\
        --annotation ${gff3_name}
    """


    stub:
    """
    touch gene_transcript_lengths.csv
    """

}

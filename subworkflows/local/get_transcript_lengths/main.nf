include { COMPUTE_GENE_TRANSCRIPT_LENGTHS              } from '../../../modules/local/compute_gene_transcript_lengths'
include { DOWNLOAD_ENSEMBL_ANNOTATION                  } from '../../../modules/local/download_ensembl_annotation'


/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow GET_TRANSCRIPT_LENGTHS {

    take:
    species

    main:

    DOWNLOAD_ENSEMBL_ANNOTATION (species)
    ch_annotation = DOWNLOAD_ENSEMBL_ANNOTATION.out.gff3

    COMPUTE_GENE_TRANSCRIPT_LENGTHS (ch_annotation)



    emit:
    csv = COMPUTE_GENE_TRANSCRIPT_LENGTHS.out.csv



}

include { COMPUTE_GENE_TRANSCRIPT_LENGTHS              } from '../../../modules/local/compute_gene_transcript_lengths'
include { DOWNLOAD_ENSEMBL_ANNOTATION                  } from '../../../modules/local/download_ensembl_annotation'


/*
========================================================================================
    SUBWORKFLOW TO GET TRANSCRIPT LENGTHS FROM GENOME ANNOTATION
========================================================================================
*/

workflow GET_TRANSCRIPT_LENGTHS {

    take:
    species
    gff_file
    gff_url

    main:

    if ( gff_file ) {
        ch_annotation = channel.fromPath( gff_file, checkIfExists: true )
    } else if ( gff_url ) {
        ch_annotation = channel.fromPath( gff_url, checkIfExists: true )
    } else {
        DOWNLOAD_ENSEMBL_ANNOTATION( species )
        ch_annotation = DOWNLOAD_ENSEMBL_ANNOTATION.out.gff3
    }

    COMPUTE_GENE_TRANSCRIPT_LENGTHS( ch_annotation )



    emit:
    csv = COMPUTE_GENE_TRANSCRIPT_LENGTHS.out.csv



}

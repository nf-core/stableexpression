include { IDMAPPING_GPROFILER                    } from '../../../modules/local/idmapping/gprofiler'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow IDMAPPING {

    take:
    ch_datasets
    species

    main:

    def ch_gene_metadata = params.gene_metadata ? Channel.fromPath( params.gene_metadata, checkIfExists: true ) : Channel.empty()

    if ( params.skip_gprofiler ) {

        def ch_gene_id_mapping = params.gene_id_mapping_file ? Channel.fromPath( params.gene_id_mapping, checkIfExists: true ) : Channel.empty()

    } else {

        // tries to map gene IDs to Ensembl IDs whenever possible
        IDMAPPING_GPROFILER(
            ch_datasets,
            species,
            params.gene_id_mapping_file ? Channel.fromPath( params.gene_id_mapping_file, checkIfExists: true ) : 'none'
        )

        IDMAPPING_GPROFILER.out.renamed.set { ch_datasets }

        ch_gene_metadata
            .mix( IDMAPPING_GPROFILER.out.metadata )
            .set { ch_gene_metadata }

        // the gene id mappings are the sum
        // of those provided by the user and those fetched from g:Profiler
        IDMAPPING_GPROFILER.out.mapping.set { ch_gene_id_mapping }

    }

    emit:
    datasets                        = ch_datasets
    gene_metadata                   = ch_gene_metadata
    gene_id_mapping                 = ch_gene_id_mapping
}

include { IDMAPPING_GPROFILER                    } from '../../../modules/local/idmapping/gprofiler'

/*
========================================================================================
    SUBWORKFLOW TO MAP ORIGINAL IDS TO ENSEMBL GENE IDS
========================================================================================
*/

workflow IDMAPPING {

    take:
    ch_datasets
    ch_species

    main:

    ch_gene_metadata = params.gene_metadata ? Channel.fromPath( params.gene_metadata, checkIfExists: true ) : Channel.empty()
    ch_gene_id_mapping = params.gene_id_mapping_file ? Channel.fromPath( params.gene_id_mapping, checkIfExists: true ) : Channel.empty()

    if ( !params.skip_gprofiler ) {

        // tries to map gene IDs to Ensembl IDs whenever possible
        IDMAPPING_GPROFILER(
            ch_datasets,
            ch_species,
            params.gene_id_mapping_file ? Channel.fromPath( params.gene_id_mapping_file, checkIfExists: true ) : Channel.value( [] )
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

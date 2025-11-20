include { COLLECT_GENE_IDS                       } from '../../../modules/local/collect_gene_ids'
include { GPROFILER_IDMAPPING                    } from '../../../modules/local/gprofiler/idmapping'
include { RENAME_GENE_IDS                        } from '../../../modules/local/rename_gene_ids'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow ID_MAPPING {

    take:
    ch_counts
    species
    skip_id_mapping
    gprofiler_target_db
    ch_custom_gene_id_mapping
    ch_custom_gene_metadata


    main:
    ch_counts.view { a -> "counts ${a}"}
    ch_gene_id_mapping = Channel.empty()

    if ( !params.skip_id_mapping ) {

        COLLECT_GENE_IDS(
            ch_counts.map{ meta, file -> file }.collect()
        )

        GPROFILER_IDMAPPING(
            COLLECT_GENE_IDS.out.gene_ids,
            species,
            gprofiler_target_db
        )
        GPROFILER_IDMAPPING.out.mapping.set { ch_gene_id_mapping }
    }

    RENAME_GENE_IDS(
        ch_counts,
        ch_gene_id_mapping,
        ch_custom_gene_id_mapping
    )

    emit:
    counts          = RENAME_GENE_IDS.out.counts
    mapping         = ch_gene_id_mapping
    metadata        = GPROFILER_IDMAPPING.out.metadata

}

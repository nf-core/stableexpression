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
    outdir


    main:

    ch_gene_id_mapping = Channel.empty()
    ch_gene_metadata = Channel.empty()

    if ( !params.skip_id_mapping ) {

        // -----------------------------------------------------------------
        // COLLECTING ALL GENE IDS FROm ALL DATASETS
        // -----------------------------------------------------------------

        COLLECT_GENE_IDS(
            ch_counts.map{ meta, file -> file }.collect()
        )

        // -----------------------------------------------------------------
        // MAPPING THESE GENE IDS TO THE CHOSEN TARGET DB
        // -----------------------------------------------------------------

        GPROFILER_IDMAPPING(
            COLLECT_GENE_IDS.out.gene_ids,
            species,
            gprofiler_target_db
        )
        GPROFILER_IDMAPPING.out.mapping.set { ch_gene_id_mapping }
        GPROFILER_IDMAPPING.out.metadata.set { ch_gene_metadata }
    }

    // -----------------------------------------------------------------
    // RENAMING GENE IDS IN ALL COUNT DATASETS
    // -----------------------------------------------------------------

    RENAME_GENE_IDS(
        ch_counts,
        ch_gene_id_mapping,
        ch_custom_gene_id_mapping
    )

    // -----------------------------------------------------------------
    // COLLECTING GLOBAL GENE ID MAPPING AND METADATA
    // -----------------------------------------------------------------

    ch_gene_id_mapping
        .mix( ch_custom_gene_id_mapping )
        .filter { it != [] } // handle no custom mappings
        .splitCsv( header: true )
        .unique()
        .collectFile(
            name: 'global_gene_id_mapping.csv',
            seed: "original_gene_id,gene_id",
            newLine: true,
            storeDir: "${params.outdir}/idmapping/"
        ) {
            item -> "${item["original_gene_id"]},${item["gene_id"]}"
        }
        .ifEmpty([])
        .set { ch_global_gene_id_mapping }

    ch_gene_metadata
        .mix( ch_custom_gene_metadata )
        .filter { it != [] } // handle no custom metadata
        .splitCsv( header: true )
        .unique()
        .collectFile(
            name: 'global_gene_metadata.csv',
            seed: "gene_id,name,description",
            newLine: true,
            storeDir: "${params.outdir}/idmapping/"
        ) {
            item -> "${item["gene_id"]},${item["name"]},${item["description"]}"
        }
        .ifEmpty([])
        .set { ch_global_gene_metadata }

    emit:
    counts          = RENAME_GENE_IDS.out.counts
    mapping         = ch_global_gene_id_mapping
    metadata        = ch_global_gene_metadata

}

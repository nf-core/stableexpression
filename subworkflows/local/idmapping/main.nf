include { CLEAN_GENE_IDS                         } from '../../../modules/local/clean_gene_ids'
include { COLLECT_GENE_IDS                       } from '../../../modules/local/collect_gene_ids'
include { GPROFILER_IDMAPPING                    } from '../../../modules/local/gprofiler/idmapping'
include { FILTER_OUT_RARE_GENES                  } from '../../../modules/local/filter_out_rare_genes'
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
    custom_gene_id_mapping
    custom_gene_metadata
    min_freq_occurrence
    outdir

    main:

    ch_gene_id_mapping      = channel.empty()
    ch_gene_metadata        = channel.empty()
    ch_valid_gene_ids       = channel.empty()

    if ( !skip_id_mapping ) {

        // -----------------------------------------------------------------
        // CLEANING GENE IDS
        // -----------------------------------------------------------------

        CLEAN_GENE_IDS ( ch_counts )
        ch_counts           = CLEAN_GENE_IDS.out.counts
        ch_cleaned_gene_ids = CLEAN_GENE_IDS.out.gene_ids

        // -----------------------------------------------------------------
        // COLLECTING ALL CLEANED GENE IDS FROM ALL DATASETS
        // -----------------------------------------------------------------

        // sorting files in order to have a consistent input and be able to retry
        COLLECT_GENE_IDS(
            ch_cleaned_gene_ids.toSortedList()
        )

        // -----------------------------------------------------------------
        // MAPPING THESE GENE IDS TO THE CHOSEN TARGET DB
        // -----------------------------------------------------------------

        GPROFILER_IDMAPPING(
            COLLECT_GENE_IDS.out.unique_gene_ids,
            species,
            gprofiler_target_db
        )
        ch_gene_id_mapping      = GPROFILER_IDMAPPING.out.mapping
        ch_gene_metadata        = GPROFILER_IDMAPPING.out.metadata

        // -----------------------------------------------------------------
        // FILTERING OUT GENE IDS THAT DO NOT HAVE ENOUGH OCCURRENCES
        // -----------------------------------------------------------------

        FILTER_OUT_RARE_GENES(
            ch_gene_id_mapping,
            COLLECT_GENE_IDS.out.gene_id_occurrences,
            ch_counts.count(),
            min_freq_occurrence
        )
        ch_valid_gene_ids = FILTER_OUT_RARE_GENES.out.valid_gene_ids
    }

    // -----------------------------------------------------------------
    // COLLECTING GLOBAL GENE ID MAPPING AND METADATA
    // -----------------------------------------------------------------

    ch_global_gene_id_mapping = ch_gene_id_mapping
                                    .mix(
                                        custom_gene_id_mapping ?
                                        channel.fromPath( custom_gene_id_mapping, checkIfExists: true ) :
                                        channel.empty()
                                    )
                                    .splitCsv( header: true )
                                    .unique()
                                    .collectFile(
                                        name: 'global_gene_id_mapping.csv',
                                        seed: "original_gene_id,gene_id",
                                        newLine: true,
                                        storeDir: "${outdir}/idmapping/",
                                        sort: true
                                    ) {
                                        item -> "${item["original_gene_id"]},${item["gene_id"]}"
                                    }

    ch_global_gene_metadata = ch_gene_metadata
                                .mix(
                                    custom_gene_metadata ?
                                    channel.fromPath( custom_gene_metadata, checkIfExists: true ) :
                                    channel.empty()
                                )
                                .splitCsv( header: true )
                                .unique()
                                .collectFile(
                                    name: 'global_gene_metadata.csv',
                                    seed: "gene_id,name,description",
                                    newLine: true,
                                    storeDir: "${outdir}/idmapping/",
                                    sort: true
                                ) {
                                    item -> "${item["gene_id"]},${item["name"]},${item["description"]}"
                                }

    // -----------------------------------------------------------------
    // RENAMING GENE IDS IN ALL COUNT DATASETS (ONLY IF NECESSARY)
    // -----------------------------------------------------------------

    if ( !skip_id_mapping || custom_gene_id_mapping ) {

        RENAME_GENE_IDS(
            ch_counts,
            ch_global_gene_id_mapping.first(),
            ch_valid_gene_ids.collect()
        )
        ch_counts = RENAME_GENE_IDS.out.counts

    }


    emit:
    counts          = ch_counts
    mapping         = ch_global_gene_id_mapping
    metadata        = ch_global_gene_metadata

}

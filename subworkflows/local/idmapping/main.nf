include { CLEAN_GENE_IDS                         } from '../../../modules/local/clean_gene_ids'
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
    custom_gene_id_mapping
    custom_gene_metadata
    outdir

    main:

    ch_gene_id_mapping = Channel.empty()
    ch_gene_metadata = Channel.empty()

    if ( !skip_id_mapping ) {

        // -----------------------------------------------------------------
        // COLLECTING ALL GENE IDS FROM ALL DATASETS
        // -----------------------------------------------------------------

        // here we cannot use directly COLLECT_GENE_IDS for runs comprising a huge number of files (eg. human)
        // so that we proceed by chunks, and perform a final merging step using the Java VM

        CLEAN_GENE_IDS ( ch_counts )
        ch_counts = CLEAN_GENE_IDS.out.counts

        // TRICK:
        // the buffer operator creates non-deterministic chunks
        // which prevents resuming the pipeline
        // so we sort the list of files before buffering them
        ch_chunck_counts = ch_counts
                            .map{ meta, file -> file }
                            .collect( sort: true ) // get all files and sort them
                            .flatten() // needed to convert the list back to individual channel items
                            .buffer( size: 100, remainder: true )

        COLLECT_GENE_IDS( ch_chunck_counts )

        ch_gene_ids = COLLECT_GENE_IDS.out.gene_ids
                        .splitText()
                        .unique()
                        .collectFile(
                            name: 'original_gene_ids.txt',
                            storeDir: "${outdir}/idmapping/"
                        )

        // -----------------------------------------------------------------
        // MAPPING THESE GENE IDS TO THE CHOSEN TARGET DB
        // -----------------------------------------------------------------

        GPROFILER_IDMAPPING(
            ch_gene_ids,
            species,
            gprofiler_target_db
        )
        GPROFILER_IDMAPPING.out.mapping.set { ch_gene_id_mapping }
        GPROFILER_IDMAPPING.out.metadata.set { ch_gene_metadata }
    }

    // -----------------------------------------------------------------
    // COLLECTING GLOBAL GENE ID MAPPING AND METADATA
    // -----------------------------------------------------------------

    ch_gene_id_mapping
        .mix( custom_gene_id_mapping ? Channel.fromPath( custom_gene_id_mapping, checkIfExists: true ) : Channel.empty() )
        .splitCsv( header: true )
        .unique()
        .collectFile(
            name: 'global_gene_id_mapping.csv',
            seed: "original_gene_id,gene_id",
            newLine: true,
            storeDir: "${outdir}/idmapping/"
        ) {
            item -> "${item["original_gene_id"]},${item["gene_id"]}"
        }
        .set { ch_global_gene_id_mapping }

    ch_gene_metadata
        .mix( custom_gene_metadata ? Channel.fromPath( custom_gene_metadata, checkIfExists: true ) : Channel.empty() )
        .splitCsv( header: true )
        .unique()
        .collectFile(
            name: 'global_gene_metadata.csv',
            seed: "gene_id,name,description",
            newLine: true,
            storeDir: "${outdir}/idmapping/"
        ) {
            item -> "${item["gene_id"]},${item["name"]},${item["description"]}"
        }
        .set { ch_global_gene_metadata }

    // -----------------------------------------------------------------
    // RENAMING GENE IDS IN ALL COUNT DATASETS (ONLY IF NECESSARY)
    // -----------------------------------------------------------------

    if ( !skip_id_mapping || custom_gene_id_mapping ) {

        RENAME_GENE_IDS(
            ch_counts,
            ch_global_gene_id_mapping.first()
        )
        ch_counts = RENAME_GENE_IDS.out.counts

    }


    emit:
    counts          = ch_counts
    mapping         = ch_global_gene_id_mapping
    metadata        = ch_global_gene_metadata

}

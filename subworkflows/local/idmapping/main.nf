include { GPROFILER_IDMAPPING                    } from '../../../modules/local/gprofiler/idmapping'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow ID_MAPPING {

    take:
    ch_counts
    species
    ch_gene_id_mapping
    ch_gene_metadata


    main:

    ch_counts
        .map {
            meta, file ->
                def platform_taxon = meta.platform_taxon ?: species
                meta.platform_taxon = platform_taxon
                [ meta, file ]
        }
        .set { ch_counts }

    GPROFILER_IDMAPPING(
        ch_counts,
        species,
        ch_gene_id_mapping,
        ch_gene_metadata
    )

    GPROFILER_IDMAPPING.out.counts.set { ch_counts }
    GPROFILER_IDMAPPING.out.mapping.set { ch_gene_id_mapping }
    GPROFILER_IDMAPPING.out.metadata.set { ch_gene_metadata }

    emit:
    counts          = GPROFILER_IDMAPPING.out.counts
    mapping         = GPROFILER_IDMAPPING.out.mapping
    metadata        = GPROFILER_IDMAPPING.out.metadata

}

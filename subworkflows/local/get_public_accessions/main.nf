include { EXPRESSIONATLAS_GETACCESSIONS as EXPRESSION_ATLAS          } from '../../../modules/local/expressionatlas/getaccessions'
include { GEO_GETACCESSIONS as GEO                                   } from '../../../modules/local/geo/getaccessions'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow GET_PUBLIC_ACCESSIONS {

    take:
    species
    skip_fetch_eatlas_accessions
    fetch_geo_accessions
    platform
    keywords
    ch_accessions
    ch_accessions_file
    ch_excluded_accessions
    ch_excluded_accessions_file
    random_sampling_size
    random_sampling_seed
    outdir

    main:

    ch_fetched_eatlas_accessions = channel.empty()
    ch_fetched_geo_accessions    = channel.empty()
    ch_sampling_quota            = channel.of( "ok" )

    // -----------------------------------------------------------------
    // GET EATLAS ACCESSIONS
    // -----------------------------------------------------------------

    // fetching Expression Atlas accessions if applicable
    if ( !skip_fetch_eatlas_accessions ) {

        // getting Expression Atlas accessions given a species name and keywords
        // keywords can be an empty string
        EXPRESSION_ATLAS(
            species,
            keywords,
            platform?: [],
            random_sampling_size?: [],
            random_sampling_seed?: []
        )

        ch_fetched_eatlas_accessions = EXPRESSION_ATLAS.out.accessions.splitText()
        ch_sampling_quota            = EXPRESSION_ATLAS.out.sampling_quota

    }

    // ------------------------------------------------------------------------------------
    // GET GEO ACCESSIONS
    // ------------------------------------------------------------------------------------

    // fetching GEO accessions if applicable
    if ( fetch_geo_accessions ) {

        // all Expression Atlas accessions starting with E-GEOD- are imported from GEO
        // we do not want to collect these GEO data if we already get them from Expression Atlas
        ch_excluded_eatlas_accessions_file = ch_fetched_eatlas_accessions
            .filter { accession -> accession.startsWith("E-GEOD-") }
            .map { accession -> accession.replace("E-GEOD-", "GSE") }
            .collectFile(
                name: 'excluded_geo_accessions.txt',
                storeDir: "${outdir}/geo/",
                sort: true,
                newLine: true
            )
            .ifEmpty( [] )

        // trick to avoid fetching accessions from GEO when the sampling quota is already exceeded
        ch_species = channel.of( species )
                        .combine( ch_sampling_quota )
                        .filter { species_name, quota -> quota == "ok" }
                        .map { species_name, quota -> species_name }

        // getting GEO accessions given a species name and keywords
        // keywords can be an empty string
        GEO(
            ch_species,
            keywords,
            platform?: [],
            ch_excluded_eatlas_accessions_file,
            random_sampling_size?: [],
            random_sampling_seed?: []
        )

        ch_fetched_geo_accessions = GEO.out.accessions.splitText()
    }

    // -----------------------------------------------------------------
    // MERGING AND EXCLUDING UNWANTED ACCESSIONS
    // -----------------------------------------------------------------

    // getting accessions to exclude and preparing in the right format
    ch_excluded_accessions = ch_excluded_accessions
                                .mix( ch_excluded_accessions_file.splitText() )
                                .unique()
                                .map { acc -> acc.trim() }
                                .toList()
                                .map { lst -> [lst] } // list of lists : mandatory when combining in the next step

    ch_fetched_public_accessions = ch_fetched_eatlas_accessions
                                    .mix( ch_fetched_geo_accessions )
                                    .map { acc -> acc.trim() }
                                    .filter { acc ->
                                        (acc.startsWith('E-') || acc.startsWith('GSE')) && !acc.startsWith('E-PROT-')
                                    }
                                    .combine ( ch_excluded_accessions )
                                    .filter { accession, excluded_accessions -> !(accession in excluded_accessions) }
                                    .map { accession, excluded_accessions -> accession }

    // -----------------------------------------------------------------
    // ADDING USER PROVIDED ACCESSIONS
    // -----------------------------------------------------------------

    ch_input_accessions = ch_accessions
                            .mix( ch_accessions_file.splitText() )
                            .unique()
                            .map { acc -> acc.trim() }

    // appending to accessions provided by the user
    // ensures that no accessions is present twice (provided by the user and fetched from E. Atlas)
    // removing E-PROT- accessions because they are not supported in subsequent steps
    // removing excluded accessions
    ch_all_accessions = ch_input_accessions
                        .mix( ch_fetched_public_accessions )
                        .unique()
                        .map { acc -> acc.trim() }

    emit:
    accessions          = ch_all_accessions

}

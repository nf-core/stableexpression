include { EXPRESSIONATLAS_GETACCESSIONS          } from '../../../modules/local/expressionatlas/getaccessions'
include { EXPRESSIONATLAS_GETDATA                } from '../../../modules/local/expressionatlas/getdata'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow EXPRESSIONATLAS_FETCHDATA {

    take:
    ch_species


    main:

    ch_eatlas_datasets = Channel.empty()
    ch_fetched_accessions = Channel.empty()

    ch_eatlas_accessions_file = params.eatlas_accessions_file ? Channel.fromPath(params.eatlas_accessions_file, checkIfExists: true) : Channel.empty()

    Channel.fromList( params.eatlas_accessions.tokenize(',') )
        .mix( ch_eatlas_accessions_file.splitText() )
        .unique()
        .map { it -> it.trim() }
        .set { ch_input_accessions }

    // fetching Expression Atlas accessions if applicable
    if ( !params.skip_fetch_eatlas_accessions || params.eatlas_keywords ) {

        // getting Expression Atlas accessions given a species name and keywords
        // keywords can be an empty string
        EXPRESSIONATLAS_GETACCESSIONS(
            ch_species,
            params.eatlas_keywords
        )
        EXPRESSIONATLAS_GETACCESSIONS.out.txt.set { ch_fetched_accessions }

    }

    ch_exclude_eatlas_accessions_file = params.exclude_eatlas_accessions_file ? Channel.fromPath(params.exclude_eatlas_accessions_file, checkIfExists: true) : Channel.empty()

    // getting accessions to exclude and preparing in the right format
    Channel.fromList( params.exclude_eatlas_accessions.tokenize(',') )
        .mix( ch_exclude_eatlas_accessions_file.splitText() )
        .unique()
        .map { it -> it.trim() }
        .toList()
        .map { lst -> [lst] } // list of lists : mandatory when combining in the next step
        .set { ch_excluded_accessions }

    // appending to accessions provided by the user
    // ensures that no accessions is present twice (provided by the user and fetched from E. Atlas)
    // removing E-PROT- accessions
    // removing excluded accessions
    ch_input_accessions
        .mix( ch_fetched_accessions.splitText() )
        .unique()
        .map { it -> it.trim() }
        .filter { it.startsWith('E-') && !it.startsWith('E-PROT-') }
        .combine ( ch_excluded_accessions )
        .filter { accession, excluded_accessions -> !(accession in excluded_accessions) }
        .map { accession, excluded_accessions -> accession }
        .set { ch_accessions }

    if ( !params.accessions_only ) {

        // Downloading Expression Atlas data for each accession in ch_accessions
        EXPRESSIONATLAS_GETDATA( ch_accessions )

        // adding dataset id (accession + data_type) in the file meta
        ch_etlas_design = addDatasetIdToMetadata( EXPRESSIONATLAS_GETDATA.out.design.flatten() )
        ch_eatlas_counts = addDatasetIdToMetadata( EXPRESSIONATLAS_GETDATA.out.counts.flatten() )

        // adding design files to the meta of their respective count files
        ch_eatlas_datasets = groupFilesByDatasetId( ch_etlas_design, ch_eatlas_counts )

        // adding normalisation state in the meta
        augmentToMetadata( ch_eatlas_datasets )

    }

    emit:
    downloaded_datasets = ch_eatlas_datasets

}



/*
========================================================================================
    FUNCTIONS
========================================================================================
*/


//
// Get Expression Atlas Batch ID (accession + data_type) from file stem
//
def addDatasetIdToMetadata( ch_files ) {
    return ch_files
            .map {
                file ->
                    def meta = [dataset: file.getSimpleName()]
                    [meta, file]
            }
}

//
// Groups design and data files by accession and data_type
// Design and count files have necessarily the same dataset ID (same file stem)
//
def groupFilesByDatasetId(ch_design, ch_counts) {
    return ch_design
        .concat( ch_counts ) // puts counts at the end of the resulting channel
        .groupTuple() // groups by dataset ID; design files are necessarily BEFORE count files
        .filter {
            it.get(1).size() == 2 // only groups with two files
        }
        .filter { // only groups with first file as design file and second one as count fileWARN: java.net.ConnectException: Connexion refusée
            meta, files ->
                files.get(0).name.endsWith('.design.csv') && !files.get(1).name.endsWith('.design.csv')
        }
        .map { // putting design file in meta
            meta, files ->
                def new_meta = meta + [design: files[0]]
                [new_meta, files[1]]
        }
}

def getNthPartFromEnd(String s, int n) {
    def tokens = s.tokenize('.')
    return tokens[tokens.size() - n]
}

//
// Add normalised: true / false in meta
//
def augmentToMetadata( ch_files ) {
    return ch_files
            .map {
                meta, file ->
                    if ( getNthPartFromEnd(file.name, 3) == 'raw' ) {
                        meta.normalised = false
                    } else {
                        meta.normalised = true
                    }
                    meta.platform = getNthPartFromEnd(file.name, 4)
                    [meta, file]
            }
}

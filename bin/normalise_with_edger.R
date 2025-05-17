#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

library(edgeR)
library(optparse)

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

get_args <- function() {

    option_list <- list(
        make_option("--counts", dest = 'count_file', help = "Path to input count file"),
        make_option("--design", dest = 'design_file', help = "Path to input design file")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Normalize counts using edgeR"
        ))

    return(args)
}

remove_all_zero_columns <- function(df) {
    # remove columns which contain only zeros
    df <- df[, colSums(df) != 0]
    return(df)
}

check_samples <- function(count_matrix, design_data) {
    # check if the column names of count_matrix match the sample names
    if (!all(colnames(count_matrix) == design_data$sample)) {
        stop("Sample names in the count matrix do not match the design data.")
    }
}

prefilter_counts <- function(count_matrix) {
    # remove genes having zeros for all counts
    # it is advised to remove them analysis
    non_zero_rows <- rownames(count_matrix[apply(count_matrix!=0, 1, any),])
    filtered_count_matrix <- count_matrix[rownames(count_matrix) %in% non_zero_rows, ]
    return(filtered_count_matrix)
}

replace_zero_counts_with_pseudocounts <- function(count_data_matrix) {
    # Add a small pseudocount of 0.01 to avoid zero counts
    count_data_matrix[count_data_matrix == 0] <- 0.01
    return(count_data_matrix)
}

filter_out_lowly_expressed_genes <- function(dge) {
    # filter the dataframe to exclude rows where the mean is 0
    # filter out lowly expressed genes
    keep <- filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    return(dge)
}


get_cpm_counts <- function(dge) {
    cpm_counts <- cpm(dge, normalised.lib.sizes = TRUE)
    return(cpm_counts)
}


get_normalised_cpm_counts <- function(count_file, design_file) {

    print(paste('Normalizing counts in:', count_file))

    count_data <- read.csv(args$count_file, row.names = 1)
    design_data <- read.csv(design_file)

    count_matrix <- as.matrix(count_data)
    # in some rare datasets, columns can contain only zeros
    # we do not consider these columns
    count_matrix <- remove_all_zero_columns(count_matrix)

    # getting design data
    design_data <- design_data[design_data$sample %in% colnames(count_matrix), ]

    # check if the column names of count_matrix match the sample names
    check_samples(count_matrix, design_data)

    # pre-filter genes with low counts
    count_matrix <- prefilter_counts(count_matrix)

    # Add a small pseudocount to avoid zero counts
    count_matrix_pseudocount <- replace_zero_counts_with_pseudocounts(count_matrix)

    group <- factor(design_data$condition)
    dge <- DGEList(counts = count_matrix_pseudocount, group = group)
    rownames(dge) <- rownames(count_matrix)
    colnames(dge) <- colnames(count_matrix)

    dge <- filter_out_lowly_expressed_genes(dge)

    # if the dataframe is now empty, stop the process
    if (nrow(dge) == 0) {
        message("No genes left after pre-filtering.")
        quit(save = "no", status = 100)
    }

    # normalisation
    dge <- calcNormFactors(dge, method="TMM")

    cpm_counts <- get_cpm_counts(dge)

    return(cpm_counts)
}

export_data <- function(cpm_counts, filename) {
    filename <- sub("\\.csv$", ".cpm.csv", filename)
    print(paste('Exporting normalised counts per million to:', filename))
    write.table(cpm_counts, filename, sep = ',', row.names = TRUE, col.names = NA, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

cpm_counts <- get_normalised_cpm_counts(args$count_file, args$design_file)

export_data(cpm_counts, basename(args$count_file))

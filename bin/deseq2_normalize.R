#!/usr/bin/env Rscript

library(DESeq2)
library(optparse)

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


get_args <- function() {

    option_list <- list(
        make_option("--counts", dest = 'count_file', help = "Path to input count file"),
        make_option("--design", dest = 'design_file', help = "Path to input design file"),
        make_option("--accession", help = "Accession number of expression atlas experiment. Example: E-MTAB-552")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Normalize counts using DESeq2"
        ))

    return(args)
}

get_normalized_counts <- function(count_file, design_file) {

    print(paste('Normalizing counts in:', count_file))

    count_data <- read.csv(count_file, row.names = 1)
    design_data <- read.csv(design_file)

    design_data <- design_data[design_data$sample %in% colnames(count_data), ]

    count_matrix <- as.matrix(count_data)

    col_data <- data.frame(
        row.names = design_data$sample,
        condition = factor(design_data$condition)
    )

    # Check if the column names of count_matrix match the row names of col_data
    if (!all(colnames(count_matrix) == rownames(col_data))) {
        stop("Sample names in the count matrix do not match the design data.")
    }

    # Add a small pseudocount of 1 (it has to be integer...) to avoid zero counts
    count_matrix[count_matrix == 0] <- 1

    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)

    # Filter genes with low counts
    # Keep genes that have at least 5 counts in at least two samples
    dds_filtered <- dds[rowSums(counts(dds) >= 5) >= 2, ]

    # perform normalization
    dds <- estimateSizeFactors(dds)

    normalized_counts <- counts(dds, normalized = TRUE)

    # Calculate total counts per sample (library size)
    library_sizes <- colSums(count_data)

    # Convert normalized counts to CPM
    cpm_counts <- t(t(normalized_counts) / library_sizes * 1e6)

    cpm_counts <- log2(cpm_counts + 1)

    return(cpm_counts)
}

export_data <- function(cpm_counts, accession) {
    filename <- paste0(accession, ".log_cpm.csv")
    print(paste('Exporting normalized counts per million to:', filename))
    write.table(cpm_counts, filename, sep = ',', row.names = TRUE, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

cpm_counts <- get_normalized_counts(args$count_file, args$design_file)

export_data(cpm_counts, args$accession)

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
        make_option("--design", dest = 'design_file', help = "Path to input design file")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Normalize counts using DESeq2"
        ))

    return(args)
}

get_variation_coefficient <- function(count_data, design_data) {

    print('Getting coefficients of variation')

    count_matrix <- as.matrix(count_data)
    col_data <- data.frame(
        row.names = design_data$sample,
        condition = factor(design_data$condition),
        batch = factor(design_data$batch)
    )

    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)

    # estimating size factors for normalization and dispersions
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    gene_disp <- dispersions(dds)

    # calculating the square root of the dispersion as a CV-like measure
    cv_df <- sqrt(gene_disp)

    return(cv_df)
}

export_data <- function(cv_df) {
    outfilename <- 'variation_coefficients.csv'
    print(paste('Exporting normalized counts to:', outfilename))
    write.table(cv_df, sep=",", file=outfilename, row.names = TRUE, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

count_data <- read.csv(args$count_file, header=TRUE, row.names = 1)
design_data <- read.csv(args$design_file, header=TRUE)

cv_df <- get_variation_coefficient(count_data, design_data)

export_data(cv_df)



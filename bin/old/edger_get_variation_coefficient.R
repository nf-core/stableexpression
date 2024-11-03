#!/usr/bin/env Rscript

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

get_variation_coefficient <- function(count_file, design_file) {

    print(paste('Normalizing counts in:', count_file))

    count_data <- read.csv(args$count_file, row.names = 1)
    design_data <- read.csv(design_file)

    design_data <- design_data[design_data$sample %in% colnames(count_data), ]

    group <- factor(design_data$group)

    count_data_matrix <- as.matrix(count_data)

    # Add a small pseudocount of 0.01 to avoid zero counts
    count_data_matrix[count_data_matrix == 0] <- 0.01

    dge <- DGEList(counts = count_data_matrix, group = group)

    rownames(dge) <- rownames(count_data_matrix)
    colnames(dge) <- colnames(count_data_matrix)

    # filter out lowly expressed genes
    keep <- filterByExpr(dge)
    ddge <- dge[keep, , keep.lib.sizes=FALSE]

    #normalization
    dge <- calcNormFactors(dge, method="TMM")
    normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)

    # estimating dispersion
    dge <- estimateCommonDisp(dge)
    dge <- estimateTagwiseDisp(dge)

    # getting the coefficient of variation as the sqrt of the tagwise dispersion
    cv_df <- data.frame(
        gene = rownames(dge),
        variation_coefficient = sqrt(dge$tagwise.dispersion),
        average_log_cpm = dge$AveLogCPM
    )

    cv_df <- cv_df[order(cv_df$variation_coefficient, decreasing = FALSE), ]

    return(cv_df, normalized_counts)
}

export_data <- function(cv_df, normalized_counts) {
    cv_filename <- "variation_coefficient.csv"
    print(paste('Exporting coefficients of variation to:', cv_filename))
    write.table(cv_df, cv_filename, sep = ',', row.names = TRUE, quote = FALSE)

    normalized_filename <- "all_counts.normalized.csv"
    print(paste('Exporting normalized counts to:', normalized_filename))
    write.table(normalized_counts, normalized_filename, sep = ',', row.names = TRUE, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

cv_df, normalized_counts <- get_variation_coefficient(args$count_file, args$design_file)

export_data(cv_df, normalized_counts)

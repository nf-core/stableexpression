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

get_normalized_counts <- function(count_file, design_file) {

    print(paste('Normalizing counts in:', count_file))

    count_data <- read.csv(args$count_file, row.names = 1)
    design_data <- read.csv(design_file)

    design_data <- design_data[design_data$sample %in% colnames(count_data), ]

    group <- factor(design_data$condition)

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

    cpm_counts <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)

    return(cpm_counts)
}

export_data <- function(cpm_counts, filename) {
    filename <- sub("\\.csv$", ".log_cpm.csv", filename)
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

export_data(cpm_counts, basename(args$count_file))

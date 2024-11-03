#!/usr/bin/env Rscript

library(optparse)

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


get_args <- function() {

    option_list <- list(
        make_option("--counts", dest = 'count_file', help = "Path to input count file (normalized)")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Get variation coefficient from count data for each gene"
        ))

    return(args)
}

average_log2 <- function(row) {
    # the dataframe has already been filtered to exclude rows where mean is 0
    return(mean(log2(row + 1))) # adds 1 to avoid log(0) and to stabilize variance
}

get_variation_coefficient <- function(count_data) {

    print('Getting coefficients of variation')

    # filter the dataframe to exclude rows where the mean is 0
    count_data <- count_data[rowMeans(count_data) != 0, ]

    # filter the dataframe to exclude rows where row mean is in the top 5% or bottom 5%
    # determine the percentile thresholds
    row_means <- rowMeans(count_data)
    lower_threshold <- quantile(row_means, 0.05)
    upper_threshold <- quantile(row_means, 0.95)

    count_data <- count_data[row_means >= lower_threshold & row_means <= upper_threshold, ]

    # calculate the coefficient of variation
    row_means <- rowMeans(count_data)
    row_sds <- apply(count_data, 1, sd)
    cv <- row_sds / row_means

    av_log_cpm <- apply(count_data, 1, average_log2)
    # combine results into a dataframe
    df <- data.frame(
        gene = rownames(count_data),
        variation_coefficient = cv,
        average_log_cpm = av_log_cpm
        )

    df <- df[order(df$variation_coefficient, decreasing = FALSE), ]

    return(df)
}

export_data <- function(cv_df) {
    outfilename <- 'variation_coefficients.csv'
    print(paste('Exporting normalized counts to:', outfilename))
    write.table(cv_df, sep=",", file=outfilename, row.names = FALSE, col.names = TRUE, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

count_data <- read.csv(args$count_file, header=TRUE, row.names = 1)

cv_df <- get_variation_coefficient(count_data)

export_data(cv_df)

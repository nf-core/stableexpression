#!/usr/bin/env Rscript

library("DESeq2")

count_file <- '$count_file'
count_data <- read.csv(count_file, row.names = 1)

count_data_matrix <- as.matrix(count_data)
colData <- DataFrame(row.names = colnames(count_data_matrix))

# Add a small pseudocount of 1 (it has to be integer...) to avoid zero counts
count_data_matrix[count_data_matrix == 0] <- 1

dds <- DESeqDataSetFromMatrix(countData = count_data_matrix, colData = colData, design = ~ 1)

# perform normalization
dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)

normalized_filename <- sub(".csv", "_normalized.csv", basename(count_file))
write.csv(normalized_counts, normalized_filename, row.names = TRUE)

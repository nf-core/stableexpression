#!/usr/bin/env Rscript

library(edgeR)

count_file <- '$count_file'

count_data <- read.csv(count_file, header=TRUE, row.names = 1)
count_data_matrix <- as.matrix(count_data)

# Add a small pseudocount of 0.01 to avoid zero counts
count_data_matrix[count_data_matrix == 0] <- 0.01

dge <- DGEList(counts = count_data_matrix)
rownames(dge) <- rownames(count_data_matrix)
colnames(dge) <- colnames(count_data_matrix)

#normalization
dge <- calcNormFactors(dge, method="TMM")
normalized_counts <- cpm(dge)

normalized_filename <- sub(".csv", "_normalized.csv", basename(count_file))
write.table(normalized_counts, normalized_filename,, row.names = TRUE, quote = FALSE)

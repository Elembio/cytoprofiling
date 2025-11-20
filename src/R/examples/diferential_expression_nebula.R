library(cytoprofiling)
library(nebula)

calculate_scaling_factors <- function(df, normalization_targets) {
    normed_df = df
    norm_values <- rowSums(normed_df[, normalization_targets], na.rm = TRUE)
    for (target in normalization_targets) {
        normed_df[target] <- normed_df[target] / norm_values
    }
    
    # aggregate by mean count per well, ignoring NAs
    well_counts <- aggregate(normed_df[, normalization_targets], by = list(df$Well), FUN = mean, na.rm = TRUE)
    
    # log transform well_counts
    well_counts[, -1] <- log(well_counts[, -1])
    
    # center each column of well_counts by mean
    well_counts[, -1] <- sweep(well_counts[, -1], 2, colMeans(well_counts[, -1], na.rm = TRUE))
    
    # find the median of each row of well_counts
    median_counts <- apply(well_counts[, -1], 1, median, na.rm = TRUE)
    
    # calculate the scaling factor as the exponential of the median log ratio
    scaling_factors <- exp(median_counts)
    
    return(scaling_factors)
}

parquet_file = ""
apply_library_size_correction = TRUE

# load and filter the data the data
df <- load_cytoprofiling(parquet_file)
df <- filter_cells(df)

well_names <- unique(df$Well)
batch_names <- get_barcoding_batches(df)

# Initialize empty data frame for results
final_results <- data.frame()
for (batch_name in batch_names){
        print(batch_name)

        normalization_targets <- get_default_normalization_targets(df, batch_name)

        # sum the counts per each cell across the normalization targets
        cell_sums <- as.vector(rowSums(df[, normalization_targets], na.rm = TRUE))

        if (apply_library_size_correction) {
            adjusted_cell_sums <- cell_sums * calculate_scaling_factors(df, normalization_targets)[match(df$Well, well_names)]
        } else {
            adjusted_cell_sums <- cell_sums
        }

        pred <- model.matrix(~ WellLabel, data = df)

        # Create the count matrix
        count_matrix <- t(as.matrix(df[normalization_targets]))

        # Fit model
        res <- nebula(count = count_matrix,
                    id = df$Well,
                    pred = pred,
                    offset = adjusted_cell_sums)

        # Add batch information to results
        res$summary$batch <- batch_name
        
        # Append to final results
        final_results <- rbind(final_results, res$summary)
}

# Write final summary data frame to CSV file
write.csv(final_results, file = "nebula_output.csv", row.names = FALSE, quote = FALSE)
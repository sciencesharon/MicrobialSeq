#this function takes the input directly from a metaphlan4 merged abundance table output 
#it generates read count tables at all taxanomic levels - Kingdom, Phylum, Class, Order, Family, Genus, Species, and Strain
#it automatically removes unmapped and unknown taxa levels
#it also agglomerates the taxa so that none are repeated, summing the read counts per identical taxa
#finally, it calculates relative abundance and outputs relabundance tables at all taxanomic levels
#each output is saved to the given directory and also output to the global environment

merged_abundance_processing <- function(merged, output_dir) {
  # Load required packages
  required_packages <- c("dplyr", "tidyr")
  lapply(required_packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  })
  
  # Step 1: Make a temporary tax table 
  merged_abundance_to_tax <- function(merged) {
    tax <- as.data.frame(merged[, c(1)])
    tax <- as.data.frame(tax[-c(1:2), ])
    rownames(tax) <- NULL
    rownames(tax) <- sub("^", "OTU", rownames(tax))
    colnames(tax)[1] <- "V1"
    tax <- tax %>% separate(V1, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), "\\|")
    tax$Kingdom <- gsub("k__", "", tax$Kingdom)
    tax$Phylum <- gsub("p__", "", tax$Phylum)
    tax$Class <- gsub("c__", "", tax$Class)
    tax$Order <- gsub("o__", "", tax$Order)
    tax$Family <- gsub("f__", "", tax$Family)
    tax$Genus <- gsub("g__", "", tax$Genus)
    tax$Species <- gsub("s__", "", tax$Species)
    tax$Strain <- gsub("t__", "", tax$Strain)
    
    tax$Species <- gsub("_", " ", tax$Species)
    tax$Strain <- gsub("_", " ", tax$Strain)
    
    return(as.data.frame(tax))
  }
  
  # Step 2: Make a temporary OTU table
  merged_abundance_to_otu <- function(merged) {
    colnames(merged) <- merged[1, ]
    merged <- merged[-c(1:2), ]
    merged <- merged[!grepl('Archaea|Viruses|Eukaryota', merged$clade_name), ]
    otu <- merged[, -c(1)]
    rownames(otu) <- NULL
    rownames(otu) <- sub("^", "OTU", rownames(otu))
    otu <- dplyr::mutate_all(otu, as.numeric)
    return(as.matrix(otu))
  }
  
  # Step 3: Subset the tax table by !is.na for each column and save as a new dataframe
  subset_by_na <- function(df) {
    df_list <- list()
    for (col_name in names(df)) {
      subset_df <- df[!is.na(df[[col_name]]), col_name, drop = FALSE]
      assign(col_name, subset_df, envir = .GlobalEnv)
      df_list[[col_name]] <- subset_df
    }
    return(df_list)
  }
  
  # Step 4: Subset OTU table by rownames of each new tax table 
  subset_by_rownames <- function(df_list, dataframe2) {
    df_subset_list <- list()
    for (df_name in names(df_list)) {
      subset_df <- df_list[[df_name]]
      subset_dataframe2 <- dataframe2[rownames(dataframe2) %in% rownames(subset_df), ]
      if (!all(rownames(subset_dataframe2) %in% rownames(subset_df))) {
        warning(paste("Row names do not match for", df_name))
      }
      assign(df_name, subset_dataframe2, envir = .GlobalEnv)
      df_subset_list[[df_name]] <- subset_dataframe2
    }
    return(df_subset_list)
  }
  
  # Step 5: Now we combine the otu and tax tables 
  subset_and_cbind <- function(df_list, dataframe2) {
    df_combined_list <- list()
    for (df_name in names(df_list)) {
      subset_df <- df_list[[df_name]]
      subset_dataframe2 <- dataframe2[rownames(dataframe2) %in% rownames(subset_df), ]
      if (!all(rownames(subset_dataframe2) %in% rownames(subset_df))) {
        warning(paste("Row names do not match for", df_name))
      }
      combined_df <- cbind(subset_df, subset_dataframe2[rownames(subset_df), , drop = FALSE])
      assign(df_name, combined_df, envir = .GlobalEnv)
      df_combined_list[[df_name]] <- combined_df
    }
    return(df_combined_list)
  }
  
  # Step 6: Combine rows which are identical taxa 
  summarize_and_save_to_csv <- function(df_list, output_dir) {
    summarized_list <- list()
    for (df_name in names(df_list)) {
      df <- df_list[[df_name]]
      df <- df[!grepl('Unclassified|GB|unclassified', df[[1]], ignore.case = TRUE), ]
      unique_values <- unique(df[[1]])
      result_df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
      colnames(result_df) <- colnames(df)
      for (value in unique_values) {
        subset_df <- df[df[[1]] == value, ]
        summed_row <- colSums(subset_df[, -1], na.rm = TRUE)
        new_row <- data.frame(value, t(summed_row))
        colnames(new_row) <- colnames(df)
        result_df <- rbind(result_df, new_row)
      }
      rownames(result_df) <- result_df[[1]]
      result_df <- result_df[, -1]
      new_df_name <- paste0(df_name, "_ReadCount")
      assign(new_df_name, result_df, envir = .GlobalEnv)
      summarized_list[[new_df_name]] <- result_df
      write.csv(result_df, file.path(output_dir, paste0(new_df_name, ".csv")), row.names = TRUE)
    }
    return(summarized_list)
  }
  
  # Step 7: Convert to Relative Abundance
  convert_to_percentages_list <- function(df_list, output_dir) {
    df_percent_list <- list()
    for (df_name in names(df_list)) {
      df <- df_list[[df_name]]
      df_percent <- df
      for (col_name in colnames(df)) {
        col_sum <- sum(df[[col_name]], na.rm = TRUE)
        df_percent[[col_name]] <- (df[[col_name]] / col_sum) * 100
      }
      new_df_name <- paste0(sub("_.*", "", df_name), "_RelAbun")
      assign(new_df_name, df_percent, envir = .GlobalEnv)
      df_percent_list[[new_df_name]] <- df_percent
      write.csv(df_percent, file.path(output_dir, paste0(new_df_name, ".csv")), row.names = TRUE)
    }
    return(df_percent_list)
  }
  
  tax <- merged_abundance_to_tax(merged)
  otu <- merged_abundance_to_otu(merged)
  df_list <- subset_by_na(tax)
  df_list2 <- subset_by_rownames(df_list, otu)
  df_combined_list <- subset_and_cbind(df_list, otu)
  summarized_list <- summarize_and_save_to_csv(df_combined_list, output_dir)
  percentage_list <- convert_to_percentages_list(summarized_list, output_dir)
  
  return(percentage_list)
}

#example usage
processed_list <- merged_abundance_processing(merged, "/home/cornell/Documents/Temp")

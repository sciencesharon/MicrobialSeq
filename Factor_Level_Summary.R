#this function takes the input of a metadata table with samples as rows and columns of character data
#the output is a summary table of the number of samples in each factor level per column
#it is formatted for ease of use in publication
#it also compares the groups of samples by fishers exact test to determine if there are differing numbers of samples in each group
#it is useful for making a demographic table for publication
#comparing experimental groups and their other demographic characteristics to find differences between groups

factor_levels_summaries <- function(datasets, dataset_names, columns, file_path) {
  results_list <- list()
  
  for (column in columns) {
    # Ensure the columns are factors
    datasets <- map(datasets, ~mutate(.x, across(all_of(column), as.factor)))
    
    # Get the summary of factor levels and their counts for each dataset
    summaries <- map(seq_along(datasets), function(i) {
      datasets[[i]] %>%
        group_by(across(all_of(column))) %>%
        summarise(!!dataset_names[i] := n(), .groups = 'drop') %>%
        rename_with(~ "FactorLevel", all_of(column))
    })
    
    # Merge all summaries
    combined <- reduce(summaries, function(x, y) dplyr::full_join(x, y, by = "FactorLevel"))
    
    # Fill missing values with 0
    combined <- combined %>%
      mutate(across(all_of(dataset_names), ~replace_na(.x, 0))) %>%
      mutate(Group = column)
    
    # Ensure all count columns are numeric
    combined <- combined %>%
      mutate(across(all_of(dataset_names), as.numeric))
    
    # Perform Fisher's exact test for each dataset pair
    fisher_p_values <- map_dbl(seq_len(nrow(combined)), function(i) {
      p_values <- NULL
      for (j in seq(2, ncol(combined) - 2)) {
        for (k in seq(j + 1, ncol(combined) - 2)) {
          if (startsWith(names(combined)[j], dataset_names[j-1]) && startsWith(names(combined)[k], dataset_names[k-1])) {
            count1 <- combined[[j]][i]
            count2 <- combined[[k]][i]
            if (all(is.finite(count1), is.finite(count2), count1 >= 0, count2 >= 0)) {
              p_value <- fisher.test(matrix(c(count1, sum(combined[[j]]) - count1,
                                              count2, sum(combined[[k]]) - count2),
                                            nrow = 2))$p.value
              p_values <- c(p_values, p_value)
            } else {
              p_values <- c(p_values, NA)
            }
          }
        }
      }
      min(p_values, na.rm = TRUE)  # Return the smallest p-value for significance
    })
    
    combined <- combined %>%
      mutate(p_value = format(fisher_p_values, digits = 2, scientific = TRUE),
             Significant = ifelse(as.numeric(fisher_p_values) < 0.05, "*", ""),
             p_value_combined = if_else(is.na(p_value), "", paste0(p_value, Significant)))
    
    # Debugging print statements
    print(paste("Combined after adding p_value_combined for column:", column))
    print(colnames(combined))
    print(head(combined))
    
    # Add the column name as a row in the FactorLevel column with matching column names
    new_row <- tibble(FactorLevel = column)
    combined <- bind_rows(new_row, combined)
    
    # Remove the Group, p_value, and Significant columns using base R
    combined <- combined[, !names(combined) %in% c("Group", "p_value", "Significant")]
    
    # Debugging print statements
    print(paste("Combined after removing columns for column:", column))
    print(colnames(combined))
    print(head(combined))
    
    # Add the combined summary for this column to the results list
    results_list[[column]] <- combined
  }
  
  # Combine all dataframes in the list into a single dataframe
  final_combined <- bind_rows(results_list)
  
  # Debugging print statements
  print("Final combined before ensuring p_value_combined:")
  print(colnames(final_combined))
  print(head(final_combined))
  
  # Ensure the p_value_combined column is created
  if (!"p_value_combined" %in% colnames(final_combined)) {
    final_combined$p_value_combined <- ""
  }
  
  # Debugging print statements
  print("Final combined after ensuring p_value_combined:")
  print(colnames(final_combined))
  print(head(final_combined))
  
  # Replace NAs with empty strings in specific columns
  final_combined <- final_combined %>%
    mutate(across(all_of(dataset_names), ~replace_na(as.character(.x), "")))
  
  # Replace NAs with empty strings in the p_value column using base R
  final_combined[is.na(final_combined)] <- ""
  
  # Debugging print statements
  print("Final combined before writing to CSV:")
  print(colnames(final_combined))
  print(head(final_combined))
  
  # Rename p_value_combined to p_value using base R
  colnames(final_combined)[colnames(final_combined) == "p_value_combined"] <- "p_value"
  
  # Debugging print statements
  print("Final combined after renaming p_value_combined to p_value:")
  print(colnames(final_combined))
  print(head(final_combined))
  
  # Write the final dataframe to a .csv file
  write.csv(final_combined, file_path, row.names = FALSE)
  
  return(final_combined)
} 

metadata <- read.csv("/path/to/your/file/Example_Metadata.csv")

# Convert specified columns to factors
character_cols <- colnames(metadata)[1:5]
metadata <- metadata %>% mutate(across(all_of(character_cols), as.factor))

#separate out groups of metadata which you would like to compare
Cherry <- metadata[metadata$Fruit %in% c("Cherry"),]
Banana <- metadata[metadata$Fruit %in% c("Banana"),]
Apple <- metadata[metadata$Fruit %in% c("Apple"),]
Date <- metadata[metadata$Fruit %in% c("Date"),]
Elderberry <- metadata[metadata$Fruit %in% c("Elderberry"),]

#list the dataframes
datasets <- list(Cherry, Banana, Apple, Date, Elderberry) 

#set the names of the dataframes
dataset_names <- c("Cherry", "Banana", "Apple", "Date", "Elderberry")

#set the destination file path and file name
file_path <- ("/path/to/summary.csv")

#set the columns to compare
columns <- colnames(metadata)[1:5]

#run the function
summary <- factor_levels_summaries(datasets = datasets, 
                                   dataset_names = dataset_names, 
                                   columns = columns, 
                                   file_path = file_path)


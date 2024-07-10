#This function takes input of samples in rows and numerical data in columns
#You also need character data column(s) by which to separate and compare the data
#For example, a column of sample gender with levels Male and Female
#Or experimental groups like control and variable1, variable2, variable3 with numerical data you wish to compare 
#This function makes a nice little table with the mean of a numerical column input for each group and the standard deviation 
#You can also optionally input either "ttest" or "wilcox" to perform statistical testing - this outputs the p-value
#You can also optionally adjust the p-value by either bonferroni or BH - this outputs the p-adjust 
#Then finally you can also optionally add a column called "stars" which will show significance stars based on the p-value or p-adjust 

data <- read.csv("/path/to/your/data/Example_Data_Table.csv")

Male <- data[data$Gender %in% c("Male"),]
Female <- data[data$Gender %in% c("Female"),]

columns <- colnames(data)[1:10]

data_list <- list(Male = Male, Female = Female)

calculate_mean_sd <- function(data_list, columns, file_name, directory, decimal_places = 2, statistical_test = NULL, p_adjust_method = NULL, stars = FALSE) {
  # Initialize an empty data frame to store results
  results <- data.frame(
    Column = columns,
    stringsAsFactors = FALSE
  )
  
  # Iterate over each data frame in the list
  for (name in names(data_list)) {
    data <- data_list[[name]]
    
    # Initialize an empty vector to store the mean ± SD strings
    mean_sd_vec <- character(length(columns))
    
    for (i in seq_along(columns)) {
      column <- columns[i]
      
      if (!is.numeric(data[[column]])) {
        stop(paste(column, "is not numeric. Please provide numeric columns."))
      }
      
      mean_val <- mean(data[[column]], na.rm = TRUE)
      sd_val <- sd(data[[column]], na.rm = TRUE)
      mean_sd_vec[i] <- paste0(round(mean_val, decimal_places), " ± ", round(sd_val, decimal_places))
    }
    
    # Add the mean ± SD strings as a new column in the results data frame
    results[[name]] <- mean_sd_vec
  }
  
  # Perform statistical tests if requested
  if (!is.null(statistical_test) && length(data_list) == 2) {
    p_values <- numeric(length(columns))
    
    for (i in seq_along(columns)) {
      column <- columns[i]
      group1 <- data_list[[1]][[column]]
      group2 <- data_list[[2]][[column]]
      
      if (is.numeric(group1) && is.numeric(group2)) {
        if (statistical_test == "ttest") {
          test_result <- t.test(group1, group2, var.equal = TRUE)
        } else if (statistical_test == "wilcox") {
          test_result <- wilcox.test(group1, group2)
        } else {
          stop("Unsupported statistical test. Choose 'ttest' or 'wilcox'.")
        }
        p_values[i] <- test_result$p.value
      } else {
        p_values[i] <- NA
      }
    }
    
    # Adjust p-values if requested
    if (!is.null(p_adjust_method)) {
      p_values <- p.adjust(p_values, method = p_adjust_method)
      results$P_Adjust <- ifelse(p_values < 0.05,
                                 formatC(p_values, format = "e", digits = decimal_places),
                                 formatC(p_values, format = "f", digits = decimal_places))
    } else {
      results$P_Value <- ifelse(p_values < 0.05,
                                formatC(p_values, format = "e", digits = decimal_places),
                                formatC(p_values, format = "f", digits = decimal_places))
    }
    
    # Add significance stars if requested
    if (stars) {
      significance_stars <- character(length(p_values))
      
      for (i in seq_along(p_values)) {
        p_val <- p_values[i]
        if (p_val < 0.000005) {
          significance_stars[i] <- "*****"
        } else if (p_val < 0.00005) {
          significance_stars[i] <- "****"
        } else if (p_val < 0.0005) {
          significance_stars[i] <- "***"
        } else if (p_val < 0.005) {
          significance_stars[i] <- "**"
        } else if (p_val < 0.05) {
          significance_stars[i] <- "*"
        } else {
          significance_stars[i] <- ""
        }
      }
      
      results$Significance <- significance_stars
    }
  } else if (!is.null(statistical_test) && length(data_list) != 2) {
    stop("Statistical tests can only be performed when exactly two groups are provided.")
  }
  
  # Construct the full file path
  file_path <- file.path(directory, file_name)
  
  # Save the result as a CSV file
  write.csv(results, file = file_path, row.names = FALSE)
  
  return(results)
}

dir <- ("/path/to/save/your/output")

means <- calculate_mean_sd(data_list, 
                           columns, 
                           file_name = "mean_SD.csv", 
                           directory = dir, 
                           decimal_places = 2,
                           statistical_test = "ttest",
                           p_adjust_method = "BH",
                           stars = TRUE)


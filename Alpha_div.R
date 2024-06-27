#this function takes the input of a count data table
#the output of readcount table from Merged_Abundance_Processing.R can be used directly with this function
#count_data is the count data table in a dataframe format with samples as columns and taxa as rows
#met is the metadata file with the grouping/comparision variable
#the metadata file needs to have samples in rows and groups as columns
#the samples as rows in the metadata must match the columns in the Read Count table
#the comparison_var is the grouping variable by which the data will be compared
#the output_dir is where the data and the pdf file plots will be saved
#the function automatically does chao1, eveness, richness, and shannon diversity
#the output is four plots, the alpha diversity calculations as a dataframe .csv file, and the wilcox test results for each comparison
#the function also does the statistical comparison of groups by wilcox rank sum
#it plots the p-value by stars for significant comparisons and leaves out NS comparisons
#stack_plots allows for optional output of stacked richness, evenness, and complexity plot output
#you can also change the width and height of the output pdf files

alpha_div_auto <- function(count_data, met, comparison_var, output_dir = "results", pdf_width = 8, pdf_height = 5, stack_plots = FALSE) {
  
  # Check if rownames of met match column names of count_data
  if (!all(rownames(met) %in% colnames(count_data))) {
    stop("Rownames of met do not match column names of count_data")
  }
  
  # Transform metadata columns into factors
  met <- met %>% mutate(across(everything(), as.factor))
  
  # Alpha diversity calculations
  alpha_div <- function(count_data, met) {
    div_x <- diversity(count_data, index = "shannon")
    count_data <- as.data.frame(t(count_data))
    count_data <- round_df(count_data, 0, rf = "round")
    richness <- estimateR(count_data)
    shannon <- diversity(t(count_data), index = "shannon")
    evenness <- div_x / log(specnumber(count_data))
    colnames(evenness) = "evenness"
    alphadiv <- cbind(t(richness), shannon, evenness)
    alphadiv <- alphadiv %>% mutate_if(is.character, as.numeric)
    alphadiv <- cbind(met, alphadiv)
    return(as.data.frame(alphadiv))
  }
  
  alpha_data <- alpha_div(count_data, met)
  
  # Save alpha_data to global environment and as a CSV file
  assign("alpha_data", alpha_data, envir = .GlobalEnv)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  write.csv(alpha_data, file.path(output_dir, "alpha_data.csv"), row.names = FALSE)
  
  # Perform pairwise Wilcoxon tests and save results
  perform_pairwise_tests <- function(data, variable, xy, output_dir = "results") {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    
    perform_and_save <- function(data, formula, file_name) {
      result <- data %>% pairwise_wilcox_test(formula, p.adjust.method = "BH")
      result <- result %>% add_xy_position(x = xy)
      result <- result %>% mutate(across(where(is.list), ~ sapply(., toString)))
      write.csv(result, file.path(output_dir, file_name), row.names = FALSE)
      return(result)
    }
    
    pwc_shannon <- perform_and_save(data, reformulate(variable, "shannon"), "pwc_shannon.csv")
    pwc_chao <- perform_and_save(data, reformulate(variable, "S.chao1"), "pwc_chao.csv")
    pwc_eveness <- perform_and_save(data, reformulate(variable, "evenness"), "pwc_eveness.csv")
    pwc_richness <- perform_and_save(data, reformulate(variable, "S.obs"), "pwc_richness.csv")
    
    return(list(
      pwc_shannon = pwc_shannon,
      pwc_chao = pwc_chao,
      pwc_eveness = pwc_eveness,
      pwc_richness = pwc_richness
    ))
  }
  
  pairwise_results <- perform_pairwise_tests(alpha_data, comparison_var, comparison_var, output_dir)
  
  # Function to create individual plots
  plot_single_metric <- function(data, y_var, x_var, fill_var, plot_title, y_label, plot_tag) {
    pwc <- data %>% pairwise_wilcox_test(as.formula(paste(y_var, "~", x_var)), p.adjust.method = "BH")
    pwc <- pwc %>% add_xy_position(x = x_var)
    
    p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
      geom_boxplot(aes_string(fill = fill_var)) +
      labs(title = plot_title, x = ' ', y = y_label, tag = plot_tag) +
      scale_fill_grey(guide = "none") +
      theme_classic() +
      stat_pvalue_manual(pwc, label = "p.adj.signif", tip.length = 0, step.increase = 0.1, hide.ns = TRUE)
    
    return(p)
  }
  
  # Create and save plots
  create_pairwise_plots <- function(data, formula, x_var, fill_var, output_dir, width = 8, height = 5, stack_plots = FALSE) {
    plot_and_save <- function(data, y_var, x_var, fill_var, output_file, plot_title, y_label, plot_tag) {
      if (!x_var %in% colnames(data) || !y_var %in% colnames(data) || !fill_var %in% colnames(data)) {
        stop(paste("Necessary columns not found in the data"))
      }
      
      p <- plot_single_metric(data, y_var, x_var, fill_var, plot_title, y_label, plot_tag)
      
      ggsave(output_file, plot = p, width = width, height = height)
    }
    
    if (stack_plots) {
      plot1 <- plot_single_metric(data, "S.chao1", x_var, fill_var, "Richness", "Chao1", "A")
      plot2 <- plot_single_metric(data, "evenness", x_var, fill_var, "Evenness", "Evenness", "B")
      plot3 <- plot_single_metric(data, "shannon", x_var, fill_var, "Complexity", "Shannon", "C")
      
      combined_plot <- plot1 / plot2 / plot3
      
      ggsave(file.path(output_dir, "combined_plot.pdf"), plot = combined_plot, width = width, height = height * 3)
    } else {
      metrics <- list(
        list(y_var = "S.obs", plot_title = "Richness", y_label = "Richness", plot_tag = "A", file_name = "richness.pdf"),
        list(y_var = "S.chao1", plot_title = "Chao1", y_label = "Chao1 Index", plot_tag = "B", file_name = "chao1.pdf"),
        list(y_var = "evenness", plot_title = "Evenness", y_label = "Evenness", plot_tag = "C", file_name = "evenness.pdf"),
        list(y_var = "shannon", plot_title = "Shannon", y_label = "Shannon Diversity Index", plot_tag = "D", file_name = "shannon.pdf")
      )
      
      for (metric in metrics) {
        tryCatch({
          plot_and_save(data, metric$y_var, x_var, fill_var, file.path(output_dir, metric$file_name), 
                        metric$plot_title, metric$y_label, metric$plot_tag)
        }, error = function(e) {
          print(paste("Error creating plot for", metric$y_var, ":", e$message))
        })
      }
    }
  }
  
  create_pairwise_plots(alpha_data, paste("~", comparison_var), comparison_var, comparison_var, output_dir, pdf_width, pdf_height, stack_plots)
}

#this function takes a list of files as output by metaphlan4 diversity.R function
#the Bray-curtis, Jaccard, Aitchison, Unifrac, and Weighted Unifrac files directly output by the diversity.R function
#can be input directly to this function
#you can optionally subset the beta diversity table by the samples using subset_names
#prefix is the experiment prefix which will be added to the beginning of each output
#the initial output if plot = FALSE is the subsetted files to your output directory
#if plot = TRUE, the function will output PCoA plots of the files 
#this requires input of metadata (met)
#the metadata should have the samples in rows and the grouping variable(s) as columns
#the PcOA plots generated will have ellipses aroun your grouping variable (color and fill)
#to generate the PCoA plot, you need to input the colour and fill - these are each the grouping variable (one of the columns in the metadata)
#the title of the plot is optionally set by title
#and set the plot height and plot width for the dimensions of the output pdfs
#you can optionally change the color of the ellipses around your groups with fill_palette 

beta_div_auto <- function(file_paths, output_directory, prefix, subset_names = NULL, met, colour, fill, title = NULL, plot_width = 8, plot_height = 6, fill_palette = NULL, plot = TRUE) {
  for (file_path in file_paths) {
    # Read the first line to get the column names
    header <- read_tsv(file_path, col_names = FALSE, n_max = 1)
    header <- unlist(header)
    
    # Shift column names and add "SampleID" as the first column name
    col_names <- c("SampleID", header)
    
    # Read the TSV file again without column names, starting from the second line
    data <- read_tsv(file_path, col_names = FALSE, skip = 1)
    
    # Apply the adjusted column names
    colnames(data) <- col_names
    
    # Convert the "SampleID" column to row names
    data <- column_to_rownames(data, var = "SampleID")
    
    # If subset_names is provided, subset the data frame
    if (!is.null(subset_names)) {
      data <- data[subset_names, subset_names, drop = FALSE]
    }
    
    # Convert the data frame to a matrix
    data_matrix <- as.matrix(data)
    
    # Extract the original filename without extension
    original_filename <- file_path_sans_ext(basename(file_path))
    
    # Construct the output file path with prefix for the CSV
    output_file_name <- paste0(prefix, "_", original_filename, ".csv")
    output_file_path <- file.path(output_directory, output_file_name)
    
    # Write the matrix to a CSV file in the specified output directory
    write.csv(data_matrix, file = output_file_path, row.names = TRUE)
    
    # Assign the matrix to the global environment with the name prefix_filename
    assign(paste0(prefix, "_", original_filename), data_matrix, envir = .GlobalEnv)
    
    # If plotting is not enabled, stop here
    if (!plot) {
      next
    }
    
    # Perform PCoA
    pcoa_result <- pcoa(data_matrix)
    
    # Calculate percentage values
    Axis1.percent <- pcoa_result$values$Relative_eig[[1]] * 100
    Axis2.percent <- pcoa_result$values$Relative_eig[[2]] * 100
    
    # Create the PCoA data frame
    pcoa_data <- data.frame(
      Sample = rownames(pcoa_result$vectors),
      X = pcoa_result$vectors[, 1],
      Y = pcoa_result$vectors[, 2]
    )
    
    # Combine with metadata
    pcoa_data <- cbind(pcoa_data, met)
    
    # Save the PCoA data frame as a CSV file
    pcoa_csv_name <- paste0(prefix, "_pcoa_", original_filename, ".csv")
    pcoa_csv_path <- file.path(output_directory, pcoa_csv_name)
    write.csv(pcoa_data, pcoa_csv_path, row.names = FALSE)
    
    # Assign the PCoA data to the global environment with the name prefix_filename_pcoa_data
    assign(paste0(prefix, "_", original_filename, "_pcoa_data"), pcoa_data, envir = .GlobalEnv)
    
    # Set fill and color scales based on the provided fill_palette or default to discrete colors
    if (is.null(fill_palette)) {
      fill_colors <- scale_fill_discrete()
      color_scale <- scale_color_discrete()
    } else {
      fill_colors <- scale_fill_manual(values = fill_palette)
      color_scale <- scale_color_manual(values = fill_palette)
    }
    
    # Generate the ggplot
    p <- ggplot(data = pcoa_data, aes(x = X, y = Y)) +
      geom_point(aes_string(colour = colour), size = 1.5) +
      xlab(paste("PCoA1 - ", round(Axis1.percent, 2), "%", sep = "")) +
      ylab(paste("PCoA2 - ", round(Axis2.percent, 2), "%", sep = "")) +
      theme_bw() +
      theme(
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 10),
        legend.key.size = unit(1, 'lines')
      ) +
      stat_ellipse(geom = "polygon", aes_string(fill = fill), alpha = 0.25) +
      fill_colors + color_scale
    
    # Add the title if provided
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    }
    
    # Construct the output file path with prefix for the PDF
    pdf_file_name <- paste0(prefix, "_", original_filename, ".pdf")
    pdf_file_path <- file.path(output_directory, pdf_file_name)
    
    # Save the plot as a PDF
    ggsave(pdf_file_path, plot = p, width = plot_width, height = plot_height)
  }
}

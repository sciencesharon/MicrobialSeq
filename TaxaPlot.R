# Install required packages
install.packages(c("dplyr", "tidyr", "reshape2", "tibble", "ggplot2", "stringr", "SQMtools"))

# Load required packages
library(dplyr)
library(tidyr)
library(reshape2)
library(tibble)
library(ggplot2)
library(stringr)
library(SQMtools)

#read your merged abundance table
merged <- read.table("/path/to/file/merged_abundance_example.txt", header = FALSE, sep = "\t") 


#set output directory
output_dir <- ("/path/to/file/output")

#utilize merged_abundance_processing.R
relabun_list <- merged_abundance_processing(merged, output_dir)

#read your metadata
met <- read.csv("/path/to/file/metadata_example.csv")
met <- tibble::column_to_rownames(met, var = "X")

#subset the data and make the metadata groups by which to facet_wrap
#change this to however you want to subset for your example data
samples <- met %>%
  filter(str_detect(SampleID, regex("\\bSample", ignore_case = TRUE)))

#extract sample ID's for the metadata group
subset_samples <- samples$SampleID

#subset the merged abundance output by your metadata group and save
#df_list is the output of merged_abundance_processing
#colname_list is the subset of samples in your column names found in the output of merged_abundance_processing
#set the output_dir for where you want the files saved
#set a prefix for the new output files
#this prefix will also be added to each output subset file
#the outputs are TaxaLevel_ReadCount and TaxaLevel_RelAbun for the read count and relative abudance files, respectively (TaxaLevels are Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain)
#the outputs are saved both to .csv files and to the global environment
subset_by_colnames_and_save <- function(df_list, colname_list, prefix, output_dir) {
  subset_list <- list()
  
  for (df_name in names(df_list)) {
    df <- df_list[[df_name]]
    
    # Subset the dataframe by the given list of column names
    subset_df <- df[, colnames(df) %in% colname_list, drop = FALSE]
    
    # Create the name for the new dataframe
    new_df_name <- paste0(prefix, "_", df_name)
    
    # Assign the subset dataframe to the global environment with the new name
    assign(new_df_name, subset_df, envir = .GlobalEnv)
    
    # Add the subset dataframe to the list with the new name
    subset_list[[new_df_name]] <- subset_df
    
    # Save the dataframe to a CSV file in the specified directory
    write.csv(subset_df, file.path(output_dir, paste0(new_df_name, ".csv")), row.names = TRUE)
  }
  
  return(subset_list)
}

new_dir <- ("/path/to/file/output/subset")

subset_relabun_list <- subset_by_colnames_and_save(relabun_list, subset_samples, prefix = "Subset", output_dir = new_dir)

#subset the most abundance taxa and subset the rest to "Other"
Genus_otu_most = mostAbundant(Subset_Genus_RelAbun, N = 15, items = NULL, others = TRUE, rescale = TRUE)
Species_otu_most = mostAbundant(Subset_Species_RelAbun, N = 15, items = NULL, others = TRUE, rescale = TRUE)

#set metagroups for the facet_wrap for plots 
#metagroups for plots 
sample_metagroup <- list('All_Samples' = subset_samples)

#if you have multiple groups to define, define the samples in each group, then the metagroup
five_samples <- subset(samples, Timepoint == "5_month")
five_samples <- five_samples$SampleID
two_samples <- subset(samples, Timepoint == "2_week")
two_samples <- two_samples$SampleID
time_samples <- (c(two_samples, five_samples))

timepoint_metagroup <- list('5_month' = five_samples, '2_week' = two_samples)

#if your metagroups are a subset, make sure to subset the relative abundance as well
subset_genus <- Genus_otu_most[,subset_samples]
subset_species <- Species_otu_most[,subset_samples]

timepoint_genus <- Genus_otu_most[,time_samples]
timepoint_species <- Species_otu_most[,time_samples]

#set colors for the taxa in the plot
color_taxa <- function(data){
  taxa <- as.data.frame(data)
  taxa_to_color <- as.factor(rownames(data))
  names <- levels(taxa_to_color)
  colors = c('#bdbcbc', '#85A1EF', '#b2d689', '#FFEA00', '#f7766d', '#196619', '#aaf0d1','#0047AB', '#f5bcbc', '#AB6666',
             '#f59a40', '#c8b1d5','#D22B2B', '#6b3e96','#FF7518', '#FF69B4')
  colorme <- setNames(colors, unique(as.character(names)))
  return(colorme)
}

color_taxa_gen <- color_taxa(Genus_otu_most)
color_taxa_spec <- color_taxa(Species_otu_most)

#set the labels for the taxa which will be used in the plot
label_fill_gen <- rownames(Genus_otu_most)
label_fill_spec <- rownames(Species_otu_most)

#taxa_plot
taxa_plot = function(data, label_x = 'Samples', label_y = 'Abundances', label_fill = 'Features', color = NULL, base_size = 11, max_scale_value = NULL, metadata_groups = NULL, nrow = 2, ncol = 2)
{
  sample = item = abun = NULL # to appease R CMD check (they are later used by ggplot2's aes in the context of data_melt, but that syntax bother's R CMD check)
  if (!is.data.frame(data) & !is.matrix(data)) { stop('The first argument must be a matrix or a data frame') }
  if(!is.null(max_scale_value) & !is.numeric(max_scale_value)) { stop('max_scale_value must be numeric') }
  data=t(data)
  data_melt = reshape2::melt(as.matrix(data), value.name = 'abun')
  colnames(data_melt) = c('sample', 'item', 'abun')
  data_melt$sample = as.factor(data_melt$sample)
  data_melt$item = as.factor(data_melt$item)
  data_melt$abun = as.numeric(data_melt$abun)
  
  if (!is.null(metadata_groups)) 
  {
    if(sum(sapply(metadata_groups, length)) != length(unique(data_melt$sample)))
    {stop('metadata_groups must contain the same samples that data')}	    
    data_melt$group = apply(data_melt, 1, function(s) unlist(lapply(names(metadata_groups), function(x) if( s['sample'] %in% metadata_groups[[x]]){x})))
  }
  
  #PLOT DATA
  p = ggplot2::ggplot(data_melt, ggplot2::aes(x = sample, y = abun, fill = item))
  p = p + ggplot2::geom_col()
  # There are two types of bar charts: geom_bar() and geom_col().
  # geom_bar() makes the height of the bar proportional to the number of cases in each group (or if the weight aesthetic is supplied, the sum of the weights).
  # If you want the heights of the bars to represent values in the data, use geom_col() instead.
  # geom_bar() uses stat_count() by default: it counts the number of cases at each x position.
  # geom_col() uses stat_identity(): it leaves the data as is.
  # geom_col == geom_bar(stat = 'identity')
  
  if(is.null(label_x)) { label_x = '' }
  if(is.null(label_y)) { label_y = '' }
  if(is.na(label_x)  ) { label_x = '' }
  if(is.na(label_y)  ) { label_y = '' }
  
  if(nchar(label_x)==0) { theme_x = ggplot2::element_blank()
  }else{ theme_x = ggplot2::element_text() }
  if(nchar(label_y)==0) { theme_y = ggplot2::element_blank()
  }else{ theme_y = ggplot2::element_text() }
  
  p = p + ggplot2::theme_light(base_size = base_size)
  p = p + ggplot2::theme(axis.title.x = theme_x, axis.title.y = theme_y)
  p = p + ggplot2::labs(x = label_x, y = label_y, fill = label_fill)
  if (!is.null(color) & length(color) >= length(unique(as.character(data_melt$item))))
  {
    p = p + ggplot2::scale_fill_manual( values = setNames(color, unique(as.character(data_melt$item))) )
  }
  if (!is.null(max_scale_value))
  {
    p = p + ggplot2::ylim(0, max_scale_value)
  }
  if (!is.null(metadata_groups))
  {
    p = p + ggplot2::facet_wrap(.~group, scales = 'free', nrow = nrow, ncol = ncol) + ggplot2::theme(strip.text = ggplot2::element_text(colour = 'black'))
  }
  
  if (length(unique(data_melt$sample)) > 1)
  {
    p = p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }	
  return(p)
}

#make the plot and save it
#set the color with your colors we made above
#set the labels with the labels from above
#base size is the font size
#label_y is the y_lab
#metadata_groups are the groups we set above which will result in facet_wrap
#your nrow and ncol are the number of columns and rows in the facet_wrap - it must equal or exceed the number of groups you have made
pdf("/path/to/file/status.genus.pdf", width = 8, height = 5) #change the height and width of the pdf as you desire
taxa_plot(data = timepoint_genus, 
          label_y = "Relative Abundance", 
          color = color_taxa_gen, 
          label_fill = label_fill_gen , 
          base_size = 12, 
          nrow = 1,
          ncol = 2,
          metadata_groups = timepoint_metagroup)
dev.off()

#save your data used in the plot as well 
write.csv(timepoint_genus, "/path/to/file/timepoint.genus.csv")



#to plot more than 15 species and set their colors and their order in the plot
#take out most abundant taxa and condense the rest to other
#first I subset by N = 30 which will return the top 30 most abundant species
Species_otu_most = mostAbundant(Subset_Species_RelAbun, N = 30, items = NULL, others = TRUE, rescale = TRUE)

#I then make the rownames into a column in order to reorder it
Species_otu_most <- Species_otu_most %>% rownames_to_column(var = "Species")

#I set a custom order for the taxa in the order I'd like them to appear on the plot
custom_order <- c("Other", 
                  "Bifidobacterium bifidum", 
                  "Bifidobacterium breve", 
                  "Bifidobacterium dentium", 
                  "Bifidobacterium longum", 
                  "Bifidobacterium pseudocatenulatum",
                  "Klebsiella michiganensis", 
                  "Klebsiella oxytoca", 
                  "Klebsiella pneumoniae", 
                  "Klebsiella variicola",
                  "Escherichia coli",
                  "Veillonella parvula", 
                  "Veillonella ratti", 
                  "Veillonella rogosae",
                  "Bacteroides fragilis", 
                  "Bacteroides ovatus", 
                  "Citrobacter freundii", 
                  "Citrobacter sp RHBSTW 00671", 
                  "Phocaeicola dorei", 
                  "Phocaeicola vulgatus", 
                  "Streptococcus salivarius",
                  "Parabacteroides distasonis", 
                  "Parabacteroides merdae", 
                  "Staphylococcus epidermidis", 
                  "Erysipelatoclostridium ramosum", 
                  "Enterococcus faecalis", 
                  "Clostridium butyricum",
                  "Hungatella hathewayi", 
                  "Megasphaera sp MJR8396C", 
                  "Ruminococcus gnavus", 
                  "Ruminococcus torques")

#I reorder the relative abundance table by my custom order
Species_otu_most <- Species_otu_most[match(custom_order, Species_otu_most$Species), ]

#then I reset the rownames as the taxa
rownames(Species_otu_most) = NULL
Species_otu_most <- Species_otu_most %>% column_to_rownames(var = "Species")

#I set custom colors for the species
color_taxa_spec <- c("Other" = "#bdbcbc", 
                  "Bifidobacterium bifidum" = "#85A1EF", 
                  "Bifidobacterium breve" = "#B7D8E8", 
                  "Bifidobacterium dentium" = "#7EC8E3", 
                  "Bifidobacterium longum" = "#0F52BA", 
                  "Bifidobacterium pseudocatenulatum" = "#00FFFF",
                  "Klebsiella michiganensis"  = "#b2d689", 
                  "Klebsiella oxytoca" = "#B4C424", 
                  "Klebsiella pneumoniae" = "#7CFC00", 
                  "Klebsiella variicola" = "#DFFF00",
                  "Escherichia coli" = '#FFEA00',
                  "Veillonella parvula" = "#f7766d", 
                  "Veillonella ratti" = "#b54e4b", 
                  "Veillonella rogosae"= "#732728",
                  "Bacteroides fragilis" = "#196619", 
                  "Bacteroides ovatus" = "#123D12", 
                  "Citrobacter freundii" = "#D4FBE7", 
                  "Citrobacter sp RHBSTW 00671"= "#75D7A2", 
                  "Phocaeicola dorei" = "#0000FF", 
                  "Phocaeicola vulgatus"  = "#00008B", 
                  "Streptococcus salivarius" = "#f5bcbc" ,
                  "Parabacteroides distasonis" = "#AB6666", 
                  "Parabacteroides merdae" = "#F5DEB3" , 
                  "Staphylococcus epidermidis"= "#f59a40" , 
                  "Erysipelatoclostridium ramosum"= "#c8b1d5", 
                  "Enterococcus faecalis" = "#D22B2B", 
                  "Clostridium butyricum" = "#FF7518",
                  "Hungatella hathewayi"= "#FF10F0", 
                  "Megasphaera sp MJR8396C"= "#9F2B68", 
                  "Ruminococcus gnavus" ="#5D3FD3", 
                  "Ruminococcus torques" = "#800080")


#set the labels
label_fill_spec <- rownames(Species_otu_most)

#subset the relative abundance table as needed
timepoint_species <- Species_otu_most[,time_samples]

#make the plot
pdf("/path/to/file/timepoint.species_extended.pdf", width = 12, height = 5) #change the height and width as needed
taxa_plot(timepoint_species, 
          label_y = "Relative Abundance", 
          color = color_taxa_spec, 
          label_fill = label_fill_spec , 
          base_size = 12, 
          nrow = 1,
          ncol = 2,
          metadata_groups = timepoint_metagroup)
dev.off()

#save the data used to make the plot
write.csv(timepoint_species, "/path/to/file/timepoint.species_extended.csv")

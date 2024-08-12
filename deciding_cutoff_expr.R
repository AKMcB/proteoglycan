####################
###load libraries###
####################
library(tidyverse)
library(data.table)
library(ggpubr)
library(mclust)
library(optparse)

########################
# Define option parser #
########################

option_list <- list(
  make_option(c("-e", "--expression"), type="character", default=NULL,
              help="expression data. Log2+1 values. Obligatory"), 
  
  make_option(c("-c", "--clinical"), type="character", default=NULL,
              help="clinical data. Obligatory"),
  
  make_option(c("-l", "--list"), type="character", default=NULL,
              help="list of genes common between ccle and tcga with high variance. Obligatory"),
  
  make_option(c("-t", "--type"), type="character", default=NULL,
              help="Cancer type of interest.Obligatory"),
  
  make_option(c("-g", "--gene"), type="character", default=NULL,
              help="gene of interest. Obligatory"),
  
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file path for tile plot for all genes"))


# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


EXPR <- "rawdata/expression_data_tcga_log2.csv"

CLIN <- "rawdata/clinical_tcga_data.csv"

LIST <- "figures/tileplot/filtered_variance_genes_cancer_type_tcga.csv"

GENE <- "SPOCK1"

CANCER <- "BRCA"

OUT <- "figures/tileplot/"
### load tcga data 
tcga <- as.data.frame(fread(EXPR, header = TRUE))
colnames(tcga)[1] <- "gene"

### load clinical tcga data 
clin<- as.data.frame(fread(CLIN))
clin <- clin[,-1]

### load matched gene data 
gene_list <- as.data.frame(fread(LIST, header = TRUE))
gene_list <- gene_list[,-1]

############################
#### filter based on trim ##
############################
# filter tcga 
gene_tcga <- subset(tcga, tcga$gene %in% gene_list$name)

# to get id on column 1 and trims on rows
rownames(gene_tcga) <- gene_tcga[,1]
gene_tcga$gene <- NULL
gene_tcga <- as.data.frame(t(gene_tcga))
gene_tcga <- rownames_to_column(gene_tcga, "id")

merged <- merge(gene_tcga, clin, by.x = "id", by.y = "sample")
merged <- select(merged, c(id, gene_list$name, `tcga code`))

#Transform to long format
merged_long <- pivot_longer(merged, cols = gene_list$name)
saveRDS(merged_long, file = "expression_values.rds")

rownames(merged_long) <- merged_long[,1]
merged_long$id <- NULL

variance_data <- merged_long %>%
  group_by(`tcga code`, name) %>%
  summarize(VarianceExpression = var(value, na.rm = TRUE)) %>%
  ungroup()

merged_with_variance <- merged_long %>%
  left_join(variance_data, by = c("tcga code", "name"))

filt <- merged_with_variance$VarianceExpression > 0.745 #cutoff value based on housekeeping genes
merged.filt <- merged_with_variance[filt,]

################################################
## Create histogram and caluclate expr cutoff ##
################################################
# Create an empty list to store cutoff values
cutoff_values <- list()

# Function to calculate cutoff value for a gene and cancer type
calculate_cutoff <- function(merged.filt, GENE, CANCER) {
  data <- merged.filt %>% filter(name == GENE & `tcga code` == CANCER) %>%
    pull(value)
  
  fit <- Mclust(data, G = 2)
  means <- fit$parameters$mean
  v <- mean(means)
  
  # Store the cutoff value for this gene and cancer type
  cutoff_values[[paste(GENE, CANCER, sep = "_")]] <<- v
  
  hist_GMM <- ggplot(merged.filt %>% filter(name == GENE & `tcga code` == CANCER), aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), fill = "blue", color = "black", alpha = 0.5, bins = 30) +
    geom_density(color = "red", size = 1) +
    geom_vline(xintercept = v, color = "red", linetype = "dashed", linewidth = 1.5) +  # Add vertical line
    labs(y = "Frequency", x = "Expression (log2+1)") +
    theme(legend.position = "none",
          plot.title = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))  
  
  return(hist_GMM)
}

# Loop through each gene and cancer type
for (GENE in unique(merged.filt$name)) {
  for (CANCER in unique(merged.filt$`tcga code`)) {
    # Check if the gene is present in this cancer type
    if (any(merged.filt$name == GENE & merged.filt$`tcga code` == CANCER)) {
      # Calculate cutoff value and store in list
      hist_plot <- calculate_cutoff(merged.filt, GENE, CANCER)
      plot_output_filename <- paste(GENE, CANCER, ".png", sep = "_")
      ggsave(plot_output_filename, hist_plot, width = 8, height = 6)
    }
  }
}

# Convert cutoff values list to a data frame
cutoff_values_df <- do.call(rbind, lapply(names(cutoff_values), function(x) {
  data.frame(GENE = gsub(".*_(.*)", "\\1", x), CANCER = gsub("_(.*$)", "", x), cutoff_value = cutoff_values[[x]], stringsAsFactors = FALSE)
}))



colnames(cutoff_values_df) <- c("cancer", "gene", "cutoff_value")
output_file_table1 <- paste(OUT,"expr_cutoff_genes_cancer_type_tcga", ".csv", sep="")
fwrite(cutoff_values_df, output_file_table1, row.names = TRUE)








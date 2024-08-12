#############
# Libraries #
#############

library(tidyverse)
library(data.table)
library(ggpubr)
library(optparse)
library(mclust)
library(edgeR)

########################
# Define option parser #
########################

option_list <- list(
  make_option(c("-e", "--expression"), type="character", default=NULL,
              help="Preprocessed tcga expression file.
              Genes should have columname as gene.Obligatory"),
  
  make_option(c("-c", "--clinical"), type="character", default=NULL,
              help="Preprocessed clinical data. Obligatory"),
  
  make_option(c("-g", "--gene"), type="character", default=NULL,
              help="gene list of interest. Just one column of gene names, 
              the full name is optional. Obligatory"), 
  
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file path for tile plot for all genes"))


# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


EXPR <- "rawdata/expression_data_tcga_log2.csv"

CLIN <- "rawdata/clinical_tcga_data.csv"

GENE <- "rawdata/proteoglycan_erik_list.csv"

OUT <- "figures/tileplot/"

########################
# Read expression data #
########################

expr <- as.data.frame(fread(EXPR, header=TRUE))
colnames(expr)[1] <- "gene"


######################
# Read clinical data #
######################

clin <- as.data.frame(fread(CLIN, header = TRUE))
clin <- clin[,-1]

##################
# Read gene file # 
##################

gene <- as.data.frame(fread(GENE))


#################################
# subset based on clinical data # 
#################################
#Filter expression file based on the gene list
gene_expr <- merge(expr, gene, by = "gene")

gene_expr <- column_to_rownames(gene_expr, "gene")

#column to row and row to column 

gene_expr <- as.data.frame(t(gene_expr))
gene_expr <- rownames_to_column(gene_expr, "id")

merged <- merge(gene_expr, clin, by.x = "id", by.y = "sample")

sort(unique(merged$`cancer type abbreviation`))


#####################
# Clean merged file #
#####################
table(merged$`tcga code`) #to check if all codes are included --> 33 codes, all codes are included 

#####Test for duplicates 
dup <- as.data.frame(duplicated(merged$id))

dup_rows <- merged[dup$`duplicated(merged$id)`, ]
#no duplicates

#########################
# Calculate mean/median #
#########################

#Transform to long format
merged_long <- pivot_longer(merged, cols = gene$gene) 
 

variance_data <- merged_long %>%
  group_by(`tcga code`,name) %>%
  summarize(VarianceExpression = var(value)) 

filt <- variance_data$VarianceExpression > 0.745 #cutoff value based on housekeeping genes
variance.filt <- variance_data[filt,]


# Create a complete grid of all possible combinations
complete_df <- expand.grid(
  `tcga code` = unique(variance.filt$`tcga code`),
  name = unique(variance.filt$name)
)

# Merge the complete grid with the actual data to identify missing combinations
complete_df <- complete_df %>%
  left_join(variance.filt, by = c("tcga code", "name"))

# Count the number of missing values per cancer type
cancer_missing <- complete_df %>%
  group_by(`tcga code`) %>%
  summarise(missing_count = sum(is.na(VarianceExpression))) %>%
  arrange(desc(missing_count))

# Count the number of missing values per gene
gene_missing <- complete_df %>%
  group_by(name) %>%
  summarise(missing_count = sum(is.na(VarianceExpression))) %>%
  arrange(desc(missing_count))

# Convert cancer_type and gene columns to factors based on the missing values count
complete_df$`tcga code` <- factor(complete_df$`tcga code`, levels = cancer_missing$`tcga code`)
complete_df$name <- factor(complete_df$name, levels = rev(gene_missing$name))

complete_df <- na.omit(complete_df)

tileplot <- ggplot(complete_df, aes(y = fct_rev(name), x = fct_rev(`tcga code`), fill = VarianceExpression)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),  
        strip.background = element_blank())

output_file_tile1 <- paste(OUT,"tileplot_filtered_variance_genes_tcga", ".pdf", sep="")

pdf(output_file_tile1, width = 15, height = 20, onefile = F)
print(tileplot)
dev.off()


output_file_table1 <- paste(OUT,"filtered_variance_genes_tcga", ".csv", sep="")
fwrite(complete_df, output_file_table1, row.names = TRUE)

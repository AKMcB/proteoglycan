#############
# Libraries #
#############

library(tidyverse)
library(data.table)
library(ggpubr)
library(optparse)

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

GENE <- "rawdata/housekeeping_genes.csv"

OUT <- "figures/"

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


#Filter expression file based on the gene list
gene_expr <- merge(expr, gene, by = "gene")

gene_expr <- column_to_rownames(gene_expr, "gene")
gene_expr <- as.data.frame(t(gene_expr))
gene_expr <- rownames_to_column(gene_expr, "id")
merged <- merge(gene_expr, clin, by.x = "id", by.y = "sample")

#Transform to long format
merged_long <- pivot_longer(merged, cols = gene$gene) 


#Calculate the median expression for each trim protein in each of the TCGA cancer code 
variance_data <- merged_long %>%
  group_by(`tcga code`,name) %>%
  summarize(VarianceExpression = var(value)) 


tileplot <- ggplot(variance_data, aes(x = `tcga code`, y = name, fill =VarianceExpression)) +
  geom_tile() + 
  scale_fill_viridis_c() +
  theme_bw() + 
  labs(y = " ", x = "Cancer Type", fill = "Variance") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),  
        strip.background = element_blank())

tileplot

pdf("tileplot_housekeeping_genes.pdf", width = 15, height = 20, onefile = F)
print(tileplot_median)
dev.off()

####################################
# Decide cut-off value of variance #
####################################

hist(variance_data$VarianceExpression)
lines(density(variance_data$VarianceExpression))
shapiro.test(variance_data$VarianceExpression)
#data:  variance_data$VarianceExpression
#W = 0.93861, p-value = 1.897e-10 #not normally distributed

c <- mean(variance_data$VarianceExpression)
c_1 <- median(variance_data$VarianceExpression)
c
c_final <- c+c














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
  make_option(c("-c", "--cox"), type="character", default=NULL,
              help="results from the cos univariate analysis.obligatory"),
  
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file path for tile plot for all genes"))


# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


RESULTS <- "figures/survival/cox_univariate_pfi_results.csv"

OUT <- "figures/survival/"

########################
# Read expression data #
########################

# Load data
res <- as.data.frame(fread(RESULTS, header=TRUE))
res <- res[,-1]
colnames(res)[7] <- "Pr(>|z|)"

# Filter data
res.filt <- subset(res, res$`Pr(>|z|)` <= 0.05)
res.filt$exp.coef. <- log2(res.filt$exp.coef.)
res.filt$exp.coef. <- abs(res.filt$exp.coef.)
res.filt$coef_sign <-  ifelse(res.filt$coef > 0, '>0', "<0")
res.filt$coef_sign <- factor(res.filt$coef_sign, levels = c('>0', "<0"))


# Create a complete grid of all possible combinations
complete_df <- expand.grid(
  Cancer = unique(res.filt$Cancer),
  Gene = unique(res.filt$Gene)
)

# Merge the complete grid with the actual data to identify missing combinations
complete_df <- complete_df %>%
  left_join(res.filt, by = c("Cancer", "Gene"))

# Count the number of non-missing values per cancer type
cancer_count <- complete_df %>%
  group_by(Cancer) %>%
  summarise(value_count = sum(!is.na(exp.coef.))) %>%
  arrange(desc(value_count))

# Count the number of non-missing values per gene
gene_count <- complete_df %>%
  group_by(Gene) %>%
  summarise(value_count = sum(!is.na(exp.coef.))) %>%
  arrange(desc(value_count))

# Convert cancer_type and gene columns to factors based on the non-missing values count
res.filt$Cancer <- factor(res.filt$Cancer, levels = cancer_count$Cancer)
res.filt$Gene <- factor(res.filt$Gene, levels = rev(gene_count$Gene))

# Define size breaks and labels
size_breaks <- c(0.5,1,2,4)
size_labels <- as.character(size_breaks)

# Plot
p <- ggplot(res.filt, aes(y = Gene, x = Cancer)) +
  geom_point(aes(size = exp.coef., colour = coef_sign)) +
  theme_bw() +
  scale_size_continuous(
    range = c(0.5, 4),  # Adjust the range as needed
    breaks = size_breaks,
    labels = size_labels
  ) +
  scale_color_manual(values = c(">0" = "red", "<0" = "blue")) +
  labs(
    x = "Cancer Types", 
    y = "Genes", 
    size = "Hazard Ratio",  # Title for the size legend
    colour = "Coefficient Sign",  # Title for the color legend
    fill = "Coefficient Value"
  ) +
  ggtitle("Coefficients of Survival:PFI") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_blank(),  
    strip.background = element_blank()
  )

p

table(res.filt$coef_sign)

pdf("dotplot_survival_pfi.pdf", width = 7, height = 7, onefile = F)
print(p)
dev.off()


fwrite(res.filt, "cofficient_sign_pfi.csv", row.names = T)







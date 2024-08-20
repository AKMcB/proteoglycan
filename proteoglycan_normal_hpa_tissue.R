library(tidyverse)
library(data.table)
library(ggthemes)
library(ComplexHeatmap)
library(dendsort)
library(circlize)
library(RColorBrewer)
###############
## Read data ##
###############

data <- fread("Normal tissue expression/rna_single_cell_type_tissue.tsv")
unique(data$Tissue)

#Subset the data based on the gene of interest 
genes <- as.data.frame(fread("rawdata/proteoglycan_erik_list.csv"))


data.subset <- subset(data, data$`Gene name` %in% genes$gene)
data.subset <- subset(data.subset, data.subset$Tissue == "kidney")
unique(data.subset$`Cell type`)

attach(data.subset)

data.subset$log2 <- log2(data.subset$nTPM +1)
min(data.subset$log2) #should not be any minus values 
max(data.subset$log2) #Should not be higher than 18 


mean_data <- data.subset %>%
  group_by(`Cell type`,`Gene name`) %>%
  summarize(MeanExpression = mean(log2)) 

test <- pivot_wider(mean_data, 
                    names_from = `Cell type`, 
                    values_from = MeanExpression)

test <- column_to_rownames(test, "Gene name")
#test <- as.data.frame(t(test))
hs <- test
# Compute row variances
row_variances <- apply(hs, 1, var)

# Identify rows with zero variance
zero_variance_rows <- which(row_variances == 0)

# Display zero variance rows
print(zero_variance_rows)
rows_to_remove <- c("ACAN" ,   "ASPN",    "CHAD", "COL15A1" ,   "EPYC"  ,  "KERA" ,   "NCAN"  ,  "OGN"  , "OPTC" , "PTPRZ1" , "SPOCK3" )

#hs <- hs[-grep('PSMA8', row.names(hs)),, drop = FALSE]
hs <- hs[!(row.names(hs) %in% rows_to_remove),]
# Compute column variances
column_variances <- apply(hs, 2, var)

# Identify columns with zero variance
zero_variance_columns <- which(column_variances == 0)

# Display zero variance columns
print(zero_variance_columns)


h <- as.matrix(hs)
#h <- hs_scaled <-t(scale(t(test)))

h_breaks <- seq(min(h), max(h), length.out = 10)
h_breaks

#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

h_breaks <- quantile_breaks(h, n = 11)
h_breaks

row_dist <- as.dist(1-cor(t(h), method = "pearson"))
row_hc <- hclust(row_dist, method = "complete")

#We can use the dendsort package and reorder the clustering
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
#c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000")


col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0467d2", "#0356af","#023369", "#130202", "#6b0002","#b20004","#ff0000"))

?colorRamp2

pdf("proteoglycan_hpa_cell_types_kidney_heatmap.pdf",  width=13,height=16)

ht <- Heatmap(h,col = col_fun,
              cluster_columns = F,
              name = "Expression Values",
              show_heatmap_legend = T ,
              #top_annotation = ha,
              #left_annotation = hr,
              show_column_names = T,
              show_row_names = T,
              cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=12),
              column_title_gp = gpar(fontsize=12),
              height = unit(15, "cm"),
              width = unit(22, "cm"),
              column_dend_reorder = F,
              row_dend_reorder = T,
              show_row_dend = T,
              border = T,
              column_dend_height = unit(2, "cm"),
              column_names_rot = 90,
              #legend_height = unit(4, "cm"),
              column_names_gp = grid::gpar(fontsize = 10),
              row_split = 2,
              row_gap = unit(c(0.7,0.7), "mm"),
              heatmap_legend_param = list(title="Expression (log2+1)", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=10)))


draw(ht,column_title = "Proteoglycan in Normal Kidney Tissue Cell Type", column_title_gp = gpar(fontsize = 16), 
     merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))

dev.off()
###################################################################################

tileplot <- ggplot(data.subset, aes(y = (`Gene name`), x = (`Cell type`), fill = log2)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "white") +
  theme_classic() +
  labs(x = "Cancer Type", y = "Proteoglycan", fill = "Variance Expression", title = "Variance Expression of Proteoglycan in Normal Tissue Cell Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,, vjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),  
        strip.background = element_blank())


tileplot

############################
## Calculate the variance ##
############################
attach(data.subset)

variance_data <- data.subset %>%
  group_by(`Cell type`,`Gene name`) %>%
  summarize(VarianceExpression = var(log2)) 

mean_data <- data.subset %>%
  group_by(`Cell type`,`Gene name`) %>%
  summarize(MeanExpression = mean(log2)) 



filt <- variance_data$VarianceExpression > 0.745 #cutoff value based on housekeeping genes
variance.filt <- variance_data[filt,]

variance.filt <- mean_data

# Create a complete grid of all possible combinations
complete_df <- expand.grid(
  `Cell type` = unique(variance.filt$`Cell type`),
  `Gene name` = unique(variance.filt$`Gene name`)
)

# Merge the complete grid with the actual data to identify missing combinations
complete_df <- complete_df %>%
  left_join(variance.filt, by = c("Cell type", "Gene name"))

# Count the number of missing values per cancer type
cancer_missing <- complete_df %>%
  group_by(`Cell type`) %>%
  summarise(missing_count = sum(is.na(MeanExpression))) %>%
  arrange(desc(missing_count))

# Count the number of missing values per gene
gene_missing <- complete_df %>%
  group_by(`Gene name`) %>%
  summarise(missing_count = sum(is.na(MeanExpression))) %>%
  arrange(desc(missing_count))

# Convert cancer_type and gene columns to factors based on the missing values count
complete_df$`Cell type` <- factor(complete_df$`Cell type`, levels = cancer_missing$`Cell type`)
complete_df$`Gene name` <- factor(complete_df$`Gene name`, levels = rev(gene_missing$`Gene name`))

complete_df_1 <- na.omit(complete_df)

tileplot <- ggplot(complete_df_1, aes(y = fct_rev(`Gene name`), x = fct_rev(`Cell type`), fill = MeanExpression)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "white") +
  theme_classic() +
  labs(x = "Cancer Type", y = "Proteoglycan", fill = "Variance Expression", title = "Variance Expression of Proteoglycan in Normal Tissue Cell Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,, vjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),  
        strip.background = element_blank())
tileplot

pdf("tileplot_proteoglycans_celltype_HPA.pdf", width = 15, height = 20, onefile = F)
print(tileplot)
dev.off()


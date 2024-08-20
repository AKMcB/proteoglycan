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

data <- fread("Normal tissue expression/rna_single_cell_type_germ_layer.csv")
data <- data[,-1]
length(unique(data$`Cell type`))

#Subset the data based on the gene of interest 
genes <- as.data.frame(fread("rawdata/proteoglycan_erik_list.csv"))

#Add gene of interes 
add_genes <- as.data.frame(c("SNAI1", "SNAI2", "ZEB1", "ZEB2","TWIST1", "TWIST2"))
colnames(add_genes) <- "gene"

genes_1 <- rbind(genes, add_genes)

data.subset <- subset(data, data$`Gene name` %in% genes_1$gene)
#data.subset <- subset(data.subset, data.subset$Tissue == "kidney")
#length(unique(data.subset$`Cell type`))

attach(data.subset)

####################
## log2 transform ##
####################
data.subset <- data.subset %>%
  group_by(`Cell type`,`Gene name`) %>%
  summarize(MeanExpression = mean(nTPM)) 


data.subset$log2 <- log2(data.subset$MeanExpression +1)
min(data.subset$log2) #should not be any minus values 
max(data.subset$log2) #Should not be higher than 18

data.subset <- data.subset %>% select(c("Gene name", "Cell type", "log2"))

data.subset_1 <- pivot_wider(data.subset, 
                    names_from = `Cell type`, 
                    values_from = log2)

test <- data.subset_1
#make annotation of the genes of interest 
ann1 <- subset(test, test$`Gene name` %in% add_genes$gene)
ann1 <- distinct(ann1, ann1$`Gene name`,.keep_all = T)
ann1$`ann1$\`Gene name\`` <- NULL
ann1$origin <- NULL
ann1 <- column_to_rownames(ann1, "Gene name")
ann1 <- as.data.frame(t(ann1))
ann1 <- rownames_to_column(ann1, "type")

ann2 <- data
ann2 <- data %>% select(c("Cell type", origin))
ann2 <- distinct(ann2, ann2$`Cell type`, .keep_all = T)
ann2$`ann2$\`Cell type\`` <- NULL
ann2 <- merge(ann2, ann1, by.x = "Cell type", by.y = "type")
ann2 <- ann2[,-c(3:8)]
ann2 <- column_to_rownames(ann2, "Cell type")

ann1 <- column_to_rownames(ann1, "type")
#remove the gene of interest NULL#remove the gene of interest in the following analysis 
test <- subset(test, test$`Gene name` %in% genes$gene)
test$origin <- NULL
test <- distinct(test, test$`Gene name`, .keep_all = T)
test$`test$\`Gene name\`` <- NULL
test <- column_to_rownames(test, "Gene name")

all(rownames(ann1) == colnames(test))
all(rownames(ann2) == colnames(test))
#####################################################


h <- as.matrix(test)
#h <- hs_scaled <-t(scale(t(test)))

# Compute row variances
row_variances <- apply(h, 1, var)

# Identify rows with zero variance
zero_variance_rows <- which(row_variances == 0)

# Display zero variance rows
print(zero_variance_rows)
#rows_to_remove <- c("ACAN","ASPN","CHAD","COL15A1","EPYC","KERA","NCAN","OPTC","PTPRZ1","SPOCK3")

#hs <- hs[-grep('PSMA8', row.names(hs)),, drop = FALSE]
#h <- h[!(row.names(h) %in% rows_to_remove),]
# Compute column variances
column_variances <- apply(h, 2, var)

# Identify columns with zero variance
zero_variance_columns <- which(column_variances == 0)

# Display zero variance columns
print(zero_variance_columns)


h_breaks <- seq(min(h), max(h), length.out = 10)
h_breaks

ann1$SNAI1 <- as.numeric(ann1$SNAI1)
a <- ann1$SNAI1 
a_breaks <- seq(min(a), max(a), length.out = 10)
a_breaks

ann1$SNAI2 <- as.numeric(ann1$SNAI2)
b <- ann1$SNAI2 
b_breaks <- seq(min(b), max(b), length.out = 10)
b_breaks

ann1$ZEB1 <- as.numeric(ann1$ZEB1)
d <- ann1$ZEB1 
d_breaks <- seq(min(d), max(d), length.out = 10)
d_breaks

ann1$ZEB2 <- as.numeric(ann1$ZEB2)
e <- ann1$ZEB2
e_breaks <- seq(min(e), max(e), length.out = 10)
e_breaks

ann1$TWIST1 <- as.numeric(ann1$TWIST1)
f <- ann1$TWIST1 
f_breaks <- seq(min(f), max(f), length.out = 10)
f_breaks

ann1$TWIST2 <- as.numeric(ann1$TWIST2)
g <- ann1$TWIST2
g_breaks <- seq(min(g), max(g), length.out = 10)
g_breaks

#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

h_breaks <- quantile_breaks(h, n = 11)
h_breaks

a_breaks <- quantile_breaks(a, n = 11)
a_breaks

b_breaks <- quantile_breaks(b, n = 11)
b_breaks

d_breaks <- quantile_breaks(d, n = 11)
d_breaks

e_breaks <- quantile_breaks(e, n = 11)
e_breaks

f_breaks <- quantile_breaks(f, n = 11)
f_breaks

g_breaks <- quantile_breaks(g, n = 11)
g_breaks
#Make clusters based on their pearson correlation
row_dist <- as.dist(1-cor(t(h), method = "pearson"))
col_dist <- as.dist(1-cor(h, method = "pearson"))
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")


#We can use the dendsort package and reorder the clustering
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")

 
#c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000")


#col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5","#023369", "#130202", "#6b0002","#b20004","#ff0000"))

#pdf("proteoglycan_hpa_cell_types_kidney_heatmap.pdf",  width=13,height=16)

col_fun = colorRamp2(c(0,15), c("blue", "#ff0000"))
col_fun_2 = colorRamp2(c(0,15), c("blue", "#ff0000"))
colors <- c("red", "green", "orange", "yellow", "grey")
ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))

ha <- HeatmapAnnotation("SNAI1"= ann1$SNAI1,
                        "SNAI2" = ann1$SNAI2,
                        "ZEB1" = ann1$ZEB1,
                        "ZEB2" = ann1$ZEB2,
                        "TWIST1" = ann1$TWIST1,
                        "TWIST2" = ann1$TWIST2,
                        "Germ Layer" = ann2$origin,
                        col = list("SNAI1"= col_fun_2,
                                   "SNAI2" = col_fun_2,
                                   "ZEB1" = col_fun_2,
                                   "ZEB2" = col_fun_2,
                                   "TWIST1" = col_fun_2,
                                   "TWIST2" = col_fun_2,
                                   "Germ Layer" = c(
                                     "Ectoderm" = "#3ac4f8",
                                     "Endoderm" = "#99EDCC",
                                     "Mesoderm" = "#9A275A",
                                     "Trophoblast" = "#CB958E", 
                                     "Varies" = "#E36588"
                                   )),
                        simple_anno_size = unit(0.30, "cm"),
                        annotation_name_side = "right",
                        border = T,
                        annotation_name_gp= gpar(fontsize = 12,fontface="bold"),
                        show_legend = T,
                        annotation_legend_param = list(title_gp = gpar(fontsize=12, fontface="bold"), labels_gp=gpar(fontsize=12)))


pdf("Normal tissue expression/proteoglycan_hpa_cell_types_normal_heatmap_emt_tf_germ_layer.pdf",  width=13,height=16)

ht <- Heatmap(h, col = col_fun,
        cluster_columns = Colv,
        name = "Expression Values",
        show_heatmap_legend = T,
        top_annotation = ha,
        show_column_names = T,
        show_row_names = T,cluster_rows = Rowv,
        row_title_gp = gpar(fontsize=12),
        column_title_gp = gpar(fontsize=12),
        height = unit(15, "cm"),
        width = unit(22, "cm"),
        column_dend_reorder = T,
        row_dend_reorder = T,
        show_row_dend = T,
        border = T,
        column_dend_height = unit(2, "cm"),
        column_names_rot = 90,
        row_names_gp = grid::gpar(fontsize= 10),
        column_names_gp = grid::gpar(fontsize = 10),
        heatmap_legend_param = list(title=" Expression (log2+1)", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=12, fontface="bold"),labels_gp = gpar(fontsize=12)))

draw(ht,column_title = "Proteoglycan in Normal Tissue Cell Types", column_title_gp = gpar(fontsize = 16), 
     merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))

dev.off()

###############
## Libraries ##
###############

library(tidyverse)
library(data.table)
library(ggpubr)
library(optparse)
library(survival)
library(survminer)

##################
## Kaplan-Meier ##
##################

create_km <- function(survival_fit, GENE_NAME, CANCER, analysis_type, OUT) {
  output_file <- paste(OUT, "survival_", analysis_type, "_" , CANCER, "_", GENE_NAME, ".pdf", sep = "")
  p <- ggsurvplot(fit = survival_fit, 
                  pval = TRUE, 
                  surv.median.line = "hv", legend = c(0.1,0.2),
                  xlab = paste(analysis_type, "(Years)"),
                  ylab = paste(analysis_type, "Probability"),
                  title = paste(GENE_NAME, analysis_type, "for", CANCER),
                  ylim = c(0.0, 1), 
                  xlim = c(0, 20),
                  palette = c("#E69F00", "#0072B2"),
                  pval.coord = c(0.05, 0.8), 
                  break.x.by = 5,         
                  conf.int = TRUE,
                  risk.table = TRUE, risk.table.title = "",
                  risk.table.height = 0.15,
                  ncensor.plot = TRUE,
                  ncensor.plot.height = 0.15,
                  legend.labs = c("High", "Low"),
                  legend.title = paste(GENE_NAME), 
                  tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                                       axis.title.x = element_blank(), axis.title.y = element_blank(),
                                       axis.text.y = element_text(size = 14)), font.x = 14, font.y = 14, font.tickslab = 14)

pdf(output_file, width = 8, height = 8, onefile = FALSE)
print(p)
dev.off()

}

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

GENE <- "figures/tileplot/expr_cutoff_genes_cancer_type_tcga.csv"

GENE_NAME <- "ACAN"

CANCER <- "ACC"

OUT <- "figures/survival/"

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
gene <- gene[,-1]
gene_list <- unique(gene$gene)

#gene_expr <- subset(expr, expr$gene %in% gene_list)
#rownames(gene_expr) <- gene_expr[,1]
#gene_expr$gene <- NULL

#gene_expr <- as.data.frame(t(gene_expr))
#gene_expr <- rownames_to_column(gene_expr, "id")
#merged <- merge(gene_expr, clin, by.x = "id", by.y = "sample")
#merged <- merged[,c(1:45,70,71,74,75)]
#saveRDS(merged,"survival_tcga.rds") #data for app
##################
## Prepare data ##
##################

# Initialize results data frames
dss_results <- data.frame()
pfi_results <- data.frame()


#Loop for each gene and cancer type 
for (GENE_NAME in gene_list) {
  for (CANCER in unique(clin$`cancer type abbreviation`)) {
  
    gene_expr <- subset(expr, expr$gene %in% gene_list)
    rownames(gene_expr) <- gene_expr[,1]
    gene_expr$gene <- NULL
    
    gene_expr <- as.data.frame(t(gene_expr))
    gene_expr <- rownames_to_column(gene_expr, "id") 
    
    gene_expr <- gene_expr[, c("id", GENE_NAME)]
    colnames(gene_expr)[2] <- "gene"
    
    merged <- merge(gene_expr, clin, by.x = "id", by.y = "sample")
    merged <- subset(merged, merged$`cancer type abbreviation` == CANCER)
    
    if (nrow(merged) == 0) next
    
    #Subset info file to obtain DSS and PFI 
    info_surv_DSS <- select(merged, c("id","gene", "cancer type abbreviation", "DSS", "DSS.time"))
    info_surv_PFI <- select(merged, c("id","gene", "cancer type abbreviation", "PFI", "PFI.time"))
    info_surv_DSS$cancer_type <- as.factor(info_surv_DSS$`cancer type abbreviation`)
    info_surv_PFI$cancer_type <- as.factor(info_surv_PFI$`cancer type abbreviation`)
    
    #Extract cutoff based on the gene and cancer type of interest
    cutoff <-  gene %>% filter(gene == GENE_NAME & cancer == CANCER) %>%
      pull(cutoff_value)
    
    if (length(cutoff) == 0) next
    
    #Define gene expression based on cutoff
    info_surv_DSS$gene_expression <- ifelse(info_surv_DSS$gene >= cutoff, 'High', "Low")
    info_surv_PFI$gene_expression <- ifelse(info_surv_PFI$gene >= cutoff, 'High', "Low")
    
    if (all(is.na(info_surv_DSS$gene_expression)) || all(is.na(info_surv_PFI$gene_expression))) next
    
    
    info_surv_DSS <- subset(info_surv_DSS, info_surv_DSS$DSS.time > 0)
    info_surv_PFI <- subset(info_surv_PFI, info_surv_PFI$PFI.time > 0)
    
    if (nrow(info_surv_DSS) == 0 || nrow(info_surv_PFI) == 0) next
    
    #convert time into years
    info_surv_DSS$years <- info_surv_DSS$DSS.time/365
    info_surv_PFI$years <- info_surv_PFI$PFI.time/365
    
    info_surv_DSS <- info_surv_DSS[complete.cases(info_surv_DSS),]
    info_surv_PFI <- info_surv_PFI[complete.cases(info_surv_PFI),]
    
    # Check if the gene expression groups have sufficient observations
    if (length(unique(info_surv_DSS$gene_expression)) < 2 || length(unique(info_surv_PFI$gene_expression)) < 2) next
    
    #Define survival
    survival_dss = Surv(time= info_surv_DSS$years, event = info_surv_DSS$DSS)
    survival_pfi = Surv(time= info_surv_PFI$years, event = info_surv_PFI$PFI)
    
    # Check if the gene expression groups have sufficient observations
    if (length(unique(info_surv_DSS$gene_expression)) < 2 || length(unique(info_surv_PFI$gene_expression)) < 2) next
    
    survival_fit_dss <- survfit(formula = survival_dss ~ info_surv_DSS$gene_expression, data = info_surv_DSS)
    survival_fit_pfi <- survfit(formula = survival_pfi ~ info_surv_PFI$gene_expression, data = info_surv_PFI)
    
    #Cox hazard model
    res.cox_dss <- coxph(Surv(years, DSS) ~ gene_expression, data = info_surv_DSS)
    res.cox_pfi <- coxph(Surv(years, PFI) ~ gene_expression, data = info_surv_PFI)
    
    # Store results
    dss_summary <- summary(res.cox_dss)
    pfi_summary <- summary(res.cox_pfi)
    
    dss_results <- rbind(dss_results, data.frame(Gene = GENE_NAME, Cancer = CANCER, dss_summary$coefficients))
    pfi_results <- rbind(pfi_results, data.frame(Gene = GENE_NAME, Cancer = CANCER, pfi_summary$coefficients))
    
    #create KM plots if Cox results are significant 
    
    if(!is.na(dss_summary$coefficients[5]) && dss_summary$coefficients[5] < 0.05) {
      create_km(survival_fit_dss, GENE_NAME, CANCER, "DSS", OUT)
    }
    if(!is.na(pfi_summary$coefficients[5]) && pfi_summary$coefficients[5] < 0.05) {
      create_km(survival_fit_pfi, GENE_NAME, CANCER, "PFI", OUT)
    }
  }
}


write.csv(dss_results, paste(OUT, "cox_univariate_dss_results.csv", sep = ""), row.names = TRUE)
write.csv(pfi_results, paste(OUT, "cox_univariate_pfi_results.csv", sep = ""), row.names = TRUE)






#############
# Libraries #
#############

library(tidyverse)
library(data.table)

##########################
### TCGA preprocessing ###
##########################

# make one file for both the survival and the tcga clinical data 
#load clinical survival data 
clin_surv <- as.data.frame(fread("Survival_SupplementalTable_S1_20171025_xena_sp"))

#load clinical tcga data
clin_tcga <- as.data.frame(fread("TCGA clinical data with abbreviations.csv"))  

#merge the two files 
clin_merged <- merge(clin_surv, clin_tcga, by = "sample")

#look for duplicates --> no duplicates found
dup_test <- as.data.frame(duplicated(clin_merged$sample))

## remove normal samples 
noncancer <- subset(clin_merged, clin_merged$sample_type == "Solid Tissue Normal")
# 1413 obs of non cancer

#remove the non-cancer samples --> 11178 samples now
clin_merged <- subset(clin_merged, clin_merged$sample_type != "Solid Tissue Normal")

table(clin_merged$`cancer type abbreviation`)

#save as a file to use further 
fwrite(clin_merged, "clinical_tcga_data.csv", append=TRUE, row.names = TRUE)


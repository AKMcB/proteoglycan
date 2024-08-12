
#############
# Libraries #
#############

library(tidyverse)
library(data.table)


######################
# CCLE preprocessing #
######################

#Read clinical file 
clin_ccle <- as.data.frame(fread("clinical_data_ccle_edited.csv"))
head(clin_ccle)

#to only get the data we need for further analysis
clin_ccle <- clin_ccle[,c("ModelID","DepmapModelType", "OncotreeCode", "OncotreeSubtype","OncotreePrimaryDisease")]


#look at all unique codes. some codes are the same as TCGA, others are not
unique(clin_ccle$OncotreeCode)  

# has to reload the ccle and the tcga clinical data after editing the files in exel 
# has to match the abbreviations

#### load in clinical data tcga (edited)
clin_tcga <- as.data.frame(fread("clinical_data_TCGA_edited.csv"))

#match the clinical data in ccle and tcga 
clin_ccle <- clin_ccle %>% arrange(OncotreePrimaryDisease)
clin_tcga <- clin_tcga %>% arrange(`_primary_disease`)

# change name so it fits with ccle name
colnames(clin_tcga)[4] <- "OncotreePrimaryDisease"

#chnage tcga so we only get primary disease and the codes
clin_tcga <- clin_tcga[,c(4,5)]


#join tcga and ccle together so we get the tcga code in the ccle file
df1 <- distinct(clin_ccle, clin_ccle$OncotreePrimaryDisease,.keep_all = TRUE)

df1$`clin_ccle$OncotreePrimaryDisease` <- NULL

df2 <- distinct(clin_tcga, clin_tcga$`tcga code`, .keep_all = TRUE)

df2$`clin_tcga$\`tcga code\`` <- NULL

comb <- merge(df1, df2, by = "OncotreePrimaryDisease")

comb <- comb[,c("OncotreePrimaryDisease", "tcga code")]


clin_ccle <- merge(clin_ccle, comb, by = "OncotreePrimaryDisease")
table(clin_ccle$`tcga code`)

#save as a file to use further 
fwrite(clin_ccle, "clinical_ccle_data.csv", append=TRUE, row.names = TRUE)

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


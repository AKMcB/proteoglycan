
#############
# libraries #
#############

library(tidyverse)
library(data.table)


######################
# CCLE preprocessing #
######################

#Read expression data
expr <- as.data.frame(fread("OmicsExpressionProteinCodingGenesTPMLogp1.csv", header=TRUE))
colnames(expr)[1] <- "id"

expr <- column_to_rownames(expr, "id")

new_colnames <- gsub("\\ .*","",colnames(expr))
colnames(expr) <- new_colnames

expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, "gene")

#Check for duplicates 

dup <- as.data.frame(duplicated(expr$gene)) #All FALSE values

#save file
fwrite(expr, "expression_data_ccle_edited.csv", row.names = TRUE)



######################
# TCGA preprocessing #
######################

##Read expression file
expr <- as.data.frame(fread("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"))

#check for gene duplicates
dup <- expr[duplicated(expr$sample),]
expr <- expr[!expr$sample == "SLC35E2",]

rownames(expr) <- expr[,1]
expr$sample <- NULL


x <- as.data.frame(names(expr))

#check for duplicates 
dup <- substr(names(expr),start = 1,stop = 15)

# gets 11060 that is unique, meaning 9 out of the 11069 is duplicates
length(unique(dup))

#look at the duplicated ones
dup <- dup[duplicated(dup)]
dup #9 duplicates 

dup_1 <- expr[,c(6702:6703)]


# identify the duplicates
#remove both instances of them (18 samples)
cols_remove <- c("TCGA-21-1076-01", "TCGA-21-1076-01.1",
                 "TCGA-DD-AACA-02", "TCGA-DD-AACA-02.1",
                 "TCGA-06-0156-01", "TCGA-06-0156-01.1",
                 "TCGA-06-0211-01", "TCGA-06-0211-01.1",
                 "TCGA-DU-6404-02", "TCGA-DU-6404-02.1",
                 "TCGA-DU-6407-02", "TCGA-DU-6407-02.1",
                 "TCGA-FG-5965-02", "TCGA-FG-5965-02.1",
                 "TCGA-TQ-A7RK-02", "TCGA-TQ-A7RK-02.1",
                 "TCGA-23-1023-01", "TCGA-23-1023-01.1")

# remove the duplicates
expr <- expr[, !(names(expr) %in% cols_remove), drop = FALSE]
# has the 11051 samples now 
dup <- substr(names(expr),start = 1,stop = 15)
length(unique(dup)) #recheck duplicates

x <- as.data.frame(names(expr))
# look at the content of expr 
# -01 = 9684
# -02 = 33 
# -03 = 178
# -04 = 0
# -05 = 9
# -06 = 396
# -07 = 1
# -08 = 0
# -09 = 0
# -10 = 0
# -11 = 737 
# total = ~11051 

# remove samples with 11, because they are Solid Tissue Normal --> 737
# before removing: 11051
# should be 10323 after removing
expr <- expr[, !grepl("-11$", colnames(expr))] # now: 10314


#check for duplicated genes 
dup <- as.data.frame(duplicated(rownames(expr))) #no gene duplicates

#chnage NAs and negative values to zero 
test <- expr %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)))

test[test<0] <- 0

expr <- test

#save file 
fwrite(expr, "expression_data_tcga_log2.csv", row.names = TRUE)



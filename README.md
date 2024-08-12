Proteoglycan project 

1) Preprocess expression file and clinical file
     preprocessing_expression_data.R & preprocessing_clinical_sample_data.R
2) Decide expression cutoff to filter out potential housekeeping genes
     deciding_cutoff_tcga.R 
3) Create overview of genes in each cancer type above cutoff value
     tcga_filter_variance.R
4) Decide cutoff expression value to determine patients in high vs low expression groups
     deciding_cutoff_expr.R
5) Use the cutoff values in a Cox univariate analysis to decide if the expression has an impact on either DSS or PFI
     survial_cox_univariate.R & dotplot_survival.R


![pipeline](https://github.com/user-attachments/assets/aa6ab836-5ec3-4af8-b997-6f330f9afd29)

# **Proteoglycan project: PanCan analyses of proteoglycan expression in survival (DSS & PFI)**

1) Preprocess expression file and clinical file : *preprocessing_expression_data.R* & *preprocessing_clinical_sample_data.R*
2) Decide expression cutoff to filter out potential housekeeping genes: *deciding_cutoff_tcga.R* 
3) Create overview of genes in each cancer type above cutoff value: *tcga_filter_variance.R*
4) Decide cutoff expression value to determine patients in high vs low expression groups: *deciding_cutoff_expr.R*
5) Use the cutoff values in a Cox univariate analysis to decide if the expression has an impact on either DSS or PFI: *survial_cox_univariate.R*, *dotplot_survival.R* & *compare_pfi_dss.R*


![pipeline](https://github.com/user-attachments/assets/aa6ab836-5ec3-4af8-b997-6f330f9afd29)


The PanCan expression file can be downloaded from: <https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fucscpublic.xenahubs.net> 
Look in data folder for additional files used in the analysis


# **Human Protein Atlas**

The expression levels of genes in normal cell types and tissue can be downloaded from <https://www.proteinatlas.org/about/download>.
The following scripts *heatmap_hpa.R* and *proteoglycan_normal_hpa_tissue.R* were used to investigate the expression of the proteoglycan genes in individual cell types and tissue. 
The script *germ_layer.R* was used to categorize the cells types into their germ layer origin. 

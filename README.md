Deconvolution with BayesPrism and HR survival analysis:
mkdir result
make -f mk_deconvolution_survival rds=BC_Minor_tumor.rds outdir=result tumor=TCGA-BLCA count=TCGA-BLCA_406_gene_counts.csv meta=TCGA-BLCA.metaMatrix.RNA.csv idents=Minor deconvolution_survival

Notice: For the count file and metaMatrix file,please refer to script/TCGA_download_readme.sh


Definition of HR:
As for the other measures of association, a hazard ratio of 1 means lack of association, a hazard ratio greater than 1 suggests an increased risk, and a hazard ratio below 1 suggests a smaller risk。


References：
1.Chu, T., Wang, Z., Pe’er, D. et al. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. Nat Cancer 3, 505–517 (2022). https://doi.org/10.1038/s43018-022-00356-3
2.Chen, Z., Zhou, L., Liu, L. et al. Single-cell RNA sequencing highlights the role of inflammatory cancer-associated fibroblasts in bladder urothelial carcinoma. Nat Commun 11, 5077 (2020). https://doi.org/10.1038/s41467-020-18916-5

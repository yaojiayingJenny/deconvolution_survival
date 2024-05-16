1) Deconvolution with BayesPrism and HR survival analysis:  
mkdir result  
make -f mk_deconvolution_survival rds=BC_Minor_tumor.rds outdir=result tumor=TCGA-BLCA count=TCGA-BLCA_406_gene_counts.csv meta=TCGA-BLCA.metaMatrix.RNA.csv idents=Minor deconvolution_survival  

  Notice: For the count file and metaMatrix file,please refer to script/TCGA_download_readme.sh  


2) Definition of HR:  
As for the other measures of association, a hazard ratio of 1 means lack of association, a hazard ratio greater than 1 suggests an increased risk, and a hazard ratio below 1 suggests a smaller risk。


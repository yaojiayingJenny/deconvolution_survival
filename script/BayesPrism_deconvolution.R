args=commandArgs(T)
        if (length(args) !=5){
                print ("Rscript thiscirpt.r <bulk.count matrix（row:gene,column:sample,基因名不可以重复）> <rdsfile> <Ident> <key_celltype,NULL> <outdir>")
                print ("Example : Rscript thiscirpt.r TCGA_gene_counts.csv rds_minor.rds Minor  key_celltype ./")
                q()
        }


#######---------------------------------------------------------------
#Run BayesPrism.
#日期:2024/2/29 作者：yaojiaying

library(ggpubr) #ggboxplot
library(ggplot2)
suppressWarnings(library(BayesPrism))
#http://192.168.2.202:6006/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/TET_PUBLIC/yaojiaying/WDL/RNA_deconvolution/script/BayesPrism-main/tutorial_deconvolution.html
# recommend the use of unnormalized and untransformed count data

bulk.matrixfile=args[1]
rdsfile=args[2]
Ident=args[3]
key_celltype=args[4]
outdir=args[5]

#The parameter key is a character in cell.type.labels that corresponds to the malignant cell type. Set to NULL if there are no malignant cells or the malignant cells between reference and mixture are from matched samples, in which case all cell types will be treated equally.
if(key_celltype=='NULL'){key_celltype<-NULL}


setwd(outdir)


bulk.matrix<-read.csv(bulk.matrixfile,header=T,row.names=1,check.names =FALSE)
bk.dat<-data.frame(t(bulk.matrix))#需要行是cell/sample ID,列是gene name


rds<-readRDS(rdsfile)
sc.dat<-data.frame(t(data.frame(rds@assays$RNA@counts)))#需要行是cell/sample ID,列是gene name
rds$celltype<-rds@meta.data[,Ident]
cell.type.labels<-rds$celltype
cell.state.labels<-rds$stim

#QC of cell type and state labels
pdf('Correlation_cell.state.labels.pdf',w=6,h=6)
plot.cor.phi (input=sc.dat,
                         input.labels=cell.state.labels,
                         title="cell state correlation",
                         #specify pdf.prefix if need to output to pdf
                         #pdf.prefix="gbm.cor.cs", 
                         cexRow=0.5, cexCol=0.5)
dev.off()

pdf('Correlation_cell.type.labels.pdf',w=6,h=6)
plot.cor.phi (input=sc.dat, 
                         input.labels=cell.type.labels, 
                         title="cell type correlation",
                         #specify pdf.prefix if need to output to pdf
                         #pdf.prefix="gbm.cor.ct",
                         cexRow=0.5, cexCol=0.5,
                         )
dev.off()

#Filter outlier genes
# pdf('variable_genes_outlier_scRNA.pdf',w=6,h=4)
# sc.stat <- plot.scRNA.outlier(
  # input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  # cell.type.labels=cell.type.labels,
  # species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  # return.raw=TRUE #return the data used for plotting. #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
# )
# dev.off()


# pdf('variable_genes_outlier_bulkRNA.pdf',w=6,h=4)
# bk.stat <- plot.bulk.outlier(
  # bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
    # sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  # cell.type.labels=cell.type.labels,
  # species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  # return.raw=TRUE#pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
# )
# dev.off()

#c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY")
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix", #"GEP","count.matrix"
                                    species="hs", 
                                    gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1") ,#
                                    exp.cells=5)
									

#only work for human
pdf('plot.bulk.vs.sc_corr.pdf',w=12,h=4)
plot.bulk.vs.sc(sc.input = sc.dat.filtered,
                            bulk.input = bk.dat
                            #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)
dev.off()

sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                        gene.type = "protein_coding")


myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", #"GEP","count.matrix"
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.type.labels,
  key=key_celltype,
  outlier.cut=0.01,
    outlier.fraction=0.1,
)

bp.res <- run.prism(prism = myPrism, n.cores=10)

saveRDS(bp.res,file='bp.res_result.rds')


# extract posterior mean of cell type fraction theta
theta <- data.frame(get.fraction (bp=bp.res,
            which.theta="final",
            state.or.type="type"))
# meta<-rds@meta.data
# theta$group<-as.vector(meta[rownames(theta),]$group)

# celllist<-names(table(rds$celltype))
# theta<-theta[,celllist]

write.csv(theta,file='fraction_result.csv',quote=F)


celllist<-colnames(theta)
dataset<-data.frame()
for (cell in celllist){
print(cell)
da<-data.frame(Sample=rownames(theta),Ratio=theta[,c(cell)])
da$celltype<-cell
dataset<-rbind(dataset,da)
}

my36colors <-c('#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87','#E95C59','#E59CC4','#AB3282','#23452F','#BD956A','#8C549C','#585658','#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3', '#E4C755','#F7F398', '#AA9A59','#E63863','#E39A35','#C1E6F3', '#6778AE','#91D0BE','#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6','#625D9E','#68A180','#3A6963','#968175')
celltypes_colors<-my36colors

p0<-theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
        legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size=15),
        axis.text.x = element_text(color="black",size=8,angle=30),
        axis.text.y = element_text(color="black",size=8),
        axis.title.x = element_text(face="plain", color="black",size=15),
        axis.title.y = element_text(face="plain", color="black",size=10))+theme(legend.position='none')
		
p<-ggboxplot(dataset, 'celltype', "Ratio",
          color = 'celltype', add = "jitter",size=0.01,
          xlab = " ", ylab = 'Predicted cell proportion',palette = celltypes_colors) + rotate_x_text()+p0
pdf('predicted_proportion_celltype_boxplot.pdf',w=6+0.2*length(celllist),h=2.5)
print(p)
dev.off()

write.table(dataset,file='predicted_proportion_by_celltype.xls',quote=F,sep='\t',row.names=F)

# extract posterior mean of cell type-specific gene expression count matrix Z  
#可以特定提出某一群细胞的表达矩阵
# Z.tumor <- get.exp (bp=bp.res,
                          # state.or.type="type",
                          # cell.name="Malignant cell")
# write.csv(Z.tumor,file='predicted_Malignant_tmp.csv',quote=F)

print("Finish BayesPrism deconvolution! ")
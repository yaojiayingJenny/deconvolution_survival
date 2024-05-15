step1: 查询需要下载的癌症类型的TCGA缩写（如TCGA-COAD，结直肠癌）： https://zhuanlan.zhihu.com/p/344353688
step2: 203登录节点运行参考以下脚本下载TCGA数据(RNAseq文件和meta临床信息文件)，直接进入交互界面/annoroad/data1/bioinfo/PMO/yaojiaying/anaconda3/envs/r4-base/bin/R
step3: 下载完数据后，参考整理提取TCGA count数据,生存count matrix矩阵文件。

####----------------------------------------------------------------------------------------
#step2：下载TCGA数据
library(GDCRNATools)
library(DT)
# COAD
# Colon adenocarcinoma
# 结肠癌
setwd('/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/B_TET-086/Data/public/TCGA')
#238登录节点下载完
id<-'TCGA-COAD' #TCGA 缩写： https://zhuanlan.zhihu.com/p/344353688

#下载bulk RNA单个样本的数据
gdcRNADownload(project.id = id,data.type  = 'RNAseq', directory=paste0(id,'/RNAseq'))

###生成癌症患者的metadata矩阵
metaMatrix.RNA <- gdcParseMetadata(project.id = id,
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)
###对矩阵样本进行筛选
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
##展示前6行
datatable(as.data.frame(metaMatrix.RNA[1:6,]), extensions = 'Scroller',
          options = list(scrollX = TRUE, deferRender = TRUE, scroller = TRUE))
write.table(metaMatrix.RNA,file=paste0(id,'.metaMatrix.RNA.csv'),quote=F,sep='\t',row.names=F)


####----------------------------------------------------------------------------------------
##step3：整理提取TCGA count/tmp数据
id<-'TCGA-COAD'
dataset<-read.csv('/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/B_TET-086/Data/public/TCGA/TCGA-COAD.metaMatrix.RNA.csv',header=T,sep='\t')
inpath='/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/B_TET-086/Data/public/TCGA/TCGA-COAD/RNAseq/'
data<-subset(dataset,sample_type=='PrimaryTumor') #SolidTissueNormal
da<-read.table(filename<-paste0(inpath,data[1,2],'/',data[1,1]),sep='\t',skip=6)
tmp<-data.frame(gene_id=da$V1,gene_name=da$V2,gene_type=da$V3)
#count
for (i in seq(1:nrow(data))){
filename<-paste0(inpath,data[i,2],'/',data[i,1])
print(i)
da<-read.table(filename,sep='\t',skip=6)
colnames(da)<-c('gene_id','gene_name','gene_type','unstranded','stranded_first','stranded_second','tpm_unstranded','fpkm_unstranded','fpkm_uq_unstranded')
#d<-as.vector(da[,'tpm_unstranded']) #tmp
d<-as.vector(da[,'unstranded']) #count
tmp<-cbind(tmp,d)
}
colnames(tmp)<-c('gene_id','gene_name','gene_type',as.vector(data$sample))
norm<-tmp
norm1<-norm[!duplicated(norm$gene_name),]
rownames(norm1)<-norm1$gene_name
norm2<-norm1[,as.vector(data$sample)] #相当于去掉'gene_id','gene_name','gene_type'
n<-dim(norm2)[2]
write.csv(norm2,file=paste0(id,'_',n,'_gene_counts.csv'),quote=F,row.names=T)





####----------------------------------------------------------------------------------------
#tmp， 可用于signature打分，生成seurat rds
da<-read.table(filename<-paste0(inpath,data[1,2],'/',data[1,1]),sep='\t',skip=6)
tmp<-data.frame(gene_id=da$V1,gene_name=da$V2,gene_type=da$V3)
for (i in seq(1:nrow(data))){
filename<-paste0(inpath,data[i,2],'/',data[i,1])
print(i)
da<-read.table(filename,sep='\t',skip=6)
colnames(da)<-c('gene_id','gene_name','gene_type','unstranded','stranded_first','stranded_second','tpm_unstranded','fpkm_unstranded','fpkm_uq_unstranded')
d<-as.vector(da[,'tpm_unstranded']) #TPM
#d<-as.vector(da[,'unstranded']) #count
tmp<-cbind(tmp,d)
}
colnames(tmp)<-c('gene_id','gene_name','gene_type',as.vector(data$sample))
norm<-tmp
norm1<-norm[!duplicated(norm$gene_name),]
rownames(norm1)<-norm1$gene_name
norm2<-norm1[,as.vector(data$sample)] #相当于去掉'gene_id','gene_name','gene_type'
n<-dim(norm2)[2]
norm3<-log2(norm2+1)
write.csv(norm3,paste0(id,'_',n,'_gene_log2TPM.csv'),quote=F,row.names=F)
genelist<-rownames(norm3)
write.table(genelist,file='genelist.xls',quote=F,row.names=F,col.names=F)

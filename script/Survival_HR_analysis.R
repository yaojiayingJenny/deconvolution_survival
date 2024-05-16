args=commandArgs(T)
        if (length(args) !=4){
                print ("Rscript thiscirpt.r <metaMatrix.RNA.csv（TCGA患者临床信息，需要有vital_status以及存活时间等）> <predicted_proportion file> <outdir> <TCGA-name>")
                print ("Example : Rscript thiscirpt.r TCGA-BLCA.metaMatrix.RNA.csv predicted_proportion_by_celltype.xls  ./ TCGA-BLCA")
                q()
        }

#######---------------------------------------------------------------
#Survival_HR_analysis.
#日期:2024/2/29 作者：yaojiaying
library("survival") 
library("survminer")
library("ggplot2") 
library("ggpubr") 
require("ggrepel") #geom_text_repel

#细胞占比生成分析

metafile=args[1]
predictedfile=args[2]
outdir=args[3]
tumor<-args[4]


mkdirs <- function(outdir,fp) {
  if(!file.exists(file.path(outdir,fp))) {
    #mkdirs(dirname(fp))
    dir.create(file.path(outdir,fp))}
  else{print(paste(fp,"Dir already exists!",sep="     "))
    unlink(file.path(outdir,fp), recursive=TRUE)
    dir.create(file.path(outdir,fp))}
}

setwd(outdir)

meta<-read.table(metafile,header=T,sep='\t')
dataset<-read.table(predictedfile,header=T,sep='\t')


data<-subset(meta,(sample %in% unique(dataset$Sample)) & vital_status %in% c('Alive','Dead') ) #过滤掉无vital_status的样本，只保留PrimaryTumor
data$status[data$vital_status=='Alive']=1
data$status[data$vital_status=='Dead']=2
data$time[data$vital_status=='Alive']=data$days_to_last_follow_up[data$vital_status=='Alive']
data$time[data$vital_status=='Dead']=data$days_to_death[data$vital_status=='Dead']
data$time<-as.numeric(data$time)

mkdirs(outdir,'survival_curve')
setwd(paste0(outdir,'/survival_curve/'))

celllist<-unique(dataset$celltype)
HRda<-data.frame()
for (name in celllist){
	da<-subset(dataset,celltype==name)
	data$Ratio<-plyr::mapvalues(x =data$sample,from = da$Sample,to =da$Ratio)
	data$Ratio<-as.numeric(data$Ratio)
	data1<-data
	#data1<-subset(data,vital_status!='Not Reported')
	data1$group[data1$Ratio<=median(data1$Ratio,na.rm=T)]='low'
	data1$group[data1$Ratio>median(data1$Ratio,na.rm=T)]='high'
	library(ggpubr) #ggboxplot
	p0<-theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
			legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size=15),
			axis.text.x = element_text(color="black",size=10,angle=30),
			axis.text.y = element_text(color="black",size=10),
			axis.title.x = element_text(face="plain", color="black",size=15),
			axis.title.y = element_text(face="plain", color="black",size=10))+theme(legend.position='none')
	w_h=c(4,3)
	dtres<-coxph(Surv(time,status)~group,data=data1)
	summary(dtres)
	#接着用Surv()函数创建生存数据对象（主要输入生存时间和状态逻辑值），再用survfit()函数对生存数据对象拟合生存函数，创建KM(Kaplan-Meier)生存曲线
	fit<- survfit(Surv(time,status)~group,data=data1)
	#计算log rank p和HR
	data.survdiff <- survdiff(Surv(time, status)~group, data=data1)
	p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
	HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
	up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
	low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
	a<-c(HR,low95,up95,p.val,name)
	HRda<-rbind(HRda,a)
	#
	ggsurvplot(fit, data = data1, 
			  pval=T,legend.title=name,
			  legend.labs=c("high","low"),
			  palette=c("#FF0000","#0000FF"),
			  xlab="Follow up time(d)")+ggtitle(paste0("Survival analysis with deconvoluted ",name,' ratio'))  		  
	ggsave(filename = paste0(name,"_deconvolution_bulk_Survival_curve_",tumor,".pdf"),width =6,height =5) 

}

setwd(outdir)
colnames(HRda)<-c('HR','low95','up95','p.val','celltype')
#绘制boxplot图

HRda$HR<-as.numeric(HRda$HR)
HRda$low95<-as.numeric(HRda$low95)
HRda$up95<-as.numeric(HRda$up95)
HRda$p.val<-as.numeric(HRda$p.val)

da<-HRda[,c('celltype','HR','low95','up95','p.val')]
write.table(da,file='deconvolution_bulk_Survival_HR_CI_dataset.xls',quote=F,sep='\t',row.names=F)
p0<-theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
			legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size=15),
			axis.text.x = element_text(color="black",size=10,angle=0),
			axis.text.y = element_text(color="black",size=10),
			axis.title.x = element_text(face="plain", color="black",size=15),
			axis.title.y = element_text(face="plain", color="black",size=10))+theme(legend.position='none')
my36colors <-c('#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87','#E95C59','#E59CC4','#AB3282','#23452F','#BD956A','#8C549C','#585658','#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3', '#E4C755','#F7F398', '#AA9A59','#E63863','#E39A35','#C1E6F3', '#6778AE','#91D0BE','#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6','#625D9E','#68A180','#3A6963','#968175')
celltypes_colors<-my36colors
pdf('deconvolution_bulk_Survival_curve_HR_CI_plot.pdf',w=4,h=3+0.2*length(celllist))
ggplot(da,aes(x=reorder(celltype, HR),y=HR,color=celltype))+geom_point()+theme_bw()+p0+xlab('')+ylab('HR(high) and 95% CI')+geom_hline(yintercept=1,linetype='dashed')+geom_errorbar(aes(x=celltype, ymax=up95, ymin=low95), width = 0.2)+geom_text_repel(data=subset(da,p.val<0.05), aes(x=celltype,label=paste0('p = ',round(p.val,4))),size=2.5,direction="both",min.segment.length = 0.05,segment.alpha=0.6,label.padding = 0.4,max.overlaps =30,size=2,angle = 0)+ggtitle(tumor)+coord_flip() +scale_color_manual(values=celltypes_colors)
dev.off()

print("Finish Survival_HR_analysis! ")

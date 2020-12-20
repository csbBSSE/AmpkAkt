#This code evaluates correlation values between AMPK and AKT phospho-levels and 
#generates scatter plots from RPPA data

setwd("A:/Ph.D/Manuscripts/AMPK-AKT_modelling/clinical_data")
rppa_file_name="TCGA-PANCAN32-L4/TCGA-PANCAN32-L4.csv"
    #"TCGA-SARC-L4/TCGA-SARC-L4.csv";
    #"TCPA-TCGA-BRCA-L4/TCGA-BRCA-L4.csv";

rppa_data=read.csv(rppa_file_name,header=T,stringsAsFactors = F)

#### scatter plot for AMPK pT172 and AKT pT308 ####
# selecting AMPK pT172 and AKT pT308
#for pan-cancer
dat1=rppa_data[,c(13,15)]
#for other cancer types 
#dat1=rppa_data[,c(14,16)]

cor_p=cor.test(dat1[,1],dat1[,2])

scat_name_png=paste(dirname(rppa_file_name),.Platform$file.sep,"Scatterplot_of_",colnames(dat1)[2],"_",colnames(dat1)[1],"_",sub(".csv",".png",basename(rppa_file_name)),sep="");
#scat_name_svg=paste(dirname(rppa_file_name),.Platform$file.sep,"Scatterplot_of_",colnames(dat1)[2],"_",colnames(dat1)[1],"_",sub(".csv",".svg",basename(rppa_file_name)),sep="");

#svg(scat_name_svg)
x11()
par(mar=c(6,5,5,4),xpd=T,cex.lab=1.5,cex.axis=1.5,cex.main=1.4,font.main=2,font.lab=2,font.axis=2)
plot(dat1[,2],dat1[,1],pch=20,col="red2",xlab = colnames(dat1)[2],ylab = colnames(dat1)[1])
title(main = "Scatter plot of AMPK pT172 and AKT pT308 levels",          
      sub=paste("pearson's R ",round(cor_p$estimate,3)," p-value ",signif(cor_p$p.value,3)," n = ",nrow(dat1),sep=""));
abline(h=0,xpd=F,col="blue",lty=3,lwd=3)
abline(v=0,xpd=F,col="blue",lty=3,lwd=3)
#dev.off() #SVG 
savePlot(scat_name_png,type = "png")
dev.off()

#### scatter plot for AMPK pT172 and AKT pS473 ####
# selecting AMPK pT172 and AKT pS473
#for pan-cancer
dat1=rppa_data[,c(12,15)]
#for other cancer types
dat1=rppa_data[,c(13,16)]

cor_p=cor.test(dat1[,1],dat1[,2])

scat_name_png=paste(dirname(rppa_file_name),.Platform$file.sep,"Scatterplot_of_",colnames(dat1)[2],"_",colnames(dat1)[1],"_",sub(".csv",".png",basename(rppa_file_name)),sep="");
#scat_name_svg=paste(dirname(rppa_file_name),.Platform$file.sep,"Scatterplot_of_",colnames(dat1)[2],"_",colnames(dat1)[1],"_",sub(".csv",".svg",basename(rppa_file_name)),sep="");

#svg(scat_name_svg)
x11()
par(mar=c(6,5,5,4),xpd=T,cex.lab=1.5,cex.axis=1.5,cex.main=1.4,font.main=2,font.lab=2,font.axis=2)
plot(dat1[,2],dat1[,1],pch=20,col="red2",xlab = colnames(dat1)[2],ylab = colnames(dat1)[1])
title(main = "Scatter plot of AMPK pT172 and AKT pS473 levels",
      sub=paste("pearson's R ",round(cor_p$estimate,3)," p-value ",signif(cor_p$p.value,3)," n = ",nrow(dat1),sep=""));
abline(h=0,xpd=F,col="blue",lty=3,lwd=3)
abline(v=0,xpd=F,col="blue",lty=3,lwd=3)
#dev.off() #SVG
savePlot(scat_name_png,type = "png")
dev.off()


#### Pan cancer correlations ####
rppa_data[,2]<-factor(rppa_data[,2])
groups<-levels(rppa_data[,2])
group_cor_data<-data.frame(matrix(NA,nrow = length(groups),ncol=6))
colnames(group_cor_data)<-c("cancer_type","n","R_AMPKpT172_AKTpT308","pvalue_AMPKpT172_AKTpT308","R_AMPKpT172_AKTpS473","pvalue_AMPKpT172_AKTpS473")
group_cor_data[,1]<-groups

for (i in 1:nrow(group_cor_data))
{
    dat=rppa_data[which(rppa_data[,2]==group_cor_data[i,1]),c(12,13,15)]
    group_cor_data[i,2]=nrow(dat)
    cor_1<-cor.test(dat[,1],dat[,3])
    cor_2<-cor.test(dat[,2],dat[,3])
    group_cor_data[i,3]<-round(cor_1$estimate,5)
    group_cor_data[i,4]<-round(cor_1$p.value,5)
    group_cor_data[i,5]<-round(cor_2$estimate,5)
    group_cor_data[i,6]<-round(cor_2$p.value,5)
}
group_cor_file_name<-paste(dirname(rppa_file_name),.Platform$file.sep,"pairwise_correlations_",paste(colnames(dat),collapse = "_"),"_",basename(rppa_file_name),sep="");
write.table(group_cor_data,group_cor_file_name,row.names = F,col.names = T,sep = ",")

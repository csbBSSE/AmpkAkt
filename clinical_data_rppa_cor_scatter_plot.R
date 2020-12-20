library(stringi)
setwd("working/directory")
rppa_file_name=c("TCPA-TCGA-BRCA-L4/TCGA-BRCA-L4.csv");
rppa_data=read.csv(rppa_file_name,header=T,stringsAsFactors = F)
patient_id=rppa_data$Sample_ID;
patient_id=sapply(patient_id,function(x) stri_join(unlist(strsplit(x,split = "-"))[1:3],sep="-",collapse="-"))
rppa_data$patient_ID=patient_id
phenotype_file_name=c("TCGA_BRCA_clinicalMatrix.tsv");
phenotype_file_name<-path.expand(phenotype_file_name)
phenotype_data = read.delim(phenotype_file_name,header=T,stringsAsFactors = F)
table(phenotype_data$AJCC_Stage_nature2012)
sel_ids=intersect(rppa_data$patient_ID,phenotype_data$X_PATIENT)
rppa_sel_data=rppa_data[which(sel_ids %in% rppa_data$patient_ID),]
phenotype_sel_data<-phenotype_data[which(sel_ids %in% phenotype_data$X_PATIENT),]
merged_data=merge.data.frame(rppa_data,phenotype_data,by.x="patient_ID",by.y = "X_PATIENT")
table(merged_data$Node_Coded_nature2012)
table(merged_data$Sample_Type,merged_data$Node_Coded_nature2012)
normal_data=merged_data[which(merged_data$Sample_Type=="Normal"),]
primary_data=merged_data[which(merged_data$Sample_Type=="Primary"),]
metastatic_data=merged_data[which(merged_data$Sample_Type=="Metastatic"),]
primary_node_neg_data=primary_data[which(primary_data$Node_Coded_nature2012 =="Negative"),]
primary_node_pos_data=primary_data[which(primary_data$Node_Coded_nature2012 =="Positive"),]



#scatter plots
dat1=primary_data[,c(15,17)]
cor_p=cor.test(dat1[,1],dat1[,2])
x11()
par(mar=c(6,5,5,4),xpd=T,cex.lab=1.3,cex.axis=1.2,cex.main=1.4,font.main=2)
plot(dat1[,2],dat1[,1],pch=20,col="red2",xlab = colnames(dat1)[2],ylab = colnames(dat1)[1])
title(main = c("Scatter plot of AMPK pT172 and AKT pT308 levels",
               paste("TCGA BRCA, n = ",nrow(dat1),sep="")),
      sub=paste("pearson's R ",round(cor_p$estimate,3)," p-value ",signif(cor_p$p.value,3),sep=""));
abline(h=0,xpd=F,col="blue",lty=3,lwd=2)
abline(v=0,xpd=F,col="blue",lty=3,lwd=2)
scat_name_png=paste(dirname(rppa_file_name),.Platform$file.sep,"Scatter_plot_of_",colnames(dat1)[2],"_",colnames(dat1)[1],"_",sub(".csv",".png",basename(rppa_file_name)),sep="");
scat_name_jpg=paste(dirname(rppa_file_name),.Platform$file.sep,"Scatter_plot_of_",colnames(dat1)[2],"_",colnames(dat1)[1],"_",sub(".csv",".jpg",basename(rppa_file_name)),sep="");
savePlot(scat_name_png,type = "png")
savePlot(scat_name_jpg,type = "jpg")
dev.off()
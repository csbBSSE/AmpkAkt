cran_packages<-c("RColorBrewer","ggplot2","ggpubr","BiocManager")
bioconductor_packages<-c("ComplexHeatmap")
to_install_cran<-cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
to_install_BioC<-cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(to_install_cran)){ install.packages(to_install_cran,repos = "https://cloud.r-project.org")}
library(BiocManager)
if(length(to_install_BioC)){BiocManager::install(to_install_BioC)}
all_packages<-c(cran_packages,bioconductor_packages)
sapply(all_packages,library,character.only = T)

#loading required parameter and final steady state files
setwd("/working/directory")
file_name="/path/to/full/parameter/sets/file";
parameter_sets=read.delim(file_name,header = T,stringsAsFactors = F);

file_name1="/path/to/full/absolute/steady/state/file";
data1=read.delim(file_name1,header = T,stringsAsFactors = F);
dat=apply(data1,2,as.numeric);
colnames(dat)<-sub("f_","",colnames(dat));
scale_dat<-dat[,2:5]
m=mean(scale_dat)
s=sd(scale_dat)
dat1=apply(scale_dat,2,function(x) (x-m)/s)
dat2=cbind(dat[,1],dat1)
colnames(dat2)<-colnames(dat)

#### scatter plots for AMPK and AKT final steady states ####
cor_p=cor.test(dat1[,1],dat1[,2])
x11()
par(mar=c(6,5,5,4),xpd=T,cex.lab=1.5,cex.axis=1.5,cex.main=1.4,font.main=2,font.lab=2,font.axis=2)
plot(dat1[,3],dat1[,1],pch=20,col="red2",#main = "Scatter plot of AMPK and AKT Z-scores",
     #xlab = colnames(dat1)[3],ylab = colnames(dat1)[1])
     xlab = paste(colnames(dat1)[2],"(Active level)"),ylab = paste(colnames(dat1)[1],"(Active level)"))
title(sub=paste("pearson's cor value ",round(cor_p$estimate,3)," p-value ",signif(cor_p$p.value,3),sep=""));
abline(h=0,xpd=F,col="blue",lty=3,lwd=3)
abline(v=0,xpd=F,col="blue",lty=3,lwd=3)
scat_name_png=paste(dirname(file_name),.Platform$file.sep,"Scatterplot_of_",colnames(dat1)[3],"_",colnames(dat1)[1],"_",sub(".txt",".png",basename(file_name)),sep="");
#scat_name_svg=paste(dirname(file_name),.Platform$file.sep,"Scatterplot_of_",colnames(dat1)[3],"_",colnames(dat1)[1],"_",sub(".txt",".svg",basename(file_name)),sep="");

savePlot(scat_name_png,type = "png")
#svg(scat_name_svg)
dev.off()

#### kac/kdac analysis ####  
kac<-parameter_sets[,6:9] # k_ac parameters 
kdac<-parameter_sets[,10:13] # k_dac parameters 
kac_kdac<-kac/kdac # ratio of k_ac/k_dac
kac_kdac<-cbind(parameter_sets[,1],kac_kdac) #combining with parameter sets
colnames(kac_kdac)<-colnames(dat)
colnames(kac_kdac)[2:5]<-paste("k_",colnames(kac_kdac)[2:5],sep = "")
k_par_sets<-merge(dat2,kac_kdac)
#categorizing final steady states
k_par_sets["groups"]<-NA
for (i in 1:nrow(k_par_sets))
{
    if(k_par_sets[i,2] <=0 & k_par_sets[i,3] <= 0)
    {
        k_par_sets[i,10]<-"LL"
    } else if(k_par_sets[i,2] >=0 & k_par_sets[i,3] <= 0)
    {
        k_par_sets[i,10]<-"HL"
    }else if(k_par_sets[i,2] <=0 & k_par_sets[i,3] >= 0)
    {
        k_par_sets[i,10]<-"LH"
    }else
    {
        k_par_sets[i,10]<-"HH"
    }
}
k_par_sets[,10]<-factor(k_par_sets[,10],ordered = T,levels = c("LL","LH","HL","HH"))

table(k_par_sets[,10])

#### heatmap ####
data2<-k_par_sets[,c(2,4,3,5,10)]
row.names(data2)<-paste("s",row.names(data2),sep="_")
data2$groups<-factor(data2$groups,ordered = T)
annotate_row=data.frame(group=factor(data2$groups,ordered = T))
plot_data=data2[,c(3,2,1,4)]
plot_data<-apply(plot_data,2,as.numeric)

b=palette.colors(n=4,"Dark2")
names(b)<- c("HH","HL","LH","LL")
color_groups<-list(group = b)

plot_names_png=paste(dirname(file_name),.Platform$file.sep,"comp_heatmap_row_annotated_",sub(".txt",".png",basename(file_name)),sep = "");
#plot_names_svg=paste(dirname(file_name),.Platform$file.sep,"comp_heatmap_row_annotated_",sub(".txt",".svg",basename(file_name)),sep = "");
png(plot_names_png)
#svg(plot_names_svg)
par(xpd=T,mar=c(5,6,4,4))
Heatmap(plot_data,cluster_rows=T,cluster_columns = F,
        show_row_names = F,show_column_names = T,
        clustering_distance_rows="euclidean",row_dend_width = unit(15, "mm"),
        col = rev(colorRampPalette(brewer.pal(5,"PuOr"))(20)),
        cell_fun = (width = unit(25,"pt")),
        column_names_rot = 0,column_names_centered=T,
        left_annotation = rowAnnotation(df = annotate_row,col=color_groups))
dev.off()

#### scatter plots- 4 kac/kdac #### 
HH<-k_par_sets[which(k_par_sets$groups=="HH"),]
HL<-k_par_sets[which(k_par_sets$groups=="HL"),]
LH<-k_par_sets[which(k_par_sets$groups=="LH"),]
LL<-k_par_sets[which(k_par_sets$groups=="LL"),]

HH_plot<-ggplot(HH,aes(x=HH[,6],y=HH[,7]))+
    geom_point(size=3,colour = "#E7298A")+
    labs(title = "k ac/k dac for HH")+
    scale_x_continuous(name="k ac/k dac AMPK", limits = c(0,10))+ 
    scale_y_continuous(name="k ac/k dac AKT", limits = c(0,10))+
    theme(plot.title =  element_text(size = 15,face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11))

HL_plot<-ggplot(HL,aes(x=HL[,6],y=HL[,7]))+
    geom_point(size=3,colour = "#7570B3")+
    labs(title = "k ac/k dac for HL")+
    scale_x_continuous(name="k ac/k dac AMPK", limits = c(0,10))+ 
    scale_y_continuous(name="k ac/k dac AKT", limits = c(0,10))+
    theme(plot.title =  element_text(size = 15,face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11))

LH_plot<-ggplot(LH,aes(x=LH[,6],y=LH[,7]))+
    geom_point(size=3,colour = "#D95F02")+
    labs(title = "k ac/k dac for LH")+
    scale_x_continuous(name="k ac/k dac AMPK", limits = c(0,10))+ 
    scale_y_continuous(name="k ac/k dac AKT", limits = c(0,10))+
    theme(plot.title =  element_text(size = 15,face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11))

LL_plot<-ggplot(LL,aes(x=LL[,6],y=LL[,7]))+
    geom_point(size=3,colour = "#1B9E77")+
    labs(title = "k ac/k dac for LL")+
    scale_x_continuous(name="k ac/k dac AMPK", limits = c(0,10))+ 
    scale_y_continuous(name="k ac/k dac AKT", limits = c(0,10))+
    theme(plot.title =  element_text(size = 15,face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11))

k_ampk_groups_plot<-ggplot(k_par_sets,aes(x=k_par_sets[,10],y=k_par_sets[,6],fill=groups))+
    geom_violin(trim=T)+scale_fill_brewer(palette="Dark2")+
    labs(title = "k ac/k dac AMPK group statistics")+
    scale_x_discrete(name="Groups")+ 
    scale_y_continuous(name="k ac/k dac AMPK", limits = c(0,12))+
    stat_summary(fun.data=mean_sdl,geom="pointrange",color="black")+
    theme(plot.title =  element_text(size = 15,face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11),
          legend.position = "none")

k_akt_groups_plot<-ggplot(k_par_sets,aes(x=k_par_sets[,10],y=k_par_sets[,7],fill=groups))+
    geom_violin(trim=T)+scale_fill_brewer(palette="Dark2")+
    labs(title = "k ac/k dac AKT group statistics")+
    scale_x_discrete(name="Groups")+ 
    scale_y_continuous(name="k ac/k dac AKT", limits = c(0,12))+ 
    stat_summary(fun.data=mean_sdl,geom="pointrange",color="black")+
    theme(plot.title =  element_text(size = 15,face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11),
          legend.position = "none")

figure1<-ggarrange(LL_plot, HL_plot, k_ampk_groups_plot, LH_plot, HH_plot, k_akt_groups_plot,
                   labels = c("A", "C", "E", "B", "D", "F"),font.label = list(size = 17),
                   ncol = 3, nrow = 2)

#figure1
x11(width = 14,height = 7)
plot(figure1)
splot_name_png=paste(dirname(file_name),.Platform$file.sep,"Scatter_and violin_plots_Kac_Kdac_AMPK_AKT_",sub(".txt",".png",basename(file_name)),sep = "");
#splot_name_svg=paste(dirname(file_name),.Platform$file.sep,"Scatter_and violin_plots_Kac_Kdac_AMPK_AKT_",sub(".txt",".svg",basename(file_name)),sep = "");

ggsave(splot_name_png,device = "png")
#ggsave(splot_name_svg,device = "svg")
dev.off()

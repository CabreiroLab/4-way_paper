library(tidyverse)
library(ggrepel)

library(ComplexHeatmap)
library(circlize)

#library(RColorBrewer)


#library(limma)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)




library(ggthemes)

theme_Publication <- function(base_size=14) {

  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            #panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
           #legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            #plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

theme_set(theme_Publication())
#theme_set(theme_light())



pole<-function(x,y) {
  r<-sqrt(x^2+y^2)
  o<- -2*atan(y/(x))/(pi/2)
  return(o) 
}

pole2<-function(x,y) {
  o<- (atan(y/(x)))/(2*pi)
  return(o) 
}







cwd<-"~/Dropbox/Projects/2015-Metformin/Biolog_Met_NGM/"
setwd(cwd)

#load('Ecoli.RData')
#save.image('Ecoli.RData')


odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



info<-read_csv('../Biolog/Biolog_metabolites_EcoCyc_Unique_PM1-PM5.csv') %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A'))

head(info)


# info<-info %>%
#   mutate(Metabolite=trimws(as.character(Metabolite)),
#          Name=trimws(as.character(Name)),
#          MetaboliteU=ifelse(Metabolite=='Negative Control','NGM',Metabolite),
#          MetaboliteU=ifelse(Metabolite %in% dupls & Plate %in% c('PM3B'),
#                                 paste(Metabolite,'N',sep='_'),MetaboliteU),
#          MetaboliteU=ifelse(Metabolite %in% dupls & Plate %in% c('PM4A') & Group=='P-Source',
#                                 paste(Metabolite,'P',sep='_'),MetaboliteU),
#          MetaboliteU=ifelse(Metabolite %in% dupls & Plate %in% c('PM4A') & Group=='S-Source',
#                                 paste(Metabolite,'S',sep='_'),MetaboliteU),
#          MetaboliteU=ifelse(Index=='PM1-A1','NGM_C1',MetaboliteU),
#          MetaboliteU=ifelse(Index=='PM2A-A1','NGM_C2',MetaboliteU),
#          MetaboliteU=ifelse(Index=='PM3B-A1','NGM_N',MetaboliteU),
#          MetaboliteU=ifelse(Index=='PM4A-A1','NGM_P',MetaboliteU),
#          MetaboliteU=ifelse(Index=='PM4A-F1','NGM_S',MetaboliteU)) %>%
#   rename(Metabolite_class=Description)
# 
# 
# head(info)



#Get data

data<-read_csv('Data/Summary.csv') %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A')) %>%
  mutate(Sample=paste(Strain,ifelse(Type=='Control','C','T'),sep='_'),
         SampleID=paste(Sample,Replicate,sep='_'),
         Sample=factor(Sample,
                        levels=c("OP50Sens_C","OP50Sens_T"),
                        labels=c("OP50Sens_C","OP50Sens_T"))) %>%
  left_join(info[,c('Index','MetaboliteU')],by='Index') %>%
  rename(Metabolite=Name,G=Int_750nm_log,GR=a_log) %>%
  select(File:Metformin_mM,Replicate,Well,Index,Sample:MetaboliteU,Metabolite:Group,G,GR) %>%
  gather(key=Measure,value=Value,G,GR) %>%
  group_by(Plate,Type,Replicate,Group,Measure) %>%
  mutate(Value_ref=Value[Metabolite=='Negative Control'],
         Value_norm=Value-Value_ref) %>%
  ungroup %>%
  mutate_at(c('SampleID','Sample','Strain','Metformin_mM'),as.factor)


data.nc<-data %>%
  filter(Metabolite!='Negative Control' & Plate !='PM5')

PM1PM2<-subset(info,Plate %in% c('PM1','PM2A'))




bioinfo<-data.nc %>%
  group_by(SampleID,Sample,Strain,Metformin_mM) %>%
  summarise %>%
  ungroup %>%
  data.frame

rownames(bioinfo)<-bioinfo$SampleID



pcashape<- data.nc %>%
  filter(Measure=='G') %>%
  select(SampleID,Index,Value_norm) %>%
  spread(Index,Value_norm) %>%
  data.frame


rownames(pcashape)<-pcashape$SampleID
pcashape$SampleID<-NULL


hca_sample<-hclust(dist(pcashape,method="euclidean"),method="ward.D2")
#hca_variable<-hclust(dist(t(log_data),method="manhattan"),method="complete")

plot(hca_sample, labels=samples,
     hang=-1, main="Cluster Dendrogram", xlab="", sub="", cex=1)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Hierarchical_Clustering.pdf",sep=''),
             width=9,height=6)

ir.pca <- prcomp(clean_data,
                 center = TRUE,
                 scale. = TRUE)

write.csv(ir.pca[2],paste(odir,"/PCA_loadings.csv",sep=''))

plot(ir.pca,type='l')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_components.pdf",sep=''),
             width=9,height=6)

pcadata<-data.frame(ir.pca$x) %>%
  rownames_to_column('SampleID') %>%
  left_join(bioinfo)





pcaresult<-summary(ir.pca)$importance
PC1prc<-round(pcaresult['Proportion of Variance',][[1]]*100,0)
PC2prc<-round(pcaresult['Proportion of Variance',][[2]]*100,0)



ellipses<-pcadata %>%
  group_by(Sample,Strain,Metformin_mM) %>%
  do(getellipse(.$PC1,.$PC2,1) ) %>%
  data.frame
  
# 
# pcadata$Metformin_mM<-as.factor(pcadata$Metformin_mM)
# ellipses$Metformin_mM<-as.factor(ellipses$Metformin_mM)

ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Sample),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=5,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  #geom_text(aes(label=Sample))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA.pdf",sep=''),
             width=12,height=9)



head(data.n)

#Heatmap
heatshape<- data.nc %>%
  filter(Measure=='G') %>%
  select(SampleID,MetaboliteU,Value_norm) %>%
  spread(SampleID,Value_norm) %>%
  data.frame


rownames(heatshape)<-heatshape$MetaboliteU

heatshape$MetaboliteU<-NULL
head(heatshape)



#Order anotation by heatmap colnames
hanot<-bioinfo[colnames(heatshape),] %>%
  select(Metformin_mM)

ha<-HeatmapAnnotation(df=hanot, col = list(Metformin_mM=c('50'='black','0'='white')))

heatshape.sc<-t(scale(t(heatshape)))

rownames(heatshape.sc)<-NULL

#make anotated heatmap
Heatmap(heatshape.sc,name = 'Z-score',
        #col=bgg4,
        column_names_side = 'top',
        clustering_method_rows='ward.D2',
        #clustering_method_columns ='ward.D2',
        top_annotation = ha,
        row_dend_reorder = TRUE,
        column_dend_reorder = TRUE,
        row_names_max_width = unit(10, "cm"))


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap.pdf",sep=''),
             width=6,height=9, useDingbats=FALSE)





#Get LM ready
#Based on selected data
#Get negative control wells
head(data.nc)




contrasts<-read.contrasts2('!Contrasts.xlsx')

contrasts.desc<-contrasts$Contrasts.table %>%
  select(Description:Strain)

contr.matrix<-contrasts$Contrasts.matrix
contr.matrix

View(data.nc)


results.all<-data.nc %>%
  group_by(Measure,Group,Plate,Well,Index,Metabolite,MetaboliteU,EcoCycID,KEGG_ID) %>%
  do(hypothesise2(.,'Value_norm~0+Sample',contr.matrix)) %>%
  getresults(contrasts.desc)






data.sum<-data.nc %>%
  filter(Plate!='PM5' & Measure=='G') %>%
  group_by(Type,Index,MetaboliteU) %>%
  summarise(G_Mean=mean(Value),G_SD=sd(Value),G_SE=G_SD/sqrt(length(Value))) %>%
  ungroup %>%
  data.frame




results<-results.all$results
head(results)

logFDRbreaks<-c(-1,1.3,2,3,14)
logFDRbins<-c('N.S.','p<0.05','p<0.01','p<0.001')



results.castfull

.c<-dcast(subset(results.m,Stat %in% c('logFC','SE','FDR','NE','PE','logFDR')),Index+Plate+Well+Metabolite+MetaboliteU+KEGG_ID+EcoCycID+Group+Metabolite_class~Measure+Contrast+Stat,value.var = 'Value')
results.cast.c<-dcast(subset(results.m,Stat %in% c('logFC','FDR')),Index+Plate+Well+Metabolite+MetaboliteU+KEGG_ID+EcoCycID+Group+Metabolite_class~Measure+Contrast+Stat,value.var = 'Value')

results.castfull.c$G_logFDR_bin<-cut(results.castfull.c$`G_T-C_logFDR`, breaks=logFDRbreaks,labels=logFDRbins)
results.castfull.c$GR_logFDR_bin<-cut(results.castfull.c$`GR_T-C_logFDR`, breaks=logFDRbreaks,labels=logFDRbins)



head(results.m)

results.castfull<-results.m %>%
  filter(Stat %in% c('logFC','SE','FDR','NE','PE','logFDR') & Measure=='G') %>%
  select(-Contrast_type,-Description,-Measure) %>%
  unite(Contrast_Stat,Contrast,Stat) %>%
  spread(Contrast_Stat,Value) %>%
  left_join(data.cast)%>%
  mutate(logFDR_bin=cut(`T-C_logFDR`, breaks=logFDRbreaks,labels=logFDRbins),
         Pole=ifelse(`C_logFC`>0,
                     pole2(`C_logFC`,`T_logFC`),
                     pole2(`C_logFC`,`T_logFC`)+0.5),
         Pole=ifelse(Pole>0.625,Pole-1,Pole),
         Pole360=Pole*360)

head(results.castfull)

  
#left_join(data.sum) %>%
head(results.castfull)
  
#results.castfull$logFDR_bin<-cut(results.castfull$`T-C_logFDR`, breaks=logFDRbreaks,labels=logFDRbins)
# 
# 
# results.castfull$Pole<-ifelse(results.castfull$`C_logFC`>0,
#                               pole2(results.castfull$`C_logFC`,results.castfull$`T_logFC`),
#                               pole2(results.castfull$`C_logFC`,results.castfull$`T_logFC`)+0.5)
# 
# results.castfull$Pole<-ifelse(results.castfull$Pole>0.625,results.castfull$Pole-1,results.castfull$Pole)
# results.castfull$Pole360<-results.castfull$Pole*360
# 


results.cast<-results.all$cast






write.csv(results.exp,paste(odir,'/Ecoli_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/Ecoli_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/Ecoli_results_sidebyside_full.csv',sep=''),row.names = FALSE)







results.cel.r<-read_csv('Celegans/Summary/Celegans_results_withNGM.csv') %>%
  filter(! Index %in% c('Controls-A1','Controls-B1'))


head(results.cel.r)
colnames(results.cel.r)

results.cel<-results.cel.r


colnames(results.cel)[11:30]

colnames(results.cel)[11:30]<-paste('Cel_',colnames(results.cel)[11:30],sep='')

colnames(results.cel)


results.cels<-results.cel[,c('Plate','Well',colnames(results.cel)[11:30])]



head(results.cels)





#Main data
results.castcomb<-merge(results.castfull,results.cels,by=c('Plate','Well'),all=TRUE)



head(results.castcomb)




#With growth rate
results.castcomb.c<-merge(results.castfull.c,results.cels,by=c('Plate','Well'))





metorder<-as.character(results.castcomb[order(results.castcomb$`T-C_logFC`,decreasing = TRUE),'MetaboliteU'])

results.castcomb$MetaboliteU<-factor(results.castcomb$MetaboliteU,levels=metorder,labels = metorder)


head(results.castcomb)
write.csv(results.castcomb,paste(odir,'/Combined_results_sidebyside_full.csv',sep=''),row.names = FALSE)





#Only filtered
metorder<-as.character(results.castcomb[order(results.castcomb$`T-C_logFC`,decreasing = TRUE),'MetaboliteU'])
results.castcomb$MetaboliteU<-factor(results.castcomb$MetaboliteU,levels=metorder,labels = metorder)




#Pick data to show

#Carboxylic acids to remove
coxy<-c('Itaconic Acid','Caproic Acid','Capric Acid','4-Hydroxy Benzoic Acid','2-Hydroxy Benzoic Acid')


selectcast<-subset(results.castcomb,!(Plate %in% c('PM3B','PM4A') & Metabolite %in% PM1PM2$Metabolite ) & !Metabolite %in% c('Negative Control',coxy) )



selectcast.c<-subset(results.castcomb.c,!(Plate %in% c('PM3B','PM4A') & Metabolite %in% PM1PM2$Metabolite ) & !Metabolite %in% c('Negative Control',coxy) )




all.groups<-as.character(unique(data.nc$Group))
strn<-'OP50Sens'
intr<-'T-C'

maincomp<-'C. elegans\nphenotype rescue\n(acs-2 GFP)'


xmin<- -3
xmax<- 4
ymin<- -3
ymax<- 4


lblsize<-4

erralpha<-1
errcolor<-'grey80'
segalpha=0.5

amp<-2
cbrks<-seq(-amp,amp,by=1)
gradcols<-c('blue4','blue','gray80','red','red4')


metorder<-as.character(results.castcomb[order(results.castcomb$`T-C_logFC`),'MetaboliteU'])
#results.castcomb$MetaboliteU<-factor(results.castcomb$MetaboliteU,levels=metorder,labels = metorder)
data.n$MetaboliteU<-factor(data.n$MetaboliteU,levels=metorder,labels=metorder)


ggplot(subset(data.n,Measure=='G'),aes(x=MetaboliteU,y=Value_norm,color=Type))+
  geom_hline(yintercept = 0,color='red',alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point()+
  ylab('Growth logFC vs NGM')+
  xlab('Metabolite')+
  coord_flip()
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ecoli_logFC_tall_raw.pdf",sep=''),
             width=9,height=50)


# 
# metorder.sel<-as.character(selectcast[order(selectcast.c$`GR_T-C_logFC`),'MetaboliteU'])
# selectcast$MetaboliteU<-factor(selectcast$MetaboliteU,levels=metorder.sel,labels = metorder.sel)
# ggplot(selectcast.c,aes(x=MetaboliteU,y=`GR_T-C_logFC`,color=GR_logFDR_bin))+
#   geom_hline(yintercept = 0,color='red',alpha=0.5)+
#   geom_errorbar(aes(ymin=`GR_T-C_NE`,ymax=`GR_T-C_PE`))+
#   geom_point()+
#   xlab('Metabolite')+
#   ylab('Doubling difference (Treatment-Control), log2(OD)/h')+
#   labs(color='Significance (FDR)')+
#   scale_colour_manual(values = c("gray40", "red4", "red3",'red'))+
#   coord_flip()
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Ecoli_logFC_tall_GrowthRate.pdf",sep=''),
#              width=9,height=50)
# #Some will be missing


metorder.sel<-as.character(selectcast[order(selectcast$`T-C_logFC`),'MetaboliteU'])
selectcast$MetaboliteU<-factor(selectcast$MetaboliteU,levels=metorder.sel,labels = metorder.sel)


ggplot(selectcast,aes(x=MetaboliteU,y=`T-C_logFC`,color=logFDR_bin))+
  geom_hline(yintercept = 0,color='red',alpha=0.5)+
  geom_errorbar(aes(ymin=`T-C_NE`,ymax=`T-C_PE`))+
  geom_point()+
  xlab('Metabolite')+
  ylab('Metabolite - metformin interaction as growth logFC vs NGM')+
  labs(color='Significance (FDR)')+
  scale_colour_manual(values = c("gray40", "red4", "red3",'red'))+
  coord_flip()
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ecoli_logFC_tall_Growth.pdf",sep=''),
             width=9,height=50)






#Volcano plots

ggplot(selectcast,aes(x=`T-C_logFC`,y=`T-C_logFDR`,color=-Cel_logFC))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=`T-C_NE`,xmax=`T-C_PE`),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=abs(`Cel_logFC`) ))+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  scale_size(range = c(0.25, 7),name=maincomp)+
  ylab('-logFDR')+
  xlab('Metabolite - metformin interaction as growth logFC vs NGM')+
  labs(color='Significance (FDR)')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Volcano_Ecoli_Growth.pdf",sep=''),
             width=12,height=9)



ggplot(selectcast.c,aes(x=`GR_T-C_logFC`,y=`GR_T-C_logFDR`,color=-Cel_logFC))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=`GR_T-C_NE`,xmax=`GR_T-C_PE`),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=abs(`Cel_logFC`) ))+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  scale_size(range = c(0.25, 7),name=maincomp)+
  ylab('-logFDR')+
  xlab('Doubling rate difference (Treatment-Control), log2(OD)/h')+
  labs(color='Significance (FDR)')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Volcano_Ecoli_GrowthRate.pdf",sep=''),
             width=12,height=9)




# ggplot(selectcast,aes(y=G_T_logFC,x=G_C_logFC,color=Metabolite_class))+
#   geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
#   geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
#   #geom_abline(aes(intercept=1,slope=1),alpha=1,color='black')+
#   geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
#   geom_errorbar(aes(ymin=G_T_NE,ymax=G_T_PE),alpha=erralpha,color=errcolor,width=0)+
#   geom_errorbarh(aes(xmin=G_C_NE,xmax=G_C_PE),alpha=erralpha,color=errcolor,height=0)+
#   ggtitle(paste('Scatterplot of metformin and metabolite supplementation effects ',grp,sep='') )+
#   #,subtitle = paste('Metabolites with FDR<',thres,' are marked',sep='')
#   geom_point(size=5)+
#   # geom_point(aes(size=abs(`Cel_logFC`) ))+
#   coord_cartesian(xlim=c(xmin,xmax),ylim = c(ymin,ymax))+
# 
#   xlab('Growth logFC vs NGM - Control')+
#   ylab('Growth logFC vs NGM - +50mM Metformin')+
#   geom_text_repel(aes(label=ifelse(Metabolite_class=='carbohydrate' , as.character(MetaboliteU),'')),
#                   size=lblsize,nudge_y = 0.3,
#                   force=1,
#                   segment.colour=errcolor,
#                   segment.alpha =segalpha)+
#   # scale_colour_gradientn(colours = gradcols,
#   #                        breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
#   # labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
#   # scale_size(range = c(0.25, 7),name=maincomp)+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Scatter_Control-Treatment_Complete_Growth_Metabolite_class.pdf",sep=''),
#              width=13,height=9, useDingbats=FALSE)





fit<-lm(GR_T_logFC~GR_C_logFC,selectcast.c)
lmeq<-lm_eqn(fit)
a<-fit$coefficients[[2]]
b<-fit$coefficients[[1]]

thres<-0.05

ggplot(selectcast.c,aes(y=GR_T_logFC,x=GR_C_logFC,color=-Cel_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(intercept=b,slope=a),alpha=1,color='red')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=GR_T_NE,ymax=GR_T_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=GR_C_NE,xmax=GR_C_PE),alpha=erralpha,color=errcolor,height=0)+
  ggtitle(paste('Scatterplot of metformin and metabolite supplementation effects ',grp,sep=''),
          subtitle = paste('Metabolites with FDR<',thres,' are marked',sep='') )+
  
  geom_point(size=3)+#aes(size=abs(Cel_logFC))
  #coord_cartesian(xlim=c(xmin,xmax),ylim = c(ymin,ymax))+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  
  xlab('Doubling rate - Control, log2(OD)/h')+
  ylab('Doubling rate - +50mM Metformin, log2(OD)/h')+
  #eval(parse(text = intfdr)) < thres & abs( eval(parse(text = intvar)) ) > 0.75 
  # geom_text_repel(aes(label=ifelse(Cel_FDR<thres, as.character(Metabolite),'')),
  #                 size=lblsize,nudge_y = 0.3,
  #                 force=1,
  #                 segment.colour=errcolor,
  #                 segment.alpha =segalpha)+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp))+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  #scale_size(range = c(0.25, 7),name=maincomp)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Scatter_Control-Treatment_Complete_GrowthRate_Celegans.pdf",sep=''),
             width=13,height=9, useDingbats=FALSE)



#Growth rate and growth comparison

head(selectcast.c)

ggplot(selectcast.c,aes(x=`G_T-C_logFC`,y=`GR_T-C_logFC`,color=-Cel_logFC) )+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  #geom_abline(aes(intercept=b,slope=a),alpha=1,color='red')+
  #geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=`GR_T-C_NE`,ymax=`GR_T-C_PE`),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=`G_T-C_NE`,xmax=`G_T-C_PE`),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=abs(Cel_logFC)))+
  geom_text_repel(aes(label=ifelse(Cel_FDR<thres, as.character(Metabolite),'')),
                  size=lblsize,nudge_y = 0.3,
                  force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)+
  ylab('Doubling rate difference (Treatment-Control), log2(OD)/h')+
  xlab('Metabolite - metformin interaction as growth logFC vs NGM')+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp))+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  scale_size(range = c(0.25, 7),name=maincomp)+
  coord_cartesian(xlim=c(-3,3),ylim = c(-2,2))


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Scatter_GrowthRate_vs_Growth_vs_Celegans.pdf",sep=''),
             width=13,height=9, useDingbats=FALSE)


#Driving effects for interaction
ggplot(selectcast,aes(y=`T-C_logFC`,x=Pole360,color=-Cel_logFC))+
  geom_point()+
  geom_point(aes(size=abs(Cel_logFC)))+
  geom_text_repel(aes(label=ifelse(Cel_FDR<0.05, as.character(Metabolite),'')),
                  size=lblsize,nudge_y = 0.3,
                  force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)+
  ggtitle('Breakdown of driving components for metabolite-metformin interaction')+
  #scale_x_discrete(labels=c("-90"="-T","-45"="-TC","0" = "C", "45" = "TC","90" = "T","135"="T-C","180"="-C","225"="-T-C","270"="-T"))+
  scale_x_continuous(breaks=seq(-135,225,by=45),labels=c("-T-C","-T","-TC","C","TC","T","T-C","-C","-T-C")  )+
  scale_y_continuous(breaks=seq(-4,4,by=1),limits = c(-3,3) )+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp))+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  xlab('Driving component')+
  ylab('Metabolite - metformin interaction as growth logFC vs NGM')+
  scale_size(range = c(0.25, 7),name=maincomp)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Celegans_Ecoli_Driving_effects_for_interaction_CelSignificant.pdf",sep=''),
             width=16,height=12, useDingbats=FALSE)






head(results)




#C elegans, E coli comparison

#theme_set(theme_Publication())


xmin<- -3
xmax<-3
ymin<- -2
ymax<- 2


errcolor<-"grey90"
sbrks<-seq(0,3,by=1)

cbrks<-seq(-1,1,by=1)
gradcols<-c('gray80','blue3','gray80','red3','gray80')

grp<-'Complete'
ggplot(selectcast,aes(x=`T-C_logFC`,
                               y=-Cel_logFC,
                               color=pole(Cel_logFC,`T-C_logFC`)))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  #geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  
  geom_errorbar(aes(ymin=-Cel_PE,ymax=-Cel_NE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=`T-C_NE`,xmax=`T-C_PE`),alpha=erralpha,color=errcolor,height=0)+
  
  geom_point(aes(size=sqrt( Cel_logFC^2 + `T-C_logFC`^2)))+
  #coord_cartesian(xlim=c(xmin,xmax),ylim = c(ymin,ymax))+
  coord_cartesian(xlim=c(xmin,xmax),ylim = c(ymin,ymax))+
  scale_x_continuous(breaks=seq(xmin,xmax,by=1))+
  scale_y_continuous(breaks=seq(ymin,ymax,by=1))+
  
  xlab('Growth rescue in E. coli (metabolite effect normalised), logFC')+
  ylab('Phenotype rescue in C. elegans (acs-2 GFP), logFC')+
  geom_text_repel(aes(label=ifelse(Cel_FDR<0.05 , as.character(MetaboliteU),'')),
                  size=lblsize,nudge_y = 0.3,force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-2,2))+
  scale_size(range = c(0.25, 7),breaks=sbrks,name='Strength of effect')+
  labs(color='Consistency of effect')+
  
  ggtitle(paste('C. elegans and E. coli phenotype rescue in metformin treatment by metabolites: ',grp,sep='') )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#& `G_T-C_FDR`<0.05

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Celegans_Ecoli_effect_correlation_Complete_CelSignificant.pdf",sep=''),
             width=13,height=9, useDingbats=FALSE)



#Venn diagram

EC.ant<-as.character(subset(selectcast,`T-C_FDR`<0.05 & `T-C_logFC`>0 )$MetaboliteU)
EC.syn<-as.character(subset(selectcast,`T-C_FDR`<0.05 & `T-C_logFC`<0 )$MetaboliteU)
EC.non<-as.character(subset(selectcast,`T-C_FDR`>0.05 )$MetaboliteU)
EC.nonstrict<-as.character(subset(selectcast,`T-C_FDR`>0.05 & `C_FDR`>0.05 &`T_FDR`>0.05)$MetaboliteU)
Cel.rescue<-as.character(subset(selectcast,`Cel_FDR`<0.05 & -`Cel_logFC`>0 )$MetaboliteU)
Cel.aggr<-as.character(subset(selectcast,`Cel_FDR`<0.05 & -`Cel_logFC`<0 )$MetaboliteU)

rescue<-list('C. elegans aggravate'=Cel.aggr,
             'E. coli synergistic'=EC.syn,
             'C. elegans rescue'=Cel.rescue,
             'E. coli antagonistic'=EC.ant)

rescue.non<-list('C. elegans aggravate'=Cel.aggr,
             'E. coli non interacting'=EC.non,
             'C. elegans rescue'=Cel.rescue)


rescue.nons<-list('C. elegans aggravate'=Cel.aggr,
                 'E. coli no effect'=EC.nonstrict,
                 'C. elegans rescue'=Cel.rescue)



intersect(EC.syn,Cel.rescue)
setdiff(Cel.rescue,EC.ant)




library('Vennerable')
plot(Venn(rescue),type='ellipses',show = list(Faces = FALSE),doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Celegans_Ecoli.pdf",sep=''),
             width=8,height=8, useDingbats=FALSE)


plot(Venn(rescue.non),show = list(Faces = FALSE),doWeights = TRUE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Celegans_Ecoli_noninteracting.pdf",sep=''),
             width=8,height=8, useDingbats=FALSE)


plot(Venn(rescue.nons),show = list(Faces = FALSE),doWeights = TRUE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Celegans_Ecoli_noeffect.pdf",sep=''),
             width=8,height=8, useDingbats=FALSE)



intersect(Cel.rescue,EC.syn)

intersect(Cel.rescue,EC.ant)


intersect(Cel.aggr,EC.ant)

intersect(Cel.aggr,EC.syn)


intersect(Cel.rescue,EC.non)
intersect(Cel.aggr,EC.non)


intersect(Cel.rescue,EC.nonstrict)
intersect(Cel.aggr,EC.nonstrict)






#Enrichment by class


head(results)

results.ec<-subset(results,Measure=='G' & Contrast=='T-C')[,c('Index','Plate','Metabolite','MetaboliteU','logFC','FDR')]

results.ec$Organism<-'Ec'
head(results.ec)
dim(results.ec)

head(results.cel.r)

results.celegans<-results.cel.r[,c('Index','Plate','Metabolite','MetaboliteU','logFC','FDR')]
results.celegans$logFC<- - results.celegans$logFC
results.celegans$Organism<-'Ce'

head(results.celegans)
dim(results.celegans)

results.enr<-rbind(results.ec,results.celegans)


results.enr.f<-subset(results.enr,!(Plate %in% c('PM3B','PM4A') & Metabolite %in% PM1PM2$Metabolite ) & !Metabolite %in% c('Negative Control',coxy ) )



dim(results.enr)
dim(results.enr.f)
#342
#337

head(results.enr.f)




#Classes preparation
classes.r<-read.xlsx2('../Biolog/Biolog_metabolites_EcoCyc.xlsx',sheetName = 'Classes',header=TRUE,endRow=481)

classes<-merge(info[,c('Index','Plate','Well','Metabolite','MetaboliteU')],classes.r[,c('Index','Class1','Class2')],by='Index')
classes<-subset(classes,Plate %in% c('PM1','PM2A','PM3B','PM4A'))


#View(subset(classes,Metabolite %in% PM1PM2$Metabolite))

classes
dim(classes)

classes.m<-melt(classes,measure.vars = c('Class1','Class2'),variable.name = 'Redundancy',value.name = 'Class')

classes.m<-subset(classes.m, Class !='' &  !(Plate %in% c('PM3B','PM4A') & Metabolite %in% PM1PM2$Metabolite ))
classes.m<-subset(classes.m,!Metabolite %in% c('Negative Control',coxy) )

dim(classes.m)

dim(subset(classes.m,Metabolite %in% PM1PM2$Metabolite & Plate %in% c('PM3B','PM4A')))

subset(classes.m,Well%in% c('A1','F1'))

class.all<-data.frame(table(as.character(classes.m$Class)))
class.all<-rename(class.all,c('Var1'='Class','Freq'='FreqAll') )
sum(class.all$FreqAll)
class.all





classes.c<-classes.m[,c('Index','Class')]

classes.data<-merge(classes.c,results.enr.f,by='Index',all.x=TRUE)


View(classes.data)


int.counts<-ddply(classes.data,.(Class,Organism),summarise,Ant=sum(logFC>0 & FDR <0.05),Syn=sum(logFC<0 & FDR <0.05),Chg=sum(FDR <0.05))
int.m<-melt(int.counts,measure.vars=c('Ant','Syn','Chg'),variable.name = 'Type',value.name = 'Count')

int.m


grp.counts<-ddply(int.m,.(Type,Organism),summarise,FreqGroup=sum(Count,na.rm = TRUE))
grp.counts$Group<-paste(met.counts$Organism,met.counts$Type,sep='_')

grp.counts

# class.counts<-ddply(int.m,.(Organism,Class),summarise,FreqGroup=sum(Count,na.rm = TRUE))
# class.counts



met.cast<-dcast(int.m,Class~Organism+Type,value.var = 'Count')
met.freq.t<-merge(class.all,met.cast,by='Class',all=TRUE)

met.freq.t

met.freq.m<-melt(met.freq.t,id.vars = c('Class','FreqAll'),variable.name = 'Group', value.name = 'Freq')
met.freq.m

met.sum<-merge(met.freq.m,grp.counts[,c('Group','FreqGroup')],by='Group')

met.sum

sum(class.all$FreqAll)
freq.All<-sum(met.freq.t$FreqAll)
freq.All


#phyper(q, m, n, k)
# pop size : 5260
# sample size : 131
# Number of items in the pop that are classified as successes : 1998
# Number of items in the sample that are classified as successes : 62
#phyper(62-1, 1998, 5260-1998, 131, lower.tail=FALSE)

met.sum$p_val<-phyper(met.sum$Freq-1,met.sum$FreqGroup,freq.All-met.sum$FreqGroup,met.sum$FreqAll,lower.tail =FALSE)

#met.sum$FDR<-p.adjust(met.sum$p_val,method = 'fdr')
met.sum$FC<- (met.sum$Freq/met.sum$FreqGroup)/(met.sum$FreqAll/freq.All)

met.sum
met.sum.m<-melt(met.sum,id.vars = c('Class','Group','FreqAll'),measure.vars = c('Freq','p_val','FC'),variable.name = 'Stat',value.name = 'Value')
met.cast<-dcast(met.sum.m,Class+FreqAll~Group+Stat,value.var = 'Value')

met.cast

write.csv(met.cast,paste(odir,'/Metabolite_class_enrichment_nocarboxy.csv',sep=''),row.names = FALSE)


subset(met.cast,Ce_Syn_p_val<0.05)


unique(as.character(selectcast$Metabolite_class))






#EcoCyc enrichment

ecli.co<-read.csv('../Biolog/Ecoli_All_class_instances_non_redundant.csv',header=TRUE,row.names = 1)

#EcoCyc to Biolog links
pmcs<-read.csv('../Biolog/EColi_net/All_Media_Mix_clean_explicit.csv',stringsAsFactors = FALSE)


head(ecli.co)


ecli.co.u<-ecli.co[!duplicated(ecli.co$Class),]
ecli.co.u[,c('Instance','IsClass')]<-NULL




#Enrichment analysis start
encast.r<-selectcast
encast.r$Res_logFC<- -encast.r$Cel_logFC
encast.r$Res_FDR<-encast.r$Cel_FDR

encast<-encast.r[,c('Plate','Well','Index','MetaboliteU','EcoCycID',
                    'T-C_logFC','T-C_FDR',
                    'Res_logFC','Res_FDR')]

colnames(encast)<-gsub('Res','Celegans',colnames(encast))
colnames(encast)<-gsub('T-C','Ecoli',colnames(encast))

ecdata<-merge(encast,pmcs[,c('Plate','Well','Instance')], by=c('Plate','Well'))

dim(ecdata)
head(ecdata)


ecli.co1<-subset(ecli.co,NoInstances>1)
ecli.co2<-subset(ecli.co,NoInstances>2)


ecli.cos<-subset(ecli.co,NoInstances>2 )



cdata<-merge(ecli.cos[,c('Class','NoInstances','subclasses','supclasses','Hierarchy','Instance')],ecdata,by='Instance')

head(cdata)
dim(cdata)
cdata.c<-cdata[!duplicated(cdata[,c('Class','Index','EcoCycID')]),]

#remove redundant assignments to class
dim(cdata.c)

write.csv(cdata.c,paste0(odir,'/EcoCyc_all_classes.csv'))




head(cdata.c)

idvariables<-c('Class','NoInstances','subclasses','supclasses','Hierarchy','Instance','Index','EcoCycID','MetaboliteU')
selstats<-c('logFC','FDR')

cdata.m<-enrichment.melt(cdata,idvariables,selstats)

unique(cdata.m$Contrast)


#remove(enrichmentmelt)

cdata.en<-enrichment(cdata.m,terms = c('Class','NoInstances','subclasses','supclasses','Hierarchy'),IDs = 'Instance',comparisons = 'Contrast',change = 'logFC',sign = 'FDR')


head(cdata.en)
subset(cdata.en,FDR<0.05)


write.csv(cdata.en,paste0(odir,'/EcoCyc_enrichment.csv'))


#Get data ready for heatmap
cdata.cast<-dcast(subset(cdata.en, Term_total>3 & Test!='All'),Class+Hierarchy~Contrast+Test,value.var = 'FDR')

comps<-c('Celegans_Up','Celegans_Down','Ecoli_Up','Ecoli_Down')


rownames(cdata.cast)<-cdata.cast$Class



#131 term
#Leave only enrichment terms with at least one significant value in a row
dim(cdata.cast)
cdata.cast<-cdata.cast[apply(cdata.cast[,comps], 1, function(x) { any(as.numeric(x) < 0.05,na.rm=TRUE)}),]
dim(cdata.cast)


cdata.cast[order(cdata.cast$Hierarchy),]



cdata.cast$Code<-apply(cdata.cast[,comps],1,function(x) paste0(as.numeric(x< 0.05),collapse = '') )


#write.csv(cdata.cast,paste0(odir,'/EcoCyc_enrichment_heatmap.csv'))



EChm<-read.xlsx(paste0(odir,'/EcoCyc_enrichment_heatmap.xlsx'),sheetName = 'Heatmap',header = TRUE)



cdata.cast.f<-subset(cdata.cast,Class %in% EChm[is.na(EChm$Remove),'Class'])


cdata.cast.f<-cdata.cast.f[order(cdata.cast.f$Class),]





reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
gyrs<-colorRampPalette(c("gray90","steelblue1","blue4"))(n = 5)
enbrks<-c(0,-log(0.05,10),2,3,4,5)



hdata<-cdata.cast.f[,comps]


rownames(hdata)<-gsub('\\|','',rownames(hdata))


hdata<- -log10(hdata)

hdatafill<-hdata
hdatafill[is.na(hdatafill)]<-0
hdatafill[hdatafill==Inf]<-30

#Use the dummy dataset to get ordering of rows (enrichment terms)
#Graphical settings here are irrelevent for the final plot
#Data needs to be supplied as a matrix - table with row/column names defined implicitly and only numeric values
hmap<-heatmap.2(data.matrix(hdatafill),
                key=TRUE,
                Colv=FALSE,
                trace='none',
                col=gyrs,
                xlab='Comparison',
                Rowv=TRUE,
                dendrogram="row",
                scale="none",
                na.color="white", 
                symkey=FALSE,
                breaks=enbrks,
                reorderfun=reorderfun_mean,
                cexRow=0.8,
                cexCol=0.5)

#Use ordering with the real data
heatmap.2(data.matrix(hdata),
          key=FALSE, #Colour key needs to be provided separately
          Colv=FALSE, #Columns should not be reordered
          Rowv=hmap$rowDendrogram, # Here comes in the ordering
          trace='none',
          col=gyrs, # Colour scale
          breaks=enbrks, # Colour breaks
          xlab='Comparison',
          dendrogram='none', #Row dendogram, but should be changed to none, as dendrogram represents our data with filled-in values
          scale="none", #Should values be normalised in rows or columns - No
          na.color="white", # What colour to use with not missing values
          symkey=FALSE, #Provided colour scale is not symetrical
          cexRow=0.7, #Some figure scaling parameters. Works only after a lot of experimentation
          cexCol=0.7,
          margin=c(10,20),
          lwid=c(0.2,0.8),
          lhei=c(0.05,0.95))

dev.copy2pdf(device=cairo_pdf,file=paste0(odir,'/Enrichment_heatmap_EcoCyc_UpDown_clean.pdf'),
             width=5,height=5,useDingbats=FALSE)



#Enrichment tables
fcodes<-unique(cdata.cast.f$Code)




all.en.all$p <- format(all.en.all$p, scientific = TRUE,digits = 2)
all.en.all$q <- format(all.en.all$q, scientific = TRUE,digits = 2)



grid.table(cdata.cast.f)






#Chebi enrichment

# source("https://bioconductor.org/biocLite.R")
# biocLite('ontoCAT')
# 
# install.packages('ontologyIndex')
# install.packages('ontoCAT')
# 
# library(ontologyIndex)
# 
# library(ontoCAT)
# efo <- getEFO()
# 
# chebi<-getOntology('../ChEbi/chebi.obo')



pmcs.c<-subset(pmcs,ChEBI_ID!='')[,c('Plate','Well','Instance','ChEBI_ID')]
pmcs.c$fChEBI_ID<-paste0('CHEBI:',pmcs.c$ChEBI_ID)

head(pmcs.c)

pmcs.counts<-ddply(pmcs.c,.(Plate,Well),summarise,Count=length(ChEBI_ID))

pmcs.co<-merge(pmcs.c,pmcs.counts,by=c('Plate','Well'))

head(pmcs.co)


head(encast)

chedata<-merge(pmcs.co,encast,by=c('Plate','Well'))


head(chedata)


chedata$Ecoli_All<-ifelse(chedata$Ecoli_FDR<0.05,1/chedata$Count,0)
chedata$Ecoli_Up<-ifelse(chedata$Ecoli_FDR<0.05 & chedata$Ecoli_logFC>0, 1/chedata$Count, 0)
chedata$Ecoli_Down<-ifelse(chedata$Ecoli_FDR<0.05 & chedata$Ecoli_logFC<0, 1/chedata$Count, 0)


chedata$Celegans_All<-ifelse(chedata$Celegans_FDR<0.05,1/chedata$Count,0)
chedata$Celegans_Up<-ifelse(chedata$Celegans_FDR<0.05 & chedata$Celegans_logFC>0, 1/chedata$Count, 0)
chedata$Celegans_Down<-ifelse(chedata$Celegans_FDR<0.05 & chedata$Celegans_logFC<0, 1/chedata$Count, 0)



head(chedata)

View(chedata)

write.csv(chedata,paste0(odir,'/ChEBI_annotation.csv'))




chebi<-data.frame()

chfls<-list.files('Summary/ChEBI_enrichment/')

for (fl in chfls) {
  print(fl)
  ed<-gsub('.txt','',fl,fixed = TRUE)
  ed<-gsub('ChEBI_','',ed)
  ce<-read.table(paste0('Summary/ChEBI_enrichment/',fl),sep='\t',header = TRUE)
  ce$Contrast<-strsplit(ed,'_')[[1]][1]
  ce$Type<-strsplit(ed,'_')[[1]][2]
  chebi<-rbind(chebi,ce)
}


View(chebi)

chebi.cast<-dcast(chebi,ChEBI_Name~Contrast+Type,value.var='Corr.PValue')







#Enrichment


dir.create(paste(odir,"/Biolog_enrichment/",sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")

enrt<-subset(data.nc,Group %in% c('C-Source','N-Source','P-Source','S-Source') & 
               !Name %in% c('Negative Control',"D,L-Malic Acid","beta-Methyl-D-Galactoside","O-Phospho-D-Tyrosine") &
               !KEGG_ID %in% c("") &
               !( Plate %in% c('PM3B','PM4A') & Name %in% PM1PM2$Name ) &
               !( Plate %in% c('PM3B','PM4A') & KEGG_ID %in% PM1PM2$KEGG_ID ) & 
               !( Plate=='PM4A' & KEGG_ID=='C00097' ) )[,c('SampleID','Sample','Plate','Index','Group','Name','KEGG_ID','Measure','Value_norm')]
head(enrt)

enrc<-dcast(enrt,SampleID+Sample+Measure~KEGG_ID,value.var = 'Value_norm')

head(enrc)

KEGGS<-as.character(unique(enrc$KEGG_ID))

#Find duplicates
#Works when aggregate function is length
dKEGGS<-KEGGS[apply(enrc[,KEGGS],2, function(x) any(x>1))]
dKEGGS
#subset(enrs,KEGG_ID %in% dKEGGS[2])
length(unique(enrt$KEGG_ID))


KEGGbac<-data.frame(unique(enrt$KEGG_ID))

write.csv(KEGGbac,paste(odir,"/Biolog_enrichment/Background.csv",sep=''),row.names = FALSE,col.names = NA,quote=FALSE)


for (gr in c('C-Source','N-Source','P-Source','S-Source','Complete')){
  print(gr)
  if (gr=='Complete') {
    enrs<-subset(enrt,Measure=='G')
  } else {
    enrs<-subset(enrt,Group==gr & Measure=='G')
  }
  enrc<-dcast(enrs,SampleID+Sample~KEGG_ID,value.var = 'Value_norm')
  write.csv(enrc,paste(odir,"/Biolog_enrichment/Biolog_Enrichment_",gr,"_T-C_",strn,".csv",sep=''),row.names = FALSE)
}





#KEGG plotting


library(org.EcK12.eg.db)
library(pathview)

#allpaths<-c('01110')




path2pathde<-getKEGGPathwayNames('eco',remove.qualifier = TRUE)
path2pathde$PathwayID<-gsub('path:','',path2pathde$PathwayID)
allpaths<-gsub('eco','',path2pathde$PathwayID)


write.csv(path2pathde,paste(odir,'/KEGG/KEGG_pathways.csv',sep=''))


allqea<-read.csv(paste(odir,'/Biolog_enrichment/pathway_results.csv',sep=''),sep=',')

allqea<-rename(allqea,c('X'='KEGG pathway'))
allqea$logFDR<- -log10(as.numeric(allqea$FDR))
allqea$Impact<-as.numeric(allqea$Impact)
allqea$Hits<-as.numeric(allqea$Hits)
allqea$`Total.Cmpd`<-as.numeric(allqea$`Total.Cmpd`)
allqea$Ratio<-allqea$Hits/allqea$`Total.Cmpd`

allqea.a<-merge(path2pathde,allqea,by.x='Description',by.y='KEGG pathway',all.y=TRUE)
subset(allqea.a,is.na(PathwayID))

head(allqea.a)

thres<-0.05
ggplot(allqea.a,aes(x=Impact,y=logFDR,size=Ratio,color=Hits/`Total.Cmpd`))+
  geom_point()+
  ylab('-log(FDR)')+
  xlab('Pathway impact')+
  geom_text_repel(aes(label=ifelse(FDR < thres & Impact > 0.25 , as.character(`Description`),'') ),
                  size=2,nudge_y = 0.3,force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)

# cnm<-c('Treatment_Sens_C','Treatment_Res_C','Treatment_Res-a_C','Treatment_Res-p_C',
#        'Treatment_Sens_N','Treatment_Res_N','Treatment_Res-a_N','Treatment_Res-p_N')
# #       'Strain_Res_C','Strain_Res-a_C','Strain_Res-p_C',
# #       'Strain_Res_N','Strain_Res-a_N','Strain_Res-p_N'
# 
# qeacast<-dcast(allqea,`Metabolite Set`~Type+Strain+Element,value.var = 'logFDR',drop = FALSE,fill = NA)
# rownames(qeacast)<-qeacast$`Metabolite Set`
# 
# 
# #Filter rows with no significance
# qeac<-qeacast[!apply(qeacast[,cnm], 1, function(x)(all(x< -log(0.05,10),na.rm=TRUE))),][,cnm]
# 
# #Create Filler
# hdatafill<-qeac
# hdatafill[is.na(hdatafill)]<-0
# 
# 
# gyrs<-colorRampPalette(c("gray90","steelblue1","blue4"))(n = 4)
# enbrks<-c(0,-log(0.05,10),2,3,4)
# 
# hmap<-heatmap.2(as.matrix(hdatafill),key=TRUE,Colv=FALSE,trace='none',
#                 col=gyrs,xlab='Comparison',Rowv=TRUE,
#                 dendrogram="row",scale="none",na.color="white",
#                 cexRow=0.8,cexCol=0.5,margin=c(8,16),
#                 lwid=c(0.2,0.8),symkey=FALSE,
#                 reorderfun=reorderfun_mean,breaks=enbrks)
# 
# 
# 
# #cairo_pdf(paste(odir,'/Metabolomics_',entype,'_heatmap.pdf',sep = ''),width=8,height=15)
# hmapr<-heatmap.2(as.matrix(qeac),key=TRUE,Colv=FALSE,trace='none',col=gyrs,
#                  xlab='Comparison',Rowv=hmap$rowDendrogram,
#                  dendrogram='row',scale="none",na.color="white",
#                  cexRow=0.7,cexCol=0.7,margin=c(16,16),
#                  lwid=c(0.2,0.8),symkey=FALSE,
#                  reorderfun=reorderfun_mean,
#                  breaks=enbrks)
# #dev.off()
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Biolog_Enrichment_CN.pdf',sep = ''),
#              width=6,height=17)
















#Time series
timedata<-read.table('Data/Data.csv',sep=',',quote = '"',header = TRUE)
timedata<-subset(timedata,Data=='750nm_f')
colnames(timedata)


colnames(timedata)<-gsub('X','',colnames(timedata))
colnames(timedata)

timem<-melt(timedata,id=colnames(timedata)[1:23],variable.name='Time_sec',value.name='OD')
timem$Time_sec<-as.numeric(as.character(timem$Time_sec))
timem$Time_min<-timem$Time_sec/60
timem$Time_h<-timem$Time_min/60

timem$Name<-as.factor(timem$Name)
timem$Replicate<-as.factor(timem$Replicate)
timem$Metformin_mM<-as.factor(timem$Metformin_mM)




tssum<-ddply(timem,.(Plate,Well,Index,Metformin_mM,Time_h,Name),summarise,Mean=mean(OD),SD=sd(OD))
tssum$Name<-relevel(tssum$Name,ref='Negative Control')

nindx<-c('PM1-A1',selindx)

ggplot(subset(tssum,Index %in% nindx),aes(x=Time_h,fill=Metformin_mM))+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD))+
  geom_line(aes(y=Mean),colour="black")+
  scale_x_continuous(breaks=seq(0,24,by=6))+
  ylab("OD")+
  xlab("Time, h")+
  labs(fill="Metformin, mM")+
  facet_grid(~Name)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ecoli_growth_Stat_Selection.pdf",sep=''),
             width=16,height=4, useDingbats=FALSE)






ggplot(subset(timem,Index %in% selindx),aes(x=Time_h,y=OD,color=Metformin_mM))+
  geom_line(aes(group=interaction(Name,Replicate,Metformin_mM)))+
  scale_x_continuous(breaks=seq(0,24,by=6))+
  ylab("OD")+
  xlab("Time, h")+
  labs(color='Metformin, mM')+
  facet_grid(.~Name)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ecoli_growth_replicates_Selection.pdf",sep=''),
             width=16,height=4, useDingbats=FALSE)








#KEGG mapping





#Only unique values with preference for PM1 PM2A
all.results.rcp<-subset(results.castfull,
                        !Metabolite %in% c('Negative Control',"D,L-Malic Acid","beta-Methyl-D-Galactoside","O-Phospho-D-Tyrosine") &
                          !KEGG_ID %in% c("") &
                          !is.na(KEGG_ID) &
                          !( Plate %in% c('PM3B','PM4A') & Metabolite %in% PM1PM2$Name ) &
                          !( Plate %in% c('PM3B','PM4A') & KEGG_ID %in% PM1PM2$KEGG_ID ) & 
                          !( Plate=='PM4A' & KEGG_ID=='C00097'))

keggxml<-'~/Dropbox/Projects/2015-Metformin/Annotations/Ecoli/KEGG_pathways'

pathcomp<-c('OP50Sens_T-OP50Sens_C_logFC')


mapkey<-'KEGG_ID'
nrow(all.results.rcp)
gdata<-all.results.rcp
#gdata<-all.results.rcp
nrow(gdata)

dupl<-duplicated(gdata[,mapkey])
print('KEGG ID duplicates')
print(table(dupl))
gdataf<-gdata[!dupl,]
rownames(gdataf)<-gdataf[,mapkey]

print(paste('Total KEGG pathways to plot:',length(allpaths)))

keggdir<-paste(cwd,odir,'/KEGG',sep='')
dir.create(keggdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
setwd(keggdir)
#comp.gr<-'Main'

limitcpd<-3

# pathcomp<-comp.list[[comp.gr]]
# outsuffx<-comp.gr
pv.out <- pathview(cpd.data = gdataf[,pathcomp,drop=FALSE],
                   pathway.id = allpaths,
                   cpd.idtype = "kegg",
                   species = "eco",
                   out.suffix = 'Complete',
                   kegg.dir = keggxml,
                   same.layer=FALSE,kegg.native = T,map.symbol=FALSE,
                   map.null=FALSE,new.signature=FALSE, plot.col.key=TRUE,
                   res=300,
                   min.nnodes = 0,
                   limit=list(gene=2,cpd=limitcpd),node.sum = "mean",
                   low = list(gene = "blue", cpd = "blue"),
                   mid=list(gene = "gray", cpd = "gray"),
                   high = list(gene = "red", cpd ="red"))
setwd(cwd)

all.kegg.mappingst<-data.frame()
paths<-names(pv.out)
for (pth in paths) {
  result<-pv.out[[pth]][['plot.data.cpd']]
  result$PathwayID<-pth
  all.kegg.mappingst<-rbind(all.kegg.mappingst,result)
}
all.kegg.mappingst$mol.col<-NULL

all.kegg.mappings<-merge(path2pathde,all.kegg.mappingst,by='PathwayID',all.y=TRUE)


write.csv(all.kegg.mappings,paste(odir,'/All_KEGG_mappings_Complete.csv',sep = ''),row.names = FALSE)




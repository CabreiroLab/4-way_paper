library(gplots)
library(ggplot2)
library(xlsx)
library(plyr)
library(reshape2)
library(ggbiplot)
library(ggrepel)
library(car)
library(heatmap3)

library(RColorBrewer)
library(grid)
library(gridExtra)
library(plot3D)
library(rgl)
library(pca3d)
library(multcomp)
library(contrast)

library(gtools)

#remove(enrichment.melt,hypothesise,mymerge)

devtools::install_github("PNorvaisas/PFun")
library(PFun)


lm_eqn = function(m) {
  fres<-summary(m)
  l <- list(a = format(abs(coef(m)[2]), digits = 2),
            b = format(coef(m)[1], digits = 2),
            r2 = format(summary(m)$r.squared, digits = 2),
            p2 = format(pf(fres$fstatistic[1], fres$fstatistic[2], fres$fstatistic[3],lower.tail = FALSE)[[1]], digits = 2,scientific=TRUE));
  
  if (coef(m)[2] >= 0)  {
    cof <- substitute(italic(y) == b + a %.% italic(x),l)
    full <- substitute(italic(y) == b + a %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  } else {
    cof <- substitute(italic(y) == b - a %.% italic(x),l) 
    full <- substitute(italic(y) == b - a %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  }
  
  stat<-substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  return(list('Coef'=as.character(as.expression(cof)),
              'Stat'=as.character(as.expression(stat)),
              'Full'=as.character(as.expression(full)),
              'Atop'=as.character(as.expression(paste(cof,'\n',stat,sep='') ))))                 
}

reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }

theme_set(theme_light())


setwd("~/Dropbox/Projects/2015-Metformin/GFP_reporters")


odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


rawdata<- read.csv("UAL_screen/UAL_all_data.csv",header = TRUE)


#Change names to avoid confusion
rawdata$Type<-recode(rawdata$Type,c("'NC'='PC';'NM'='PM';'NGC'='PGC';'NGM'='PGM'"))

head(rawdata)

#Explanation of sample names
codes<-c('PC','PM','PGC','PGM')
sampledescriptions<-c('NGM',
                      'NGM + 100mM Metformin',
                      'NGM + 0.2% Glucose',
                      'NGM + 0.2% Glucose + 100mM Metformin')


explanation<-data.frame(Code=codes,Sample_description=sampledescriptions)





#TFs

TF<-read.csv('../Annotations/Ecoli/Transcription_Factors/network_tf_gene.txt',comment.char = '#',sep = '\t',header=FALSE)
colnames(TF)<-c('TF','Gene','Effect','Evidence','Support','NAs')
TF$NAs<-NULL



#Remove duplicated effect descriptions
TFu<-TF[!duplicated(TF[,c('TF','Gene')]),]
dim(TFu)


#Annotation

annotation<-read.csv('../Annotations/Ecoli/UAL/UAL_reannotated.csv')


annotation<-rename(annotation,c("UAL_promoter"="Promoter"))
annotation$Strand<-ifelse(is.na(annotation$Strand),'1',annotation$Strand)
annotation.u<-annotation[!duplicated(annotation$Promoter),]
annotation.u$Plate<-NULL
annotation.u$Well<-NULL



#Data
data<-dcast(rawdata,Type+Replicate+Plate+Well~Reading,
                  value.var = "Value")

data$Group<-paste(data$Type,data$Replicate,sep='_')
data$Index<-paste(data$Plate,data$Well,sep='_')



data$GFPOD<-data$GFP/data$OD
data$logGFPOD<-log2(data$GFPOD)
data$logGFPOD<-ifelse(data$logGFPOD<0,0,data$logGFPOD)


data$Type<-factor(data$Type,levels=c('PC','PM','PGC','PGM'),labels=c('PC','PM','PGC','PGM'))

data.a<-merge(data,explanation,by.x='Type',by.y='Code',all.x=TRUE)
data.an<-merge(data.a,annotation,by=c('Plate','Well'),all.x=TRUE)

head(data)

data.an$Promoter<-as.factor(data.an$Promoter)
data.an$Promoter<-relevel(data.an$Promoter, ref = "U139")


head(data.an)

unique(data.an$Plate)


proms.t<-subset(data.an,!Promoter %in% c('Empty',NA))
proms.t$Plate<-factor(proms.t$Plate,levels=as.character(1:21),labels=as.character(1:21))
proms.t$Replicate<-factor(proms.t$Replicate,levels=as.character(1:4),labels=as.character(1:4))


ref.t<-subset(proms.t,Promoter %in% c('U139','U66') )
ref.t<-rename(ref.t,c('logGFPOD'='Ref_RP'))


subset(ref.t,Plate==20)

#Plate 21 does not have internal control
ref.21<-subset(ref.t,Plate==20)
ref.21$Plate<-21




Ref<-rbind(ref.t,ref.21)


ggplot(Ref,aes(x=Plate,y=Ref_RP,color=Type))+
  geom_boxplot()+
  ggtitle('Promoterless plasmid fluorescence by plate')+
  facet_wrap(~Promoter)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Boxplot_Promoterless_by_plate.pdf",sep=''),
             width=9,height=4, useDingbats=FALSE)

ggplot(Ref,aes(x=Replicate,y=Ref_RP,color=Type))+
  geom_boxplot()+
  ggtitle('Promoterless plasmid fluorescence by replicate')+
  facet_wrap(~Promoter)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Boxplot_Promoterless_by_replicate.pdf",sep=''),
             width=9,height=4, useDingbats=FALSE)

ggplot(Ref,aes(x=Type,y=Ref_RP,color=Type))+
  geom_boxplot()+
  ggtitle('Promoterless plasmid fluorescence by type')+
  facet_wrap(~Promoter)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Boxplot_Promoterless_by_type.pdf",sep=''),
             width=9,height=4, useDingbats=FALSE)




proms<-merge(proms.t,Ref[,c('Strand','Type','Plate','Replicate','Ref_RP')],
             by=c('Type','Plate','Replicate','Strand'),all.x = TRUE)

dim(proms)

head(proms)

proms$logGFPOD_norm<-proms$logGFPOD-proms$Ref_RP
proms$Promoter<-as.factor(proms$Promoter)

proms$Drug<-as.factor(ifelse(proms$Type %in% c('PC','PGC'),'Control','Metformin'))
proms$Sugar<-as.factor(ifelse(proms$Type %in% c('PC','PM'),'None','Glucose'))
proms$Drug<-factor(proms$Drug,levels=c("Control","Metformin"),labels=c("Control","Metformin"))
proms$Sugar<-factor(proms$Sugar,levels=c("None","Glucose"),labels=c("None","Glucose"))


#Investigate growth deviations

ggplot(proms,aes(x=OD,fill=Type,color=Type))+
  geom_histogram(position='identity',alpha=0.5)+
  facet_grid(Replicate~Type,labeller = label_both)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Distribution_OD.pdf",sep=''),
             width=9,height=6)


discarded<-subset(proms, OD<0.25 )
discarded


discarded.p<-as.character(unique(discarded$Promoter))
discarded.p




#Leave only growing
proms.c<-subset(proms,OD>0.25 | Promoter %in% c('U139','U66') )
proms.c$Promoter<-factor(as.character(proms.c$Promoter))



proms.c$PReplicates<-makereplicates(proms.c[,c('Promoter','Type','Replicate')])

table(subset(proms.c,PReplicates==1)$Promoter)
table(proms.c$PReplicates)



proms.m<-melt(proms.c,id.vars = c('Type','Plate','Well','Strand','Promoter','Replicate','Group','Index'),
              measure.vars = c('GFP','OD','GFPOD','logGFPOD','logGFPOD_norm'),variable.name = 'Measure',value.name = 'Value')


stats<-ddply(proms.m,.(Promoter,Type,Measure),summarise,
             Mean=mean(Value),SD=sd(Value))

stats$Type<-factor(stats$Type,levels=c('PC','PM','PGC','PGM'),labels=c('PC','PM','PGC','PGM'))
stats.m<-melt(stats,id.vars = c('Type','Promoter','Measure'),variable.name = 'Stat',value.name = 'Stat_Value')

table(stats.m$Measure)

stats.cast<-dcast(subset(stats.m,Measure=='logGFPOD_norm'),Promoter~Type+Stat,value.var = 'Stat_Value')



head(proms.c)

proms.ms<-melt(proms.c,id.vars = c('Type','Plate','Well','Promoter','Replicate','Group','Index','Strand',
                                   'Reports_genes','Gene','Operon','Operon_name','Duplicated_reporting','OD'),
               measure.vars = c('logGFPOD','logGFPOD_norm'),variable.name = 'Measure',value.name = 'Value')



ggplot(proms.ms,aes(y=OD,x=Value,color=Replicate))+
  ylim(0,2)+
  xlab('logGFPOD')+
  geom_point(size=1)+
  scale_x_continuous(breaks=seq(-4,20,by=2),limits=c(-2,16))+
  ggtitle('Raw data before and after normalisation')+
  geom_smooth(method='lm')+
  facet_grid(Measure~Type)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/GFP_OD_relationship.pdf",sep=''),
             width=9,height=4)



ggplot(proms.c,aes(x=logGFPOD_norm,fill=Replicate,color=Replicate))+
  geom_vline(xintercept = 0,color='blue')+
  geom_histogram(position='identity',alpha=0.5)+
  facet_grid(Replicate~Type,labeller = label_both)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Distribution_in_replicates.pdf",sep=''),
             width=9,height=6)




fh<-subset(stats, Measure %in% c('logGFPOD','logGFPOD_norm'))
U139.s<-subset(fh, Promoter %in% c('U139','U66'))


ggplot(fh,aes(x=Mean,fill=Type))+
  geom_histogram(position='identity',alpha=0.5)+
  xlab('logGFPOD')+
  geom_vline(data=U139.s,aes(xintercept=Mean,color=Type))+
  geom_text(data=U139.s,aes(label=Type,color=Type,x=Mean+0.5,y=(as.numeric(Type)+1)*20+750),size=2)+
  ggtitle('Promoter activity distribution before and after promoterless plasmid normalisation')+
  facet_grid(.~Measure)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Distribution_raw_normalised.pdf',sep = ''),
             width=6,height=4,useDingbats=FALSE)




testproms<-c('U66','U139','lacZ','crp','serA','wrbA','gadW','gadA','gadX','gadB','ompA','rplN')
proms.sel<-subset(proms.ms,Promoter %in% testproms)
proms.sel$Promoter<-factor(proms.sel$Promoter,levels=testproms,labels=testproms)
#proms.sel$Promoter<-relevel(proms.sel$Promoter,ref='U139')

ggplot(proms.sel,aes(y=Value,x=Promoter,color=Type))+
  geom_boxplot()+
  scale_y_continuous(breaks=seq(0,14,by=1))+
  ylab('logGFPOD')+
  facet_grid(.~Measure,labeller=label_both)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Boxplot_Selected_comparison.pdf",sep=''),
             width=12,height=6, useDingbats=FALSE)




#Duplicated reporting comparion
dupl.data<-subset(proms.ms, Measure=='logGFPOD_norm' & !is.na(Duplicated_reporting ))

dupl.data.m<-as.data.frame(splitstackshape::cSplit(dupl.data,splitCols = 'Duplicated_reporting',sep=',',direction='long',drop = FALSE))
dupl.data.m$IsOperon<-!is.na(dupl.data.m$Operon_name)

gn<-subset(dupl.data.m,!IsOperon)[,c('Promoter','Type','Index','Replicate','Duplicated_reporting','Value')]
op<-subset(dupl.data.m,IsOperon)[,c('Type','Index','Replicate','Duplicated_reporting','Value')]


dupl.gn.promoters<-as.character(unique(gn$Promoter))

dr.sbs<-merge(gn,op,by=c('Type','Replicate','Duplicated_reporting'),all=TRUE,suffixes=c("_gene","_operon"))


head(dupl.data.m)


ggplot(dupl.data.m,aes(y=Value,x=Type,color=IsOperon))+
  geom_hline(yintercept=0,color='red',alpha=0.5)+
  geom_boxplot(aes(group=interaction(Type,Promoter,IsOperon)))+
  scale_y_continuous(breaks=seq(-5,14,by=1))+
  ylab('logGFPOD_norm')+
  facet_wrap(~Duplicated_reporting,ncol = 5)#,labeller=label_both


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Duplicated_reporting.pdf",sep=''),
             width=10,height=50, useDingbats=FALSE)


ggplot(dr.sbs,aes(x=Value_gene,y=Value_operon))+
  geom_hline(yintercept=0,color='red',alpha=0.5)+
  geom_vline(xintercept=0,color='red',alpha=0.5)+
  geom_point()+
  facet_grid(Replicate~Type)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Operon_gene_activation_comparison.pdf",sep=''),
             width=10,height=10, useDingbats=FALSE)





# overview<-data.frame(type=c('raw','normalised'),varname=c('logGFPOD','logGFPOD_norm'))
# rid<-2
# 
# #for(rid in 1:nrow(overview)){
#   
# type<-as.character(overview[rid,'type'])
# varname<-as.character(overview[rid,'varname'])

type<-'normalised'
varname<-'logGFPOD_norm'
  
print(type)
pcashape<-dcast(proms,Type+Group+Sample_description~Index,value.var = varname)
rownames(pcashape)<-pcashape$Group
clean_data<-pcashape[,-c(1:3)]
samples<-pcashape$Group
Type<-pcashape$Type
sampleclasses<-pcashape$Sample_description

clean_data<-clean_data[,apply(clean_data,2,function(x) all(!is.na(x)) & !all(x==0)) ]


hca_sample<-hclust(dist(clean_data,method="manhattan"),method="complete")
#hca_variable<-hclust(dist(t(log_data),method="manhattan"),method="complete")
fname<-paste(odir,"/Hierarchical_Clustering_",type,".pdf",sep='')
print(fname)
cairo_pdf(fname,width=9,height=6)
plot(hca_sample, labels=samples,
     hang=-1, main=paste("Cluster Dendrogram - ",type,sep=''), xlab="", sub="", cex=1)
#print(dendro)
dev.off()

# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Hierarchical_Clustering_",type,".pdf",sep=''),
#              width=9,height=6)

ir.pca <- prcomp(clean_data,
                 center = TRUE,
                 scale. = TRUE)

write.csv(ir.pca[2],paste(odir,"/PCA_loadings_",type,".csv",sep=''))

fname<-paste(odir,"/PCA_components_",type,".pdf",sep='')
print(fname)
cairo_pdf(fname,width=9,height=6)
plot(ir.pca,type='l')
dev.off()


groupCluster <- kmeans(clean_data, 4, nstart = 100)
#groupCluster


clusters<-as.factor(groupCluster$cluster)


generalpca <- ggbiplot(ir.pca, obs.scale = 1,
                       var.scale = 1,
                       #groups= clusters,
                       groups = sampleclasses,
                       ellipse = TRUE,
                       circle = TRUE,
                       var.axes = 0)+
  #   geom_vline(aes(xintercept=0),color='gray80')+
  #   geom_hline(aes(yintercept=0),color='gray80')+
  scale_color_discrete(name = '')+
  ggtitle('PCA')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


fname<-paste(odir,"/PCA_",type,".pdf",sep='')
print(fname)
cairo_pdf(fname,width=9,height=6)
print(generalpca)
dev.off()




heatshape<-dcast(proms,Index~Group,value.var = varname)
rownames(heatshape)<-heatshape$Index
heatshape$Index<-NULL
heatshape<-heatshape[apply(heatshape,1,function(x) all(!is.na(x)) & !all(x==0)) ,]

#Remove numbers from sample IDs
groups<-gsub('_[[:digit:]]+', '', colnames(heatshape))
annData<-factor(groups)
#ColSideColors<-rainbow(length(unique(annData)))[as.numeric(annData)]

ColSideColors<-rainbow(length(unique(annData)))[as.numeric(annData)]
legtxt<-as.character(unique(annData))
legcol<-rainbow(length(unique(annData)))[as.numeric(unique(annData))]

fname<-paste(odir,"/Heatmap_",type,".pdf",sep='')
print(fname)
cairo_pdf(fname,width=6,height=9)

hmap<-heatmap3(as.matrix(heatshape),
         scale = 'row',
         labRow=FALSE,
         ColSideLabs = 'Group',
         ColSideColors = ColSideColors,
         legendfun=function() showLegend(legend=legtxt,col=legcol,cex=1))
#hmap
#print(hmap)
dev.off()

#}


#install.packages('kernlab')
library(kernlab)

clean_data.x<-as.matrix(clean_data)

PCdf<-data.frame(ir.pca$x)


groupCluster <- kmeans(clean_data, 4, nstart = 100)
clusters<-as.factor(groupCluster$cluster)




groupSpecGauss<-specc(clean_data.x,centers=4, kernel = "rbfdot",iterations=1000)#,iterations=100
clusters<-as.factor(groupSpecGauss)


PCdf$clusters<-as.factor(groupSpecGauss)

ggbiplot(ir.pca, obs.scale = 1,
                       var.scale = 1,
                       groups= clusters,
                       #groups = sampleclasses,
                       ellipse = TRUE,
                       circle = TRUE,
                       var.axes = 0)+
  #   geom_vline(aes(xintercept=0),color='gray80')+
  #   geom_hline(aes(yintercept=0),color='gray80')+
  scale_color_discrete(name = '')+
  ggtitle('PCA')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#Solve this to have more flexibility
ggplot(PCdf,aes(x=PC1,y=PC2,color=clusters))+
  geom_point()+
  stat_ellipse(aes(x=PC1,y=PC2,color=clusters))





#library(rgl)
#library(pca3d)

#sampleclasses


#Needs to be done manually

fname<-paste(odir,"/PCA3D_",type,".pdf",sep='')
#cairo_pdf(fname,width=9,height=9)
pca3d(ir.pca,components=1:3,group=clusters,
      radius=2,
      #fancy=TRUE,
      show.centroids = TRUE,
      #show.group.labels = TRUE,
      show.shadows = TRUE,
      show.shapes = TRUE,
      show.plane = TRUE)
#dev.off()
#snapshotPCA3d()
#rgl.postscript(fname,"pdf")





#Linear modelling starts here



as.character(unique(dupl.data$Promoter))
refproms<-c('U139','U66')


proms.lm<-subset(proms.c,!Promoter %in% refproms)
dim(proms.lm)

#Remove duplicated reporting
proms.lm<-subset(proms.lm,! (!is.na(Duplicated_reporting) & is.na(Operon_name) )  )
dim(proms.lm)


#Use data with deviant growth removed
#Some promoters are present in multiple wells, nead unique ID
lmshape<-dcast(proms.lm,Group+Type+Drug+Sugar+PReplicates~Promoter,value.var = 'logGFPOD_norm',drop=TRUE)


#All promoters
promoters<-as.character(unique(proms.lm$Promoter))
length(promoters)
length(promoters)

promoters


setdiff(annotation.u$Promoter,promoters)



promfreq<-data.frame(table(proms.c$Promoter))
colnames(promfreq)<-c('Promoter','Freq')
#View(promfreq)


contrasts<-read.contrasts('Summary/!Contrasts.xlsx','Contrasts_values',c('PC','PM','PGC','PGM'))
contrasts.table<-contrasts$Contrasts.table
contr.matrix<-contrasts$Contrasts.matrix
contrasts.table
contr.matrix


contrasts.table[grepl('_adj',contrasts.table$Contrast),]

adjustments.list<-contrasts.table[grepl('_adj',contrasts.table$Contrast),'Contrast']


contrasts.adj<-adjust.contrast(contr.matrix,contrasts.table,adjustments.list,subset(stats.cast,!Promoter %in% dupl.gn.promoters ),suffix='_Mean')



for (tp in c('Matrix','Table')){
  contrasts.adj[[tp]]['PM-PC-(PGM-PGC)_adj_m','PM']<- contrasts.adj[[tp]]['PM-PC','PM']
  contrasts.adj[[tp]]['PM-PC-(PGM-PGC)_adj_m','PC']<- contrasts.adj[[tp]]['PM-PC','PC']
  contrasts.adj[[tp]]['PM-PC-(PGM-PGC)_adj_m','PGM']<- -contrasts.adj[[tp]]['PGM-PGC_adj','PGM']
  contrasts.adj[[tp]]['PM-PC-(PGM-PGC)_adj_m','PGC']<- -contrasts.adj[[tp]]['PGM-PGC_adj','PGC']
  contrasts.adj[[tp]]['PM-PC-(PGM-PGC)_adj_m','m']<- -contrasts.adj[[tp]]['PGM-PGC_adj','m']
}






contrasts.adj$Matrix
contrasts.adj$Table





allresults<-hypothesise(lmshape,promoters,contrasts.adj$Matrix,"0+Type")



results.t<-merge(allresults$All,contrasts.table[,c('Contrast','Description','Contrast_type')],by='Contrast',all.x=TRUE)
results.t<-rename(results.t,c("Variable"="Promoter"))



results.m<-melt(results.t,measure.vars=c('logFC','p.value','SE','t.value','PE','NE','FDR','logFDR'),variable.name = 'Stat',value.name = 'Value')

results.a<-merge(annotation.u,results.m,by='Promoter',all.y=TRUE)


head(results.a)

results.castfull<-dcast(results.a,Promoter+Strand+Reports_genes+Gene+Operon+Operon_name+Duplicated_reporting~Contrast+Stat,value.var = 'Value',drop=TRUE)
results.cast<-dcast(subset(results.a,Stat %in% c('logFC','FDR')),Promoter+Strand+Reports_genes+Gene+Operon+Operon_name+Duplicated_reporting~Contrast+Stat,value.var = 'Value',drop=TRUE)

head(results.cast)

results.cast.e<-as.data.frame(splitstackshape::cSplit(results.cast,splitCols = 'Reports_genes',sep=',',direction='long',drop = FALSE))

head(results.cast.e)


results.cast.ec<-results.cast.e[,c('Reports_genes','PM-PC_logFC','PGM-PGC_adj_logFC',"PM-PC-(PGM-PGC)_adj_m_logFC")]
head(results.cast.ec)

write.table(results.cast.ec,paste(odir,'/All_results_sidebyside.tsv',sep=''),
            row.names = FALSE,col.names=FALSE,sep='\t',quote = FALSE)



head(results.t)
# results.exp<-results.t[,c('Contrast','Description','Contrast_type','Plate','Well','Group','Metabolite','logFC','SE','PE','NE','t.value','p.value','FDR')]
# head(results.exp)


write.csv(results.t,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)


write.csv(results.cast.e,paste(odir,'/All_results_sidebyside_reporting.csv',sep=''),row.names = FALSE)

#TF enrichment start

TFvenn<-list('TF-targets'=unique(TF$Gene),'Reported_genes'=unique(as.character(results.cast.e$Reports_genes)) )

setdiff(promoters,unique(TF$Gene))


library('Vennerable')
plot(Venn(TFvenn),show = list(Faces = FALSE),
     doWeights = FALSE)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_diagram_RegulonDB_Alon_overlap_after_reannotation.pdf",sep=''),width=6,height=6)



#TF filter
TFc<-subset(TF,Effect %in% c('+','-'))

proms.e<-as.data.frame(splitstackshape::cSplit(proms.lm,splitCols = 'Reports_genes',sep=',',direction='long',drop = FALSE))
proms.e$Gene<-NULL
proms.e$Operon<-NULL

dim(proms.e)
tfproms<-merge(TFc,proms.e,by.x='Gene',by.y='Reports_genes')
dim(tfproms)

head(tfproms)





#TF effect
tfproms$TFReplicate<-makereplicates(tfproms[,c('TF','Gene','Promoter','Effect','Replicate','Group')])
tfproms$logGFPOD_effect<-ifelse(tfproms$Effect=='+',tfproms$logGFPOD_norm,-tfproms$logGFPOD_norm)


#Try using comparison results instead of raw data

tfshape<-dcast(subset(tfproms,!Promoter %in% refproms),Type+Group+Replicate+Drug+Sugar+Promoter+Gene+Effect+TFReplicate~TF,value.var = 'logGFPOD_effect',drop=TRUE)

View(tfshape)


tfactors<-unique(as.character(tfproms$TF))


tfresults<-hypothesise(tfshape,tfactors,contrasts.adj$Matrix,"0+Type")

tfresults.t<-merge(tfresults$All,contrasts.table[,c('Contrast','Description','Contrast_type')],by='Contrast',all.x=TRUE)
tfresults.t<-rename(tfresults.t,c("Variable"="TF"))
tfresults.t$FDRStars<-stars.pval(tfresults.t$FDR)

head(tfresults.t)

subset(tfresults.t,TF=='CRP' & Contrast_type!='Absolute')

subset(tfresults.t,TF=='CRP' & Contrast %in% selcont)

head(tfresults.t)




tfresults.m<-melt(tfresults.t,id.vars=c('TF','Contrast'),measure.vars=c('logFC','p.value','SE','t.value','PE','NE','FDR','logFDR'),variable.name = 'Stat',value.name = 'Value')
tfresults.m

tfresults.castfull<-dcast(tfresults.m,TF~Contrast+Stat,value.var = 'Value',drop=TRUE)
tfresults.cast<-dcast(subset(tfresults.m,Stat %in% c('logFC','FDR')),TF~Contrast+Stat,value.var = 'Value',drop=TRUE)


max(tfresults.t$logFC)

min(tfresults.t$logFC)

head(tfresults.cast)

#results.cast.e<-as.data.frame(splitstackshape::cSplit(results.cast,splitCols = 'Reports_genes',sep=',',direction='long',drop = FALSE))



selconts<-c('PM-PC','PGC-PC','PGM-PGC_adj',"PM-PC-(PGM-PGC)_adj_m")



tfhdata<-dcast(subset(tfresults.m,Stat %in% c('logFC','FDR') & Contrast %in% selconts),TF~Contrast+Stat,value.var = 'Value',drop=TRUE)
rownames(tfhdata)<-tfhdata$TF


#TF effect heatmap

nstep<-10
bgg <- colorRampPalette(c("blue4", "gray90", "red4"))(n = nstep)




dim(tfhdata)
tfhdata.sel<-tfhdata[apply(tfhdata[,grep('_FDR', colnames(tfhdata))], 1, function(x) {any(as.numeric(x) < 0.05,na.rm=TRUE)}),]
dim(tfhdata.sel)


hmap<-heatmap.2(data.matrix(tfhdata.sel[,grep('_logFC', colnames(tfhdata.sel))]),
          key=FALSE, #Colour key needs to be provided separately
          Colv=FALSE, #Columns should not be reordered
          Rowv=TRUE,
          trace='none',
          col=bgg, # Colour scale
          #breaks=brks, # Colour breaks
          xlab='Comparison',
          dendrogram='row', #Row dendogram, but should be changed to none, as dendrogram represents our data with filled-in values
          scale="none", #Should values be normalised in rows or columns - No
          na.color="white", # What colour to use with not missing values
          symkey=FALSE, #Provided colour scale is not symetrical
          reorderfun=reorderfun_mean,
          cexRow=0.7, #Some figure scaling parameters. Works only after a lot of experimentation
          cexCol=0.7,
          margin=c(10,20),
          lwid=c(0.2,0.8),
          lhei=c(0.05,0.95))

TFsel<-tfhdata.sel[hmap$rowInd,'TF']

tfresults.sel<-subset(tfresults.t,TF %in% TFsel & Contrast %in% selconts)

tfresults.sel$TF<-factor(tfresults.sel$TF, levels=TFsel, labels=TFsel)

max(tfresults.sel$logFC)
min(tfresults.sel$logFC)



amp<-3
brks<-seq(-amp,amp,by=1)
ggplot(tfresults.sel,aes(x=Contrast,y=TF))+
  geom_tile(aes(fill=logFC))+
  #geom_point(aes(size=FDRc,colour=logFC),alpha=0.9)+
  theme_minimal()+
  geom_text(aes(label=as.character(FDRStars)),nudge_y=-0.1)+
  ylab('')+
  #scale_size_discrete(range = c(2,4))+#,breaks=brks
  scale_fill_gradientn(colours = bgg,
                       breaks=brks,limits=c(-amp,amp))+
  theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Heatmap_TF_effect.pdf',sep = ''),
             width=6,height=15)






#TF Enrichment

resocls<-setdiff(colnames(results.cast.e),c('Gene','Operon'))

TFa<-merge(TFu,results.cast.e[,resocls],by.x='Gene',by.y='Reports_genes',all.y=TRUE)

#Only identified and measured
TFa<-subset(TFa,! (is.na(TF) | is.na(`PM-PC_logFC`))  )

head(TFa)

idvariables<-c('TF','Gene')
selstats<-c('logFC','FDR')

TFa.m<-enrichment.melt(TFa,idvariables,selstats)

TFa.en<-enrichment(TFa.m,terms = c('TF'),IDs = 'Gene',comparisons = 'Contrast',change = 'logFC',sign = 'FDR')

write.csv(TFa.en,paste0(odir,'/TF_enrichment.csv'))


selcont<-c('PM-PC','PGC-PC','PGM-PGC_adj',"PM-PC-(PGM-PGC)_adj_m")


TFa.en.c<-dcast(subset(TFa.en,Contrast %in% selcont),TF+Test~Contrast,value.var = 'FDR')


write.csv(TFa.en.c,paste0(odir,'/TF_enrichment_side-by-side.csv'))




#Enrichment heatmap
gyrs<-colorRampPalette(c("gray90","steelblue1","blue4"))(n = 5)
enbrks<-c(0,-log(0.05,10),2,3,4,5)

seldata<-subset(TFa.en.c,Test=='All')
rownames(seldata)<-seldata$TF


selcont<-c('PM-PC','PGC-PC','PGM-PGC_adj') #,"PM-PC-(PGM-PGC)_adj_m"


hdata<- -log10(seldata[,selcont])
hdata<-hdata[!apply(hdata, 1, function(x) {all(as.numeric(x) < -log(0.05,10),na.rm=TRUE)}),]


dim(hdata)

heatmap.2(data.matrix(hdata),
          key=FALSE, #Colour key needs to be provided separately
          Colv=FALSE, #Columns should not be reordered
          Rowv=TRUE,
          trace='none',
          col=gyrs, # Colour scale
          breaks=enbrks, # Colour breaks
          xlab='Comparison',
          dendrogram='row', #Row dendogram, but should be changed to none, as dendrogram represents our data with filled-in values
          scale="none", #Should values be normalised in rows or columns - No
          na.color="white", # What colour to use with not missing values
          symkey=FALSE, #Provided colour scale is not symetrical
          reorderfun=reorderfun_mean,
          cexRow=0.7, #Some figure scaling parameters. Works only after a lot of experimentation
          cexCol=0.7,
          margin=c(10,20),
          lwid=c(0.2,0.8),
          lhei=c(0.05,0.95))

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/TF_enrichment_All.pdf',sep = ''),
             width=6,height=12)







#Ecocyc enrichment
ecyc<-read.table('../Annotations/Ecoli/EcoCyc_Patwhays.tsv',header=TRUE,sep='\t')
ecyc<-data.frame(splitstackshape::cSplit(ecyc,'Link',sep='='))
ecyc$ID<-ecyc$Link_3
ecyc<-ecyc[,-grep('Link',colnames(ecyc))]

ecyc.m<-data.frame(splitstackshape::cSplit(ecyc,'Genes',sep=';',direction = 'long'))
ecyc.m<-rename(ecyc.m,c('Genes'='Gene'))

head(ecyc.m)


selrescol<-setdiff(colnames(results.cast.e),c('Gene','Operon'))

ece<-merge(ecyc.m,results.cast.e[,selrescol],by.x='Gene',by.y='Reports_genes',all.x=TRUE)


ece.c<-subset(ece,!is.na(Promoter))

head(ece)

idvariables<-c('Pathway','ID','Gene')
selstats<-c('logFC','FDR')

ece.m<-enrichment.melt(ece.c,idvariables,selstats)


View(ece)

subset(ece,ID=='TCA')


head(ece.m)


ece.en<-enrichment(ece.m,terms = c('ID','Pathway'),IDs = 'Gene',comparisons = 'Contrast',change = 'logFC',sign = 'FDR')


head(ece.en)
subset(ece.en,FDR<0.05 & Contrast=='PM-PC')

subset(ece.en,FDR<0.05 & Contrast=='PGC-PC')



write.csv(ece.en,paste0(odir,'/EcoCyc_enrichment.csv'))










#Figures
head(results.cast)


comparisons<-c('PM-PC','PGC-PC','PGM-PGC_adj',"PM-PC-(PGM-PGC)_adj_m")


head(results)
ggplot(results.t,aes(x=-log10(p.value),y=logFDR))+
  geom_vline(xintercept=-log10(0.05),color='red')+
  geom_hline(yintercept=-log10(0.05),color='red')+
  geom_point()+
  facet_wrap(~Contrast,ncol=4)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/FDR_adjustment_All.pdf',sep = ''),
             width=9,height=9,useDingbats=FALSE)




ggplot(subset(results.t,Contrast %in% comparisons),aes(x=logFC))+
  geom_histogram(position='identity',alpha=0.5)+
  facet_wrap(~Contrast)




ggplot(results.cast,aes(x=PGC_logFC,y=PGM_logFC))+
  geom_abline(aes(intercept=0,slope=1),color='red')+
  geom_smooth(method='lm')+
  geom_point(aes(color=abs(`PGM-PGC_adj_logFC`),
                 size=abs(`PGM-PGC_adj_logFC`)))+
  geom_smooth(method='lm')+
  geom_text(aes(label=ifelse(`PGM-PGC_adj_FDR`<0.05,as.character(Reports_genes),'')))


ggplot(results.cast,aes(y=PGC_logFC,x=PC_logFC))+
  geom_abline(aes(intercept=0,slope=1),color='red')+
  geom_smooth(method='lm')+
  geom_point(aes(color=abs(`PGC-PC_adj_logFC`),
                 size=abs(`PGC-PC_adj_logFC`)))+
  geom_smooth(method='lm')+
  geom_text(aes(label=ifelse(`PGC-PC_adj_FDR`<0.05,as.character(Reports_genes),'')),size=4)




erralpha<-1
errcolor<-'grey80'


brks<-c(0,1,2,3,4)

gradcols<-c('black','purple','purple')

xvars<-c('PC','PC','PGC')
yvars<-c('PM','PGC','PGM')
intvars<-c('PM-PC','PGC-PC','PGM-PGC_adj')

comps<-data.frame('X'=xvars,'Y'=yvars,'Interaction'=intvars)


for (com.id in 1:nrow(comps)){
  x<-as.character(comps[com.id,'X'])
  y<-as.character(comps[com.id,'Y'])
  intr<-as.character(comps[com.id,'Interaction'])
  xvar<-paste('`',x,'_logFC`',sep='')
  yvar<-paste('`',y,'_logFC`',sep='')
  intvar<-paste('`',intr,'_logFC`',sep='')
  intfdr<-paste('`',intr,'_FDR`',sep='')
  
  xPE<-paste('`',x,'_PE`',sep='')
  yPE<-paste('`',y,'_PE`',sep='')
  xNE<-paste('`',x,'_NE`',sep='')
  yNE<-paste('`',y,'_NE`',sep='')
  
  fit<-lm(as.formula(paste(yvar,'~',xvar,sep='')),results.cast)
  lmeq<-lm_eqn(fit)
  a<-fit$coefficients[[2]]
  b<-fit$coefficients[[1]]
  maincomp<-intvar
  
  scat<-ggplot(results.castfull,aes_string(y=yvar,x=xvar))+
    geom_vline(xintercept = 0,color='red',alpha=0.5,linetype='longdash')+
    geom_hline(yintercept = 0,color='red',alpha=0.5,linetype='longdash')+
    geom_abline(aes(intercept=b,slope=a),alpha=1,color='red')+
    geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
    geom_errorbar(aes_string(ymin=yNE,ymax=yPE),alpha=erralpha,color=errcolor,width=0)+
    geom_errorbarh(aes_string(xmin=xNE,xmax=xPE),alpha=erralpha,color=errcolor,height=0)+
    ggtitle('Scatterplot of contrasts and trend adjusted interaction',
            subtitle = 'Promoters with FDR<0.05 and |logFC|>1 are marked')+
    geom_point(aes(size=abs(eval(parse(text = intvar))),
                   color=abs(eval(parse(text = intvar)))))+
    scale_colour_gradientn(colours = gradcols,
                           breaks=brks,limits=c(0,3),name=maincomp)+
    scale_size(range = c(0.25, 4),name=maincomp)+
    annotate('text',x = 2, y =-0.5, label = lmeq$Full, parse = TRUE,color ='red',size=3)+
    scale_x_continuous(breaks=seq(-10,10,by=1))+
    scale_y_continuous(breaks=seq(-10,10,by=1))+
    geom_text_repel(aes(label=ifelse(eval(parse(text = intfdr)) <0.05 & abs(eval(parse(text = intvar)))>1,
                                     as.character(Reports_genes),'')),
                    size=4,nudge_y = 0.3,
                    force=1,
                    segment.colour=errcolor)
  
  
  fname<-paste(odir,"/Scatter_",x,"_",y,".pdf",sep = '')
  print(fname)
  cairo_pdf(fname,width=9,height=6)
  print(scat)
  dev.off()
}



#Does not work
gradcols2<-c('black','purple','purple')

xvars<-c('PM-PC','PM-PC')
yvars<-c('PGM-PGC','PGM-PGC')
intvars<-c('PM-PC-(PGM-PGC)','PM-PC-(PGM-PGC)_adj_m')

comps<-data.frame('X'=xvars,'Y'=yvars,'Interaction'=intvars)

for (com.id in 1:nrow(comps)){
  x<-as.character(comps[com.id,'X'])
  y<-as.character(comps[com.id,'Y'])
  intr<-as.character(comps[com.id,'Interaction'])
  xvar<-paste(x,'_logFC',sep='')
  yvar<-paste(y,'_logFC',sep='')
  intvar<-paste(intr,'_logFC',sep='')
  intfdr<-paste(intr,'_FDR',sep='')
  
  xPE<-paste(x,'_PE',sep='')
  yPE<-paste(y,'_PE',sep='')
  xNE<-paste(x,'_NE',sep='')
  yNE<-paste(y,'_NE',sep='')
  
  sigproms<-results.castfull[(results.castfull[,xvar]<0.05 | results.castfull[,yvar]<0.05 ) & sqrt(results.castfull[,xvar]^2+results.castfull[,yvar]^2)>0.75,'Promoter']
  
  
  scat<-ggplot(results.castfull,aes(y=results.cast[,yvar],x=results.cast[,xvar]))+
    geom_vline(xintercept = 0,color='red',alpha=0.5,linetype='longdash')+
    geom_hline(yintercept = 0,color='red',alpha=0.5,linetype='longdash')+
    geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
    geom_errorbar(aes(ymin=results.cast[,yNE],ymax=results.cast[,yPE]),alpha=erralpha,color=errcolor,width=0)+
    geom_errorbarh(aes(xmin=results.cast[,xNE],xmax=results.cast[,xPE]),alpha=erralpha,color=errcolor,height=0)+
    ggtitle('Scatterplot illustrating difference in response to\n100mM Metformin between bacteria grown on NGM (PM, PC) and NGM + 0.2% Glucose (PGM, PGC)',
            subtitle = 'Promoters with FDR<0.05 in at least one condition and logFC>0.5 distance from center are marked')+
    geom_point(aes(size=sqrt(results.cast[,xvar]^2+results.cast[,yvar]^2),
                   color=sqrt(results.cast[,xvar]^2+results.cast[,yvar]^2)))+
    xlab(xvar)+
    ylab(yvar)+
    scale_size(range=c(0.1,6),
               breaks=c(0,1,2,3,4,5,6),
               guide="legend",
               name=maincomp)+
    scale_colour_gradientn(colours = gradcols2,
                           breaks=brks,limits=c(0,4),guide="legend",name=maincomp)+
    annotate('text',x = 2, y =-0.5, label = lmeq$Full, parse = TRUE,color ='red',size=3)+
    scale_x_continuous(breaks=seq(-14,14,by=1),limits=c(-3,3))+
    scale_y_continuous(breaks=seq(-14,14,by=1),limits=c(-3,3))+
    guides(color=guide_legend(), size = guide_legend())+
    geom_text_repel(aes(label=ifelse(Promoter %in% sigproms,
                                     as.character(Reported_genes),'')),
                    force=1,
                    segment.colour=errcolor)
  
  
  fname<-paste(odir,"/Scatter_",intvar,".pdf",sep = '')
  print(fname)
  cairo_pdf(fname,width=15,height=12)
  print(scat)
  dev.off()
}




components<-c('PC','PM','PGC','PGM')
simple<-subset(results.t,Contrast %in% components)
allcomps<-subset(results.t,Contrast %in% comparisons)


#Maybe it makes sense to take global estiamtes of variability
ggplot(simple,aes(x=logFC,y=SE))+
  geom_point()+
  facet_wrap(~Contrast)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/MVA_separate.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)


ggplot(simple,aes(x=logFC,y=SE,color=Contrast))+geom_point()
dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/MVA.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)




#Volcano plots
# xvars<-c('PC','PC','PGC','PM-PC')
# yvars<-c('PM','PGC','PGM','PGM-PGC')


maincomps<-c('PM-PC','PGC-PC','PGM-PGC','PM-PC-(PGM-PGC)','PGM-PGC_adj','PM-PC-(PGM-PGC)_adj_m')

maincomp<-'logFC'

for (x in maincomps){
  xvar<-paste('`',x,'_logFC`',sep='')
  yvar<-paste('`',x,'_logFDR`',sep='')
  xPE<-paste('`',x,'_PE`',sep='')
  xNE<-paste('`',x,'_NE`',sep='')
  
  volcano<-ggplot(results.castfull,aes_string(x=xvar,y=yvar))+
    geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
    # geom_vline(xintercept = 1,color='red',alpha=0.5,linetype='longdash')+
    # geom_vline(xintercept = -1,color='red',alpha=0.5,linetype='longdash')+
    geom_errorbarh(aes_string(xmin=xNE,xmax=xPE),
                   alpha=erralpha,color=errcolor,height=0)+
    geom_point(aes(size=abs(eval(parse(text = xvar))),
                   color=abs(eval(parse(text = xvar)))))+
    scale_x_continuous(breaks=seq(-14,14,by=1),limits=c(-3,3))+
    scale_y_continuous(breaks=seq(0,14,by=1),limits=c(0,14))+
    scale_colour_gradientn(colours = gradcols,
                           breaks=brks,limits=c(0,3),name=maincomp)+
    scale_size(range = c(0.25, 4),name=maincomp)+
    geom_text_repel(aes(label=ifelse(eval(parse(text = yvar)) > -log10(0.05)+1,
                                     as.character(Reports_genes),'')),
                    size=4,
                    nudge_y = 0.3,
                    segment.colour = errcolor)+
    coord_cartesian(ylim=c(0, 14),xlim=c(-3,3))
  
  fname<-paste(odir,"/Volcano_",x,".pdf",sep = '')
  print(fname)
  cairo_pdf(fname,width=9,height=6)
  print(volcano)
  dev.off()
}




#Heatmap

#Filtering in two steps
#Differences heatmap

gyrs<-colorRampPalette(c("gray90","steelblue1","blue4"))(n = 6)

nstep<-20

bgg <- colorRampPalette(c("blue4", "gray90", "red4"))(n = nstep)
brks<-seq(-5,5,by=10/nstep)

brks2<-seq(-4,4,by=8/nstep)



head(results.m)

results.msel<-subset(results.m,Contrast %in% comparisons & Stat %in% c('logFC','FDR'))
results.scast<-dcast(results.msel,Promoter~Contrast+Stat,value.var = 'Value')

rownames(results.scast)<-results.scast$Promoter
results.scast$Promoter<-NULL
head(results.scast)


h.fdr<-colnames(results.scast)[grep('_FDR',colnames(results.scast))]
h.logFC<-colnames(results.scast)[grep('_logFC',colnames(results.scast))]


h.afdr<-setdiff(h.fdr,h.fdr[! grep('_adj',h.fdr) ])
h.alogFC<-setdiff(h.logFC,h.logFC[grep('_adj',h.logFC)])

h.alogFC

hrows<-results.scast[,'PM-PC_FDR']<0.05 & results.scast[,'PGM-PGC_FDR']>0.05# ,]


hrows<-apply(results.scast[,h.afdr], 1, function(x)(any(x< 0.05,na.rm=TRUE)) ) &
  apply(results.scast[,h.alogFC], 1, function(x)(any(abs(x)>1,na.rm=TRUE)))

dim(results.scast)
hdata<-results.scast[hrows, h.alogFC]

dim(hdata)
min(hdata)
max(hdata)

#Leave only metabolite with at least one significant change that's: |logFC|>1

hmap<-heatmap.2(as.matrix(hdata),key=TRUE,Colv=FALSE,trace='none',col=bgg,
                xlab='Comparison',Rowv=TRUE,
                breaks = brks,
                dendrogram="row",scale="none",na.color="white",
                cexRow=0.8,cexCol=0.5,margin=c(8,16),
                lwid=c(0.2,0.8),symkey=FALSE)

# heatmap.2(as.matrix(hdata),key=TRUE,Colv=FALSE,trace='none',col=bgg,
#           xlab='Comparison',Rowv=hmap$rowDendrogram,
#           breaks = brks,
#           dendrogram="row",scale="none",na.color="oldlace",
#           cexRow=0.5,cexCol=0.5,margin=c(8,16),
#           lwid=c(0.2,0.8),symkey=FALSE)

# dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Heatmap_Comparisons_logFCover1.pdf',sep = ''),
#              width=6,height=length(rownames(hdata))/12)


# dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Heatmap_Comparisons_logFCover1_Main.pdf',sep = ''),
#              width=4,height=length(rownames(hdata))/12)
# 



amp<-4
nstep<-amp*2



bgg <- colorRampPalette(c("blue", "gray90", "red"))(n = nstep)
brks<-seq(-amp,amp,by=amp*2/nstep)


results.sel<-subset(results.t,Contrast %in% gsub('_logFC','',colnames(hdata)) & Promoter %in% rownames(hdata))

results.sel$FDRstars<-stars.pval(results.sel$FDR)
results.sel$Promoter<-factor(results.sel$Promoter,levels=rownames(hdata),labels=rownames(hdata))




length(rownames(hdata))/8


ggplot(results.sel,aes(x=Contrast,y=Promoter))+
  geom_tile(aes(fill=logFC))+
  #geom_point(aes(size=FDRc,colour=logFC),alpha=0.9)+
  theme_minimal()+
  geom_text(aes(label=as.character(FDRstars)),nudge_y=-0.1)+
  ylab('')+
  #scale_size_discrete(range = c(2,4))+#,breaks=brks
  scale_fill_gradientn(colours = bgg,
                       breaks=brks,limits=c(-amp,amp))+
  theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

# dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Heatmap_Comparisons_Significance_logFCover1.pdf',sep = ''),
#              width=4,height=length(rownames(hdata))/8)


# dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Heatmap_Comparisons_Significance_UniquetoNGM.pdf',sep = ''),
#              width=4,height=length(rownames(hdata))/8)


#Venn

M.PMPC.up<-subset(results.cast,`PM-PC_FDR`<0.05 & `PM-PC_logFC`>0)[,'Promoter']
M.PMPC.down<-subset(results.cast,`PM-PC_FDR`<0.05 & `PM-PC_logFC`<0)[,'Promoter']

M.PGMPGC.up<-subset(results.cast,`PGM-PGC_adj_FDR`<0.05 & `PGM-PGC_adj_logFC`>0)[,'Promoter']
M.PGMPGC.down<-subset(results.cast,`PGM-PGC_adj_FDR`<0.05 & `PGM-PGC_adj_logFC`<0)[,'Promoter']


PGCPC.up<-subset(results.cast,`PGC-PC_FDR`<0.05 & `PGC-PC_logFC`>0)[,'Promoter']
PGCPC.down<-subset(results.cast,`PGC-PC_FDR`<0.05 & `PGC-PC_logFC`<0)[,'Promoter']

# M.PGMPGC.up<-subset(results.cast,`PGM-PGC_FDR`<0.05 & `PGM-PGC_logFC`>0)[,'Promoter']
# M.PGMPGC.down<-subset(results.cast,`PGM-PGC_FDR`<0.05 & `PGM-PGC_logFC`<0)[,'Promoter']






#Venn diagrams
library('Vennerable')

Met_NGMvsGlu<-list('NGM Down'=M.PMPC.down,
                   'NGM + Glu down'=M.PGMPGC.down,
                   'NGM Up'=M.PMPC.up,
                   'NGM + Glu Up'=M.PGMPGC.up)

NGMGluvsNGM<-list('NGM + Glu down'=PGCPC.down,
                  'NGM + Glu up'=PGCPC.up)

plot(Venn(Met_NGMvsGlu),show = list(Faces = FALSE),
     doWeights = FALSE,type='ellipses')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_diagram_Metformin_in_medias.pdf",sep=''),
             width=6,height=6)




subset(results.cast,Promoter %in% intersect(M.PMPC.up,M.PGMPGC.up))[,c('Promoter','Reports_genes','PM-PC_logFC','PM-PC_FDR','PGM-PGC_adj_logFC','PGM-PGC_adj_FDR')]


subset(results.cast,Promoter %in% intersect(M.PMPC.down,M.PGMPGC.down))[,c('Promoter','Reports_genes','PM-PC_logFC','PM-PC_FDR','PGM-PGC_adj_logFC','PGM-PGC_adj_FDR')]




subset(results.cast,Promoter %in% setdiff(M.PMPC.up,M.PGMPGC.up))[,c('Promoter','Reports_genes','PM-PC_logFC','PM-PC_FDR','PGM-PGC_adj_logFC','PGM-PGC_adj_FDR')]




subset(results.cast,`PM-PC_FDR`<0.05 & `PM-PC_logFC` >0 )[,c('Promoter','Reports_genes','PM-PC_logFC','PM-PC_FDR','PGM-PGC_adj_logFC','PGM-PGC_adj_FDR')]




# plot(Venn(NGMGluvsNGM),show = list(Faces = FALSE),
#      doWeights = FALSE)
# 
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram_Metformin_in_medias.pdf",sep=''),width=6,height=6)
# 







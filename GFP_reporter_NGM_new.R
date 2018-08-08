library(tidyverse)

library(ComplexHeatmap)
library(circlize)


library(ggrepel)

library(plot3D)
library(rgl)
library(pca3d)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)

theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau


setwd("~/Dropbox/Projects/2015-Metformin/GFP_reporters")




#load('GFP_reporters.RData')
#save.image('GFP_reporters.RData')

odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


rawdata<- read_csv("UAL_screen/UAL_all_data.csv") %>%
  rename(Group=Type) %>%
  mutate(Group=recode(Group,'NC'='PC','NM'='PM','NGC'='PGC','NGM'='PGM'))


#Change names to avoid confusion

head(rawdata)

#Explanation of sample names
group<-c('PC','PM','PGC','PGM')
sampledescriptions<-c('NGM',
                      'NGM + 100mM Metformin',
                      'NGM + 0.2% Glucose',
                      'NGM + 0.2% Glucose + 100mM Metformin')
metf<-c("0","100","0","100")
glu<-c("0","0","0.2%","0.2%")

explanation<-data.frame(Group=group,
                        Metformin_mM=metf,
                        Glucose=glu,
                        Sample_description=sampledescriptions)





#TFs

TF<-read.csv('../Annotations/Ecoli/Transcription_Factors/network_tf_gene.txt',comment.char = '#',sep = '\t',header=FALSE)

colnames(TF)<-c('TF','Gene','Effect','Evidence','Support','NAs')
TF$NAs<-NULL

#Remove duplicated effect descriptions
TFu<-TF[!duplicated(TF[,c('TF','Gene')]),]
dim(TFu)


#Annotation

annotation<-read_csv('../Annotations/Ecoli/UAL/UAL_reannotated.csv') %>%
  rename(Promoter=UAL_promoter) %>%
  mutate(Strain=ifelse(is.na(Strand),'1',Strand))


annotation.u<-annotation%>%
  filter(!duplicated(Promoter)) %>%
  select(-c(Plate,Well))
  


#Data
data<-rawdata %>%
  unite(ID,Group,Replicate,remove = FALSE) %>%
  unite(Index,Plate,Well,remove = FALSE) %>%
  left_join(explanation) %>%
  left_join(annotation) %>%
  spread(Reading,Value) %>%
  mutate_at(c("Replicate","Promoter"),as.factor ) %>%
  mutate(GFPOD=GFP/OD,
         logGFPOD=log2(GFPOD),
         logGFPOD=ifelse(logGFPOD<0,0,logGFPOD),
         Group=factor(Group,levels=group),
         Plate=factor(Plate,levels=as.character(1:21),labels=as.character(1:21) ),
         Promoter=relevel(Promoter,ref = "U139")) %>%
  filter(!Promoter %in% c('Empty',NA)) %>%
  #As plate 21 does not have promoterless plasmid strains, use Replicate, Group average for each plasmid
  left_join(x=., y= filter(.,Promoter %in% c('U139','U66')) %>% group_by(Group,Strand,Replicate) %>% summarise(logGFPOD_ref=mean(logGFPOD),Plate="21") %>% ungroup  ) %>%
  group_by(Group,Plate,Replicate,Strand) %>%
  mutate(logGFPOD_ref=ifelse(Plate=="21",logGFPOD_ref,logGFPOD[Promoter %in%  c('U139','U66') ]),
         logGFPOD_norm=logGFPOD-logGFPOD_ref) %>%
  group_by(Promoter,Group,Replicate) %>%
  mutate(PReplicate=row_number()) %>%
  ungroup %>%
  filter( OD>0.25 | Promoter %in% c('U139','U66'))


head(data)
data %>%
  filter(Promoter %in% c('U139','U66')) %>%
  ggplot(aes(x=Plate,y=logGFPOD,color=Group))+
  geom_boxplot()+
  ggtitle('Promoterless plasmid fluorescence by plate')+
  facet_wrap(~Promoter)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Boxplot_Promoterless_by_plate.pdf",sep=''),
             width=9,height=4, useDingbats=FALSE)


data %>%
  filter(Promoter %in% c('U139','U66')) %>%
  ggplot(aes(x=Replicate,y=logGFPOD,color=Group))+
  geom_boxplot()+
  ggtitle('Promoterless plasmid fluorescence by replicate')+
  facet_wrap(~Promoter)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Boxplot_Promoterless_by_replicate.pdf",sep=''),
             width=9,height=4, useDingbats=FALSE)


data %>%
  filter(Promoter %in% c('U139','U66')) %>%
  ggplot(aes(x=Group,y=logGFPOD,color=Group))+
  geom_boxplot()+
  ggtitle('Promoterless plasmid fluorescence by type')+
  facet_wrap(~Promoter)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Boxplot_Promoterless_by_type.pdf",sep=''),
             width=9,height=4, useDingbats=FALSE)


#Investigate growth deviations


data %>%
  ggplot(aes(x=OD,fill=Group,color=Group))+
  geom_histogram(position='identity',alpha=0.5)+
  facet_grid(Replicate~Group,labeller = label_both)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Distribution_OD.pdf",sep=''),
             width=9,height=6)



data %>%
  filter(OD<0.25) %>%
  group_by(Promoter) %>%
  summarise(Count=n()) 





stats<-data %>%
  gather(Measure,Value,GFP,OD,GFPOD,logGFPOD,logGFPOD_ref,logGFPOD_norm) %>%
  group_by(Promoter,Group,Measure,Reports_genes,Gene,Operon,Operon_name,Duplicated_reporting) %>%
  summarise(Mean=mean(Value),
            SD=sd(Value)) %>%
  ungroup


  
stats.m<-multiplex(stats,c("Promoter","Group","Reports_genes","Gene","Operon","Operon_name","Duplicated_reporting"),dims=2)


stats.m %>%
  filter(x_Measure=="OD" & y_Measure %in% c("logGFPOD","logGFPOD_norm") ) %>%
  ggplot(aes(x=x_Mean,y=y_Mean,color=Group))+
  #ylim(0,2)+
  xlab('OD')+
  ylab("Measure")+
  geom_point(size=1)+
  scale_y_continuous(breaks=seq(-4,20,by=1))+
  ggtitle('Raw data before and after normalisation')+
  geom_smooth(method='lm')+
  facet_grid(y_Measure~Group,scale="free_y",space="free_y")

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/GFP_OD_relationship.pdf",sep=''),
             width=9,height=4)


data %>%
  ggplot(aes(x=logGFPOD_norm,fill=Replicate,color=Replicate))+
  geom_vline(xintercept = 0,color='blue')+
  geom_histogram(position='identity',alpha=0.5)+
  facet_grid(Replicate~Group,labeller = label_both)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Distribution_in_replicates.pdf",sep=''),
             width=9,height=6)





refdat<-stats %>%
  filter(Promoter %in% c('U139','U66') & Measure %in% c('logGFPOD','logGFPOD_norm') )


stats %>%
  filter(Measure %in% c('logGFPOD','logGFPOD_norm') ) %>%
  ggplot(aes(x=Mean,fill=Group))+
  geom_histogram(position='identity',alpha=0.5)+
  xlab('logGFPOD')+
  geom_vline(data=refdat,aes(xintercept=Mean,color=Group))+
  #geom_text(data=refdat,aes(label=Group,color=Group,x=Mean+0.5,y=(as.numeric(Group)+1)*20+750),size=2)+
  ggtitle('Promoter activity distribution before and after promoterless plasmid normalisation')+
  facet_grid(.~Measure,scale="free_x")

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Distribution_raw_normalised.pdf',sep = ''),
             width=6,height=4,useDingbats=FALSE)




testproms<-c('U66','U139','lacZ','crp','serA','wrbA','gadW','gadA','gadX','gadB','ompA','rplN')

stats %>%
  filter(Promoter %in% testproms & Measure %in% c('logGFPOD','logGFPOD_norm') ) %>%
  mutate(Promoter=factor(Promoter,levels=testproms)) %>%
  ggplot(aes(y=Mean,x=Promoter,color=Group))+
  geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD))+
  geom_point()+
  scale_y_continuous(breaks=seq(0,14,by=1))+
  ylab('logGFPOD')+
  facet_grid(.~Measure,labeller=label_both,scale="free_y")


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Boxplot_Selected_comparison.pdf",sep=''),
             width=12,height=6, useDingbats=FALSE)




#Duplicated reporting comparison



dupl.data<-stats %>%
  filter(Measure=="logGFPOD_norm") %>%
  filter(!is.na(Duplicated_reporting) ) %>%
  mutate(IsOperon=!is.na(Operon_name)) %>%
  separate_rows(Duplicated_reporting)
           
  

gn<-subset(dupl.data,!IsOperon)[,c('Promoter','Group','Duplicated_reporting','Mean','SD')]
op<-subset(dupl.data,IsOperon)[,c('Group','Duplicated_reporting','Mean','SD')]


dupl.gn.promoters<-as.character(unique(gn$Promoter))

dr.sbs<-merge(gn,op,by=c('Group','Duplicated_reporting'),all=TRUE,suffixes=c("_gene","_operon"))


head(dupl.data)

dupl.data %>%
  ggplot(aes(y=Mean,x=Group,color=IsOperon))+
  geom_hline(yintercept=0,color='red',alpha=0.5)+
  geom_boxplot(aes(group=interaction(Type,Promoter,IsOperon)))+
  scale_y_continuous(breaks=seq(-5,14,by=1))+
  ylab('logGFPOD_norm')+
  facet_wrap(~Duplicated_reporting,ncol = 5)#,labeller=label_both


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Duplicated_reporting.pdf",sep=''),
             width=10,height=50, useDingbats=FALSE)

ggplot(dr.sbs,aes(x=Mean_gene,y=Mean_operon))+
  geom_errorbar(aes(ymin=Mean_operon-SD_operon,ymax=Mean_operon+SD_operon),color="gray80")+
  geom_errorbarh(aes(xmin=Mean_gene-SD_gene,xmax=Mean_gene+SD_gene),color="gray80")+
  scale_x_continuous(breaks=-10:10)+
  scale_y_continuous(breaks=-10:10)+
  coord_cartesian()+
  geom_hline(yintercept=0,color='red',alpha=0.5)+
  geom_vline(xintercept=0,color='red',alpha=0.5)+
  geom_point()+
  facet_grid(.~Group)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Operon_gene_activation_comparison.pdf",sep=''),
             width=10,height=10, useDingbats=FALSE)





bioinfo<-data %>%
  group_by(ID,Group,Glucose,Metformin_mM) %>%
  summarise %>%
  ungroup %>%
  data.frame
rownames(bioinfo)<-bioinfo$ID


PCAres<-data %>%
  filter(!Promoter %in% c('U66','U139') ) %>%
  mutate(Index=as.character(Index)) %>%
  group_by(Index) %>%
  filter(!n()<16 & !any(is.na(logGFPOD_norm))) %>%
  ungroup %>%
  PCAprep('ID','Index','logGFPOD_norm',bioinfo)


HC<-PCAres$HC
pca<-PCAres$pca
ellipses<-PCAres$Ellipses
pcadata<-PCAres$pcadata
pcaloadings<-PCAres$Loadings
pcashape=PCAres$pcashape





plot(HC, labels=samples,
     hang=-1, main=paste0("Cluster Dendrogram - ",type), xlab="", sub="", cex=1)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Hierarchical_Clustering_",type,".pdf",sep=''))

write.csv(pcaloadings,paste(odir,"/PCA_loadings_",type,".csv",sep=''))


plot(pca,type='l')

dev.copy2pdf(device=cairo_pdf,
             file=paste0(odir,"/PCA_components_",type,".pdf"),
             width=9,height=6)


ggplot(pcadata,aes(x=PC1,y=PC2,color=Glucose))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM=="0",Glucose, NA ) ) ),size=5,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"100"=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA.pdf",sep=''),
             width=12,height=9)

#Heatmap

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap.pdf",sep=''),
             width=6,height=9, useDingbats=FALSE)





#install.packages('kernlab')
library(kernlab)

clean_data.x<-as.matrix(clean_data)

PCdf<-data.frame(ir.pca$x)


groupCluster <- kmeans(clean_data, 4, nstart = 100)
clusters<-as.factor(groupCluster$cluster)


groupSpecGauss<-specc(clean_data.x,centers=4, kernel = "rbfdot",iterations=1000)#,iterations=100
clusters<-as.factor(groupSpecGauss)


PCdf$clusters<-as.factor(groupSpecGauss)
# 
# ggbiplot(ir.pca, obs.scale = 1,
#                        var.scale = 1,
#                        groups= clusters,
#                        #groups = sampleclasses,
#                        ellipse = TRUE,
#                        circle = TRUE,
#                        var.axes = 0)+
#   #   geom_vline(aes(xintercept=0),color='gray80')+
#   #   geom_hline(aes(yintercept=0),color='gray80')+
#   scale_color_discrete(name = '')+
#   ggtitle('PCA')+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())



#sampleclasses


#Needs to be done manually

# library(rgl)
# library(pca3d)

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




#Linear modelling starts here

contrasts<-read.contrasts('Summary/!Contrasts.xlsx')

contrasts.desc<-contrasts$Contrasts.table %>%
  select(-c(Target,Reference))

contr.matrix<-contrasts$Contrasts.matrix

contrasts.desc

contr.matrix.a<-contrast.adjust(contr.matrix,adjustments)

contr.matrix.a

#Use data with deviant growth removed
#Some promoters are present in multiple wells, nead unique ID
#Take gene values from operon based on duplicated reporting observations

data %>%
  group_by(Group) %>%
  summarise(Mean=mean(logGFPOD_norm,na.rm = TRUE),
            SD=sd(logGFPOD_norm,na.rm = TRUE))


results.all<-data %>%
  filter(!Promoter %in% c("U139","U66")) %>%
  filter(! ( !is.na(Duplicated_reporting) & is.na(Operon_name)) ) %>%
  group_by(Promoter,PReplicate,Reports_genes,Gene,Operon,Operon_name,Duplicated_reporting,Amplicons) %>%
  do(hypothesise(.,'logGFPOD_norm~0+Group',contr.matrix.a)) %>%
  getresults(contrasts.desc)


#Separate results
results<-results.all$results


results<-read_csv(paste0(odir,'/All_results.csv'))


results %>%
  #group_by(Promoter) %>%
  mutate(#FDR=p.adjust(p.value,method='fdr'),
         Up=logFC>0 & FDR<0.05,
         Down=logFC<0 & FDR<0.05,
         All=FDR<0.05) %>%
  gather(Type,Pass,Up,Down,All) %>%
  group_by(Contrast,Type) %>%
  summarise(Count=sum(Pass)) %>%
  spread(Type,Count)



results.cast<-results.all$cast 

results.cast.e<-results.cast %>%
  separate_rows(Reports_genes)

results.castfull<-results.all$castfull
results.multi<-results.all$multi


groupcols<-results.multi %>%
  select(Index:Amplicons) %>%
  colnames

results.multi2<-multiplex(results,groupcols,dim=2)



adj.cont<-c("PM-PC_adj","PGM-PGC_adj","PGC-PC_adj","PM-PC-(PGM-PGC)_adj","PM-PC-(PGM-PGC)_adj_m")
adj.x<-c("PC","PGC","PC","PGM-PGC_adj","PGM-PGC")
adj.y<-c("PM","PGM","PGC","PM-PC_adj","PM-PC")

adjust<-data.frame(Contrast=adj.cont,
                    x=adj.x,
                    y=adj.y)

adjustments<-adjust%>%
  group_by(Contrast,x,y) %>%
  do(covariation(adjustments=.,data=results.multi2)) %>%
  data.frame

adjustments


write.csv(results,paste0(odir,'/All_results.csv'),row.names = FALSE)
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


#Venn

results<-read_csv(paste0(odir,"/All_results.csv"))


venn.pass<-results %>%
  mutate(Up=logFC>0 & FDR<0.05,
         Down=logFC<0 & FDR<0.05,
         All=FDR<0.05) %>%
  #gather different thresholds
  gather(Type,Pass,Up,Down,All) %>%
  filter(Pass) %>%
  group_by(Contrast,Type) %>%
  do(List=c(as.character(.$Promoter)))


metf<-venn.pass %>%
  filter(Contrast %in% c("PM-PC","PGM-PGC_adj") & Type!="All") %>%
  unite(OT,Contrast,Type) %>%
  select(OT,List) %>%
  spread(OT,List) %>%
  as.list %>%
  unlist(recursive=FALSE)


vcols<-c("red","red4","blue","blue4")
grid::grid.draw(VennDiagram::venn.diagram(metf, col=vcols,cat.col=vcols ,NULL))

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Metformin_in_medias.pdf",sep=''),
             width=3,height=3, useDingbats=FALSE)
# dev.off()#


M.PMPC.up<-subset(results.cast,`PM-PC_FDR`<0.05 & `PM-PC_logFC`>0)[,'Promoter']
M.PMPC.down<-subset(results.cast,`PM-PC_FDR`<0.05 & `PM-PC_logFC`<0)[,'Promoter']

M.PGMPGC.up<-subset(results.cast,`PGM-PGC_adj_FDR`<0.05 & `PGM-PGC_adj_logFC`>0)[,'Promoter']
M.PGMPGC.down<-subset(results.cast,`PGM-PGC_adj_FDR`<0.05 & `PGM-PGC_adj_logFC`<0)[,'Promoter']

PGCPC.up<-subset(results.cast,`PGC-PC_FDR`<0.05 & `PGC-PC_logFC`>0)[,'Promoter']
PGCPC.down<-subset(results.cast,`PGC-PC_FDR`<0.05 & `PGC-PC_logFC`<0)[,'Promoter']

# M.PGMPGC.up<-subset(results.cast,`PGM-PGC_FDR`<0.05 & `PGM-PGC_logFC`>0)[,'Promoter']
# M.PGMPGC.down<-subset(results.cast,`PGM-PGC_FDR`<0.05 & `PGM-PGC_logFC`<0)[,'Promoter']

#Venn diagrams
# library('Vennerable')
# 
# Met_NGMvsGlu<-list('NGM Down'=M.PMPC.down,
#                    'NGM + Glu down'=M.PGMPGC.down,
#                    'NGM Up'=M.PMPC.up,
#                    'NGM + Glu Up'=M.PGMPGC.up)
# 
# NGMGluvsNGM<-list('NGM + Glu down'=PGCPC.down,
#                   'NGM + Glu up'=PGCPC.up)
# 
# plot(Venn(Met_NGMvsGlu),show = list(Faces = FALSE),
#      doWeights = FALSE,type='ellipses')

# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Venn_diagram_Metformin_in_medias.pdf",sep=''),
#              width=6,height=6)




subset(results.cast,Promoter %in% intersect(M.PMPC.up,M.PGMPGC.up))[,c('Promoter','Reports_genes','PM-PC_logFC','PM-PC_FDR','PGM-PGC_adj_logFC','PGM-PGC_adj_FDR')]


subset(results.cast,Promoter %in% intersect(M.PMPC.down,M.PGMPGC.down))[,c('Promoter','Reports_genes','PM-PC_logFC','PM-PC_FDR','PGM-PGC_adj_logFC','PGM-PGC_adj_FDR')]




subset(results.cast,Promoter %in% setdiff(M.PMPC.up,M.PGMPGC.up))[,c('Promoter','Reports_genes','PM-PC_logFC','PM-PC_FDR','PGM-PGC_adj_logFC','PGM-PGC_adj_FDR')]




subset(results.cast,`PM-PC_FDR`<0.05 & `PM-PC_logFC` >0 )[,c('Promoter','Reports_genes','PM-PC_logFC','PM-PC_FDR','PGM-PGC_adj_logFC','PGM-PGC_adj_FDR')]




# plot(Venn(NGMGluvsNGM),show = list(Faces = FALSE),
#      doWeights = FALSE)
# 
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram_Metformin_in_medias.pdf",sep=''),width=6,height=6)
# 






#TF effect - Up or down
TFc<-TF %>%
  filter(Effect %in% c('+','-') )

tfproms<-results %>%
  separate_rows(Reports_genes) %>%
  select(-c(Gene,Operon)) %>%
  left_join(TFc,by = c("Reports_genes"="Gene")) %>%
  mutate(logFC_effect=ifelse(Effect=='+',logFC,-logFC)) %>%
  group_by(Contrast,TF,Reports_genes,Promoter,Effect) %>%
  mutate(TFReplicate=row_number() ) %>%
  filter(!is.na(logFC))
  
# tfproms<-merge(TFc,proms.e,by.x='Gene',by.y='Reports_genes')


dim(tfproms)

head(tfproms)


library(broom)
tfresults<-tfproms %>%
  group_by(Contrast) %>%
  do(tidy(lm(logFC_effect ~ 0+TF, data=.)) ) %>%
  mutate(FDR=p.adjust(p.value,method='fdr'),
         TF=str_replace(term,'TF',''),
         term=NULL) %>%
  select(Contrast,TF,everything())



tfresults %>%
  filter(Contrast=="PM-PC_adj" & FDR<0.05) %>%
  View

#Try using comparison results instead of raw data

tfshape<-dcast(subset(tfproms,!Promoter %in% refproms),Type+Group+Replicate+Drug+Sugar+Promoter+Gene+Effect+TFReplicate~TF,value.var = 'logGFPOD_effect',drop=TRUE)

View(tfshape)

tfactors<-unique(as.character(tfproms$TF))

tfresults<-hypothesise(tfshape,tfactors,contrasts.adj$Matrix,"0+Group")

# tfresults.t<-merge(tfresults$All,contrasts.table[,c('Contrast','Description','Contrast_type')],by='Contrast',all.x=TRUE)
# tfresults.t<-rename(tfresults.t,c("Variable"="TF"))
# tfresults.t$FDRStars<-stars.pval(tfresults.t$FDR)

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

# write.csv(TFa.en,paste0(odir,'/TF_enrichment.csv'))


selcont<-c('PM-PC','PGC-PC','PGM-PGC_adj',"PM-PC-(PGM-PGC)_adj_m")


TFa.en.c<-dcast(subset(TFa.en,Contrast %in% selcont),TF+Test~Contrast,value.var = 'FDR')


# write.csv(TFa.en.c,paste0(odir,'/TF_enrichment_side-by-side.csv'))




TFen<-read_csv(paste0(odir,'/TF_enrichment.csv')) %>%
  select(-X1)



selcont<-c('PM-PC','PGC-PC','PGM-PGC_adj')

seldesc<-as.character(contrasts.desc[contrasts.desc$Contrast %in% selcont,'Description'])



#,"PM-PC-(PGM-PGC)_adj_m"



TForder<-c("EvgA","OmpR","PgrR","CytR","PhoP","PhoB","IclR",
           "LacI","NtrC","AgaR","Nac","MarA","Cra","Mlc","CueR",
           "ArgR","TorR","H-NS","SoxS","NagC","CecR","RbsR","IHF","GatR",
           "NanR","CusR","YedW","Ada","CysB","LexA","TreR","Zur","McbR",
           "GlpR","TrpR","IscR","NarL","FlhDC","NarP","ModE","ArcA","FNR","CRP")


enrbrks<-c(0,-log(0.05,10),2,3,4,200)
enrlbls<-c('N.S.','<0.05','<0.01','<0.001','<0.0001')

TFen %>%
  left_join(contrasts.desc) %>%
  filter(Test=="All" & Contrast %in% selcont ) %>% 
  mutate(logFDR=ifelse(FDR==1 ,0,-log10(FDR) ),
         logFDRbin=cut(logFDR, breaks=enrbrks,labels=enrlbls,right = FALSE),
         Description=factor(Description,levels=seldesc)) %>%
  group_by(TF) %>%
  filter( any(FDR<0.05)) %>%
  ungroup %>%
  #clustorder('TF',c("Contrast","Test"),'logFDR',descending=FALSE) %>% 
  mutate(TF=factor(TF,levels=TForder)) %>%
  group_by(TF) %>%
  filter( sum(FDR<0.05)>1 ) %>%
  ungroup %>%
  PlotEnrichment("Description","TF","logFDRbin",ncols = length(enrbrks)) +
  ylab("Transcription factor")


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/TF_enrichment_All_new.pdf',sep = ''),
             width=3,height=6)



length(unique(TFdat$TF))
length(TForder)

setdiff(unique(TFdat$TF),TForder)

setdiff(TForder,unique(TFdat$TF))



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




#Volcano plots
cnames<-c('SM-S'='OP50-C treatment','RM-R'='OP50-MR treatment','R-S'='Strain difference','SM-S-(RM-R)'='Longevity effect')


results <- read_csv('All_results.csv') %>%
  group_by(Contrast) %>%
  arrange(desc(logFDR)) %>%
  mutate(Name=ifelse(row_number()<=10 & FDR<0.05,Promoter,NA ))


unique(data$Comparison)

results %>%
  filter(Contrast %in% c('PM-PC','PGM-PGC_adj','PGC-PC')) %>%
  ggplot(aes(x=logFC,y=logFDR))+
  geom_hline(yintercept = -log10(0.05),color='red4',alpha=0.5)+
  geom_point(size=1,color='gray50' )+
  geom_text_repel(aes(label=Name),color='red2')+
  ylab('-log10(FDR)')+
  facet_wrap(~Description,ncol = 3)

ggsave(file=paste0(odir,'/UAL_Volcano_plot_clean.pdf'),
       width=120,height=60,units='mm',scale=2,device=cairo_pdf,family="Arial")






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







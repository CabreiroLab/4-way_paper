library(tidyverse)
library(ggrepel)

library(ComplexHeatmap)
library(circlize)

#library(RColorBrewer)


#library(limma)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)

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



#Get data

data<-read_csv('Data/Summary.csv') %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A')) %>%
  mutate(Sample=paste(Strain,ifelse(Type=='Control','C','T'),sep='_'),
         SampleID=paste(Sample,Replicate,sep='_'),
         Sample=factor(Sample,
                        levels=c("OP50Sens_C","OP50Sens_T"),
                        labels=c("OP50Sens_C","OP50Sens_T") ) ) %>%
  left_join(info[,c('Index','MetaboliteU')],by='Index') %>%
  select(File:Metformin_mM,Replicate,Well,Index,Sample:MetaboliteU,Metabolite=Name,EcoCycID:Group,G=Int_750nm_log,GR=a_log) %>%
  gather(Measure,Value,G,GR) %>%
  group_by(Plate,Type,Replicate,Group,Measure) %>%
  mutate(Value_ref=Value[Metabolite=='Negative Control'],
         Value_norm=ifelse(Metabolite=="Negative Control",Value,Value-Value_ref)) %>%
  ungroup %>%
  mutate_at(c('SampleID','Sample','Strain','Metformin_mM'),as.factor)


data %>%
  filter(Metabolite=="Negative Control")


View(data)

data.nc<-data %>%
  filter(Plate !='PM5')

PM1PM2<-subset(info,Plate %in% c('PM1','PM2A'))



#Unique sample descriptors
bioinfo<-data.nc %>%
  group_by(SampleID,Sample,Strain,Metformin_mM) %>%
  summarise %>%
  ungroup %>%
  data.frame

rownames(bioinfo)<-bioinfo$SampleID



PCAres<-PCAprep(data.nc %>% filter(Measure=="G"),"SampleID","Index","Value_norm", bioinfo )

HC<-PCAres$HC
pca<-PCAres$pca
ellipses<-PCAres$Ellipses
pcadata<-PCAres$pcadata
pcaloadings<-PCAres$Loadings
pcashape=PCAres$pcashape

plot(HC, labels=rownames(pcashape),
     hang=-1, main="Cluster Dendrogram", xlab="", sub="", cex=1)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Hierarchical_Clustering.pdf",sep=''),
             width=9,height=6)


write.csv(pcaloading,paste(odir,"/PCA_loadings.csv",sep=''))

plot(pca,type='l')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_components.pdf",sep=''),
             width=9,height=6)

ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',PCAres$PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PCAres$PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Sample),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=5,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2) )+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  #geom_text(aes(label=Sample))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA.pdf",sep=''),
             width=12,height=9)


data.nc %>%
  filter(Measure=="G") %>%
  HMap("SampleID","Index","Value_norm",
       bioinfo %>% select(Metformin_mM) %>% rename(`Metformin, mM`=Metformin_mM),
       cols=list("Metformin, mM"=c('50'='black','0'='white')) )

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap.pdf",sep=''),
             width=6,height=9, useDingbats=FALSE)





#Get LM ready
#Based on selected data
#Get negative control wells
head(data.nc)

contrasts<-read.contrasts('!Contrasts.xlsx')

contrasts.desc<-contrasts$Contrasts.table %>%
  select(Description:Strain)

contr.matrix<-contrasts$Contrasts.matrix
contr.matrix

contrasts.desc


#Carboxylic acids to remove
coxy<-c('Caproic Acid','Capric Acid','4-Hydroxy Benzoic Acid','2-Hydroxy Benzoic Acid')

selmets<-info %>%
  filter(!(Plate %in% c('PM3B','PM4A') &
             Metabolite %in% PM1PM2$Metabolite ) &
           !Metabolite %in% c('Negative Control',coxy))


results.all<-data.nc %>%
  filter(Metabolite !="Negative Control") %>%
  #filter(!MetaboliteU %in% selmets$MetaboliteU) %>%
  group_by(Measure,Group,Plate,Well,Index,Metabolite,MetaboliteU,EcoCycID,KEGG_ID) %>%
  do(hypothesise(.,'Value_norm~0+Sample',contr.matrix)) %>%
  getresults(contrasts.desc,c("Measure"))





#Unfiltered results
# results.all$results %>%
#   write.csv(paste0(odir,'/Ecoli_unfiltered_results.csv'),row.names = FALSE)


#Separate results
results<-results.all$results %>%
  filter(MetaboliteU %in% selmets$MetaboliteU)


results.cast<-results.all$cast %>% filter(Measure=='G') %>%
  filter(MetaboliteU %in% selmets$MetaboliteU)

results.castfull<-results.all$castfull %>%
  filter(MetaboliteU %in% selmets$MetaboliteU) %>%
  filter(Measure=='G') %>%
  mutate(Pole=ifelse(`C_logFC`>0,
                     pole2(`C_logFC`,`T_logFC`),
                     pole2(`C_logFC`,`T_logFC`)+0.5),
         Pole=ifelse(Pole>0.625,Pole-1,Pole),
         Pole360=Pole*360)

results.multi<-results.all$multi %>% filter(Measure=='G')



#ncresults<-c('Negative Control','L-Arabinose','Acetoacetic Acid','Phosphono Acetic Acid')


# oldresults<-results.all$results %>%
#   filter(Contrast!="T-C_none" ) %>% #& Metabolite!="Negative Control"
#   group_by(Measure,MetaboliteU) %>%
#   mutate(FDR=p.adjust(p.value,method = "fdr"))%>%
#   group_by(Measure,Contrast) %>%
#   summarise(Total=n(),
#             Ant=sum(logFC>0 & FDR<=0.05,na.rm=TRUE),
#             Syn=sum(logFC<0 & FDR<=0.05,na.rm=TRUE),
#             All=sum(FDR<=0.05,na.rm=TRUE))
# 
# oldresults %>%
#   write_csv('~/Dropbox/Ecoli_old_results.csv')


results %>%
  filter(Metabolite=="Negative Control" & Measure=="G")


results %>%
  filter(Contrast=="T-C" & Measure=="G") %>%
  arrange(desc(logFC)) %>%
  select(Measure,Index,Metabolite,MetaboliteU,logFC,SE,p.value,FDR)  %>%
  filter(Metabolite=="Negative Control") %>%
  View()



write.csv(results,paste0(odir,'/Ecoli_results.csv'),row.names = FALSE)
write.csv(results.cast,paste(odir,'/Ecoli_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/Ecoli_results_sidebyside_full.csv',sep=''),row.names = FALSE)




  
results.celr<-read_csv('Celegans/Summary/Celegans_results_withNGM.csv') %>%
  filter(! Index %in% c('Controls-A1','Controls-B1')) %>%
  select(Index,logFC:logFDR_bin)



celvars<-base::setdiff(colnames(results.celr),'Index')

results.cels<-results.celr %>%
  rename_(.dots=setNames(celvars,paste0('Cel_',celvars)) )


#Main data
results.castcomb<-results.castfull %>%
  left_join(results.cels) %>%
  arrange(desc(`T-C_logFC`))

#View(results.castcomb)

#With growth rate
results.castcomb.c<-results.castfull.c %>%
  left_join(results.cels)


#View(results.castcomb)
write.csv(results.castcomb,paste(odir,'/Combined_results_sidebyside_full.csv',sep=''),row.names = FALSE)



metorder<-as.character(results.castcomb$MetaboliteU)

data.nc %>%
  filter(Measure=='G') %>%
  mutate(MetaboliteU=factor(MetaboliteU,levels=metorder)) %>%
  ggplot(aes(x=MetaboliteU,y=Value_norm,color=Type))+
  geom_hline(yintercept = 0,color='red',alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity") +
  geom_point()+
  ylab('Growth logFC vs NGM')+
  xlab('Metabolite')+
  coord_flip()


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ecoli_logFC_tall_raw.pdf",sep=''),
             width=9,height=50)





#Hasn't been updated from here on


#Generate tables for KEGG enrichment


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







#Read KEGG enrichment results
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






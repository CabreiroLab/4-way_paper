library(tidyverse)
#devtools::install_github("PNorvaisas/PFun")
library(PFun)
library(readxl)


#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau




cwd<-"~/Dropbox/Projects/2015-Metformin/Proteomics/"
setwd(cwd)
odir<-'Results'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")





#Annotation of spots
annot<-read_xlsx('Data Filipe - Metformin Full.xlsx',sheet='Combined with previous analysis') %>%
  mutate(Description=NULL,
         Uniprot=as.character(Uniprot),
         Uniprot=ifelse(Uniprot=='B1XGK9','P77581',Uniprot)) %>%
  rename(Spot=Spot_No,Description=Protein_ID) %>%
  filter(!Uniprot %in% c('NA','/'))



map<-read_delim('UniProt_gene_name.txt',delim = '\t') %>%
  select(UniProtID=From,Gene=To)



#Check for duplicates
table(duplicated(map[,c('Gene','UniProtID')]))



data<-read_xlsx('Filipe POI analysis.xlsx',sheet='Data')

ugroups<-unique(as.character(data$Group))


data.annot<-colnames(data)[1:8]
data.annot

spots<-setdiff(colnames(data),data.annot)
spots

data.m<-data %>%
  gather(Spot,log10SA,spots) %>%
  mutate(log10SA=ifelse(log10SA=='NaN',NA,log10SA)) %>%
  select(-c(Cond_1,Cond_2))



norm<-read_xlsx('Filipe POI analysis.xlsx',sheet='Volumes')
norm[,1:10]

norm$Sample<-data$Sample



norm.m<-norm %>%
  gather(Spot,Volume,spots) %>%
  mutate(Volume=ifelse(Volume=='NaN',NA,Volume)) %>%
  select(-c(`Cond 1`,`Cond 2`))


alldata<-data.m %>%
  left_join(norm.m) %>%
  left_join(annot) %>%
  filter(!Uniprot %in% c(NA,'NA','/')) %>%
  mutate(log2SA=log10SA/log10(2)) %>%
  group_by(Group,Uniprot) %>%
  mutate(Sum_Uni=sum(Volume,na.rm=TRUE),
         Mean_Uni=mean(Volume,na.rm=TRUE),
         log2SA_Uni=mean(log2SA,na.rm=TRUE)) %>%
  group_by(Group,Spot) %>%
  mutate(Sum_Spot=sum(Volume,na.rm=TRUE),
         Mean_Spot=mean(Volume,na.rm=TRUE),
         log2SA_Spot=mean(log2SA,na.rm=TRUE)) %>%
  ungroup %>%
  mutate(Sum=ifelse(is.na(Sum_Uni),Sum_Spot,Sum_Uni),
         Mean=ifelse(is.na(Mean_Uni),Mean_Spot,Mean_Uni),
         Weight=Volume/Sum,
         Coef<-Volume/Mean,
         log2SA_F=ifelse(is.na(log2SA),log2SA_Uni,log2SA),
         log2SA_F=ifelse(is.na(log2SA_F),log2SA_Spot,log2SA_F),
         Metformin_mM=as.factor(Metformin_mM)) %>%
  left_join(map,c('Uniprot'='UniProtID')) %>%
  filter(!Uniprot %in% c(NA,'NA','/'))
  
  
  

# norm.uni<-ddply(alldata.ac,.(Group,Uniprot),summarise,Sum_Uni=sum(Volume,na.rm=TRUE),Mean_Uni=mean(Volume,na.rm=TRUE))
# norm.spot<-ddply(alldata.a,.(Group,Spot),summarise,Sum_Spot=sum(Volume,na.rm=TRUE),Mean_Spot=mean(Volume,na.rm=TRUE))
# alldata.au<-merge(alldata.a,norm.uni,by=c('Group','Uniprot'),all.x=TRUE)
# alldata.as<-merge(alldata.au,norm.spot,by=c('Group','Spot'),all.x=TRUE)
# 
# alldata.as$Sum<-ifelse(is.na(alldata.as$Sum_Uni),alldata.as$Sum_Spot,alldata.as$Sum_Uni)
# alldata.as$Mean<-ifelse(is.na(alldata.as$Mean_Uni),alldata.as$Mean_Spot,alldata.as$Mean_Uni)

# alldata.as$Weight<-alldata.as$Volume/alldata.as$Sum
# alldata.as$Coef<-alldata.as$Volume/alldata.as$Mean

#alldata.as$log2SA<-alldata.as$log10SA/log10(2)
#alldata.as$log2SA_W<-alldata.as$log2SA*alldata.as$Weight

#alldata.asc<-subset(alldata.as,!Uniprot %in% c(NA,'NA','/'))
#table(alldata.asc$Uniprot)

#Fill data for PCA and heatmaps
# data.uni<-ddply(alldata.asc,.(Group,Uniprot),summarise,log2SA_Uni=mean(log2SA,na.rm=TRUE))
# data.spot<-ddply(alldata.as,.(Group,Spot),summarise,log2SA_Spot=mean(log2SA,na.rm=TRUE))

# alldata.fu<-merge(alldata.as,data.uni,by=c('Group','Uniprot'),all.x=TRUE)
# alldata.fs<-merge(alldata.fu,data.spot,by=c('Group','Spot'),all.x=TRUE)
# head(alldata.fs)
# alldata.fs$log2SA_F<-ifelse(is.na(alldata.fs$log2SA),alldata.fs$log2SA_Uni,alldata.fs$log2SA)
# alldata.fs$log2SA_F<-ifelse(is.na(alldata.fs$log2SA_F),alldata.fs$log2SA_Spot,alldata.fs$log2SA_F)
# head(alldata.fs)
# dim(alldata.fs)
# alldata.fsa<-merge(alldata.fs,map,by.x='Uniprot',by.y='UniProtID',all.x=TRUE)
# dim(alldata.fsa)


write.csv(alldata,paste(odir,'/All_raw_data.csv',sep=''),row.names = FALSE)


Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- levels(data_ts$Metformin_mM)
Metlab<-'Metformin, mM'

bioinfo<-alldata %>%
  group_by(Group,Sample,Strain,Metformin_mM) %>%
  summarise %>%
  data.frame
rownames(bioinfo)<-bioinfo$Sample



pcaresults<-alldata%>%
  filter(Strain=="OP50-C" & !Sample %in% c('C_T_1','R_T_3')) %>%
  PCAprep("Sample","Spot","log2SA_F",bioinfo)



ellipses<-pcaresults$Ellipses
pcadata<-pcaresults$pcadata

PC1prc<-pcaresults$PC1prc
PC2prc<-pcaresults$PC2prc


pcadata %>%
  ggplot(aes(x=PC1,y=PC2,color=Metformin_mM))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill= ifelse(Metformin_mM==0,as.character(Strain), NA ), color=Metformin_mM ),size=3,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_fill_manual(name = Metlab,values =Metcols,na.value=NA,guide=FALSE)+
  scale_color_manual(name = Metlab,values =Metcols)+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour=Metcols,fill=c("#FF0000",NA))))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position="top")


ggsave(file=paste0(odir,"/PCA_clean.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")






heatcomplete<-met.sel %>%
  group_by(Metabolite) %>%
  mutate(Completeness=sum(!is.na(Conc_log))*100/n() )%>%
  ungroup %>%
  #At least 50% observations per metabolite
  filter(Completeness>50)


cols<-list(Metformin_mM = c('0' = "white", '50' = "Black"))

hinfo<-bioinfo %>%
  select(Metformin_mM)



#Complete
alldata %>%
  filter(Strain=="OP50-C" & !Sample %in% c('C_T_1','R_T_3')) %>%
  HMap('Sample','Spot','log2SA_F',hinfo,cols=cols)

dev.copy2pdf(device=cairo_pdf,file=paste0(odir,"/Heatmap_clean.pdf"),
             width=10,height=10)



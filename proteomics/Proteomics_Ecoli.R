#Data transformation and analysis
library(xlsx)
library(plyr)
library(reshape2)
library(multcomp)
library(contrast)
library(car)

#Plotting and visualisation
library(gplots)
library(heatmap3)
library(ggplot2)
library(ggbiplot)
library(ggrepel)
library(RColorBrewer)
library(plot3D)

library(gage)


library(gtools)

library(ellipse)

library(Unicode)

library(ggthemes)


#Pathways
library(org.EcK12.eg.db)
library(AnnotationDbi)
library(pathview)
library(limma)

library(splitstackshape)
library(qdap)


#source("https://bioconductor.org/biocLite.R")
#biocLite("gage")


#devtools::install_github("PNorvaisas/PFun")
library(PFun)
#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau




cwd<-"~/Dropbox/Projects/2015-Metformin/Proteomics/"
setwd(cwd)
odir<-'Results'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#load('Proteomics.RData')

getinfo<-function(cof) {
  df<-data.frame(cof)
  df$Comparisons<-rownames(df)
  return(df)
}

getellipse<-function(x,y,sc=1) {
  as.data.frame(ellipse( cor(x, y),
                         scale=c(sd(x)*sc,sd(y)*sc),
                         centre=c( mean(x),mean(y)) ))
}

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}


MinMeanSDMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}


#Annotation of spots
annot.raw<-read.xlsx2('Data Filipe - Metformin Full.xlsx',sheetName='Combined with previous analysis',
                      header=TRUE,endRow=164,check.names = FALSE,
                      colClasses=c(rep("character",3),
                                   rep("numeric", 9),
                                   rep("character", 3)))
annot.raw$Description<-NULL
annot.raw<-rename(annot.raw,c('Spot_No'='Spot','Protein_ID'='Description'))

annot<-annot.raw[,c('Spot','Uniprot','Description')]
annot<-subset(annot,!Uniprot %in% c('NA','/'))

annot$Uniprot<-as.character(annot$Uniprot)
annot$Uniprot<-ifelse(annot$Uniprot=='B1XGK9','P77581',annot$Uniprot)



map<-read.csv('UniProt_gene_name.txt',sep='\t')
colnames(map)<-c('UniProtID','Gene')


map<-subset(map,UniProtID!='B1XGK9')


#Check for duplicates
table(duplicated(map[,c('Gene','UniProtID')]))


data<-read.xlsx2('Filipe POI analysis.xlsx',sheetName='Data',
                    header=TRUE,endRow=19,check.names = FALSE,
                    colClasses=c(rep("character",8), rep("numeric", 313)))
data[,1:10]



ugroups<-unique(as.character(data$Group))


data.annot<-colnames(data)[1:8]
spots<-setdiff(colnames(data),data.annot)


data.m<-melt(data,measure.vars = spots,variable.name = 'Spot',value.name = 'log10SA',na.rm = FALSE)
data.m$log10SA<-ifelse(data.m$log10SA=='NaN',NA,data.m$log10SA)


norm<-read.xlsx2('Filipe POI analysis.xlsx',sheetName='Volumes',
                    header=TRUE,endRow=19,check.names = FALSE,
                    colClasses=c(rep("character",5), rep("numeric", 313)))

norm$Sample<-data$Sample

norm.m<-melt(norm,id.vars = c('Sample','Group'),measure.vars = spots,variable.name = 'Spot',value.name = 'Volume',na.rm = FALSE)
norm.m$Volume<-ifelse(norm.m$Volume=='NaN',NA,norm.m$Volume)


alldata<-merge(data.m,norm.m,by=c('Sample','Group','Spot'),all=TRUE)

alldata.a<-merge(alldata,annot,by='Spot',all.x = TRUE)
alldata.ac<-subset(alldata.a,!Uniprot %in% c(NA,'NA','/'))


norm.uni<-ddply(alldata.ac,.(Group,Uniprot),summarise,Sum_Uni=sum(Volume,na.rm=TRUE),Mean_Uni=mean(Volume,na.rm=TRUE))
norm.spot<-ddply(alldata.a,.(Group,Spot),summarise,Sum_Spot=sum(Volume,na.rm=TRUE),Mean_Spot=mean(Volume,na.rm=TRUE))


alldata.au<-merge(alldata.a,norm.uni,by=c('Group','Uniprot'),all.x=TRUE)
alldata.as<-merge(alldata.au,norm.spot,by=c('Group','Spot'),all.x=TRUE)

alldata.as$Sum<-ifelse(is.na(alldata.as$Sum_Uni),alldata.as$Sum_Spot,alldata.as$Sum_Uni)

alldata.as$Mean<-ifelse(is.na(alldata.as$Mean_Uni),alldata.as$Mean_Spot,alldata.as$Mean_Uni)


alldata.as$Weight<-alldata.as$Volume/alldata.as$Sum
alldata.as$Coef<-alldata.as$Volume/alldata.as$Mean

alldata.as$log2SA<-alldata.as$log10SA/log10(2)

alldata.asc<-subset(alldata.as,!Uniprot %in% c(NA,'NA','/'))

#Fill data for PCA and heatmaps
data.uni<-ddply(alldata.asc,.(Group,Uniprot),summarise,log2SA_Uni=mean(log2SA,na.rm=TRUE))
data.spot<-ddply(alldata.as,.(Group,Spot),summarise,log2SA_Spot=mean(log2SA,na.rm=TRUE))


alldata.fu<-merge(alldata.as,data.uni,by=c('Group','Uniprot'),all.x=TRUE)
alldata.fs<-merge(alldata.fu,data.spot,by=c('Group','Spot'),all.x=TRUE)

alldata.fs$log2SA_F<-ifelse(is.na(alldata.fs$log2SA),alldata.fs$log2SA_Uni,alldata.fs$log2SA)
alldata.fs$log2SA_F<-ifelse(is.na(alldata.fs$log2SA_F),alldata.fs$log2SA_Spot,alldata.fs$log2SA_F)

alldata.fsa<-merge(alldata.fs,map,by.x='Uniprot',by.y='UniProtID',all.x=TRUE)



write.csv(alldata.fsa,paste(odir,'/All_raw_data.csv',sep=''),row.names = FALSE)

prot.c<-subset(alldata.fsa,!Uniprot %in% c(NA,'NA','/')) 

prot.c$ReplicateUniq<-makereplicates(prot.c[,c("Uniprot","Gene","Group")])


#Clean data

prot.c %>%
  group_by(Sample,Strain,Metformin_mM) %>%
  summarise


pcaresults<-prot.c %>%
  PCAprep("Sample","Spot","")



metslm.f<-dcast(subset(prot.c,Strain=="OP50-C"),Sample+Group+Strain+Metformin_mM~Spot,value.var = 'log2SA_F')
rownames(metslm.f)<-metslm.f$Sample


spots.clean<-unique(as.character(prot.c$Spot))


data.c<-data[,spots.clean]


miss.cols<-apply(data.c, 2, function(x) any(is.na(x)))
miss.rows<-apply(data.c, 1, function(x) any(is.na(x)))


miss.cols.count<-apply(data.c, 2, function(x) sum(is.na(x)))
miss.rows.count<-apply(data.c, 1, function(x) sum(is.na(x)))


#Missings rows
cbind(data.c[,1:8],miss.rows.count)
#Missing columns
miss.cols.count


miss.cols.clean<-apply(metslm.f[,spots.clean], 2, function(x) any(is.na(x)))
miss.cols.clean
table(miss.cols.clean)


samples<-as.character(data$Sample)
samples.missing<-as.character(data$Sample)[miss.rows.count>30]

samples.clean<-setdiff(samples,c('C_T_1','R_T_3'))#,'R_C_1','R_C_5'

ancols<-setdiff(colnames(metslm.f),spots)

metslm<-metslm.f[samples.clean,c(ancols,spots.clean)]

#Find compounds with missing values
miss.cols<-apply(metslm, 2, function(x) any(is.na(x)))
miss.rows<-apply(metslm, 1, function(x) any(is.na(x)))

missing.cols<-names(miss.cols[miss.cols==TRUE])
missing.rows<-rownames(metslm)[miss.rows==TRUE]

pca.group<-metslm$Sample

pca.dat<-metslm[,spots.clean]


ir.pca <- prcomp(pca.dat,
                 center = TRUE,
                 scale. = TRUE) 
plot(ir.pca,type='l')



pcadata<-data.frame(ir.pca$x)
pcadata[,c('Sample','Strain','Group','Metformin_mM')]<-metslm[,c('Sample','Strain','Group','Metformin_mM')]

pcaresult<-summary(ir.pca)$importance
PC1prc<-round(pcaresult['Proportion of Variance',][[1]]*100,0)
PC2prc<-round(pcaresult['Proportion of Variance',][[2]]*100,0)

ellipses<-ddply(pcadata,.(Strain,Group,Metformin_mM), summarise, x=getellipse(PC1,PC2,1)$x,y=getellipse(PC1,PC2,0.75)$y ) 


pcadata %>%
  ggplot(aes(x=PC1,y=PC2,colour=Metformin_mM))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=5,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  geom_text(aes(label=Sample))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_clean.pdf",sep=''),
             width=12,height=9)


#Heatmap
prot.cl<-subset(prot.c,Sample %in% samples.clean)

heatshape<-dcast(prot.cl,Spot~Group+Replicate,value.var = 'log2SA_F')

rownames(heatshape)<-heatshape$Spot
heatshape$Spot<-NULL

groups<- substr(colnames(heatshape),1,nchar(colnames(heatshape))-2)
str.tr<- matrix(unlist(strsplit(groups,'_')),2,byrow=FALSE)
strains<-str.tr[1,]
treat<-str.tr[2,]

strains.fac<-factor(strains)
treat.fac<-factor(treat)
#ColSideColors<-rainbow(length(unique(annData)))[as.numeric(annData)]

strains.col<-rainbow(length(unique(strains.fac)))[as.numeric(strains.fac)]
treat.col<-ifelse(treat.fac=='C','white','black')
#grayscale(length(unique(treat.fac)))[as.numeric(treat.fac)]

legtxt<-as.character(unique(strains.fac))
legcol<-rainbow(length(unique(strains.fac)))[as.numeric(unique(strains.fac))]


clab<-cbind(treat.col,strains.col)
colnames(clab)<-c('Treatment','Strain/Mutant')


hmap<-heatmap3(as.matrix(heatshape),
               scale = 'row',
               #hclustfun = function(x) hclust(dist(x,method="manhattan"),method="complete"),
               #distfun = function(x) dist(x,method="manhattan"),
               labRow=FALSE,
               #key=TRUE,
               #ColSideLabs = 'Group',
               ColSideColors = clab)
legend('topleft',legend=legtxt,fill=legcol, border=FALSE, bty="n",title='Strain/Mutant')
legend('left',legend=c('Treatment','Control'),fill=c('black','white'), border=TRUE, bty="n",title='Drug')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_clean.pdf",sep=''),
             width=12,height=9, useDingbats=FALSE)



#Use cleaned data

prot.clean<-subset(prot.c,Sample %in% samples.clean) 

ggplot(prot.clean,aes(x=log2SA,fill=Group))+
  geom_histogram(position='identity',alpha=0.5)


ggplot(prot.clean,aes(x=Gene,y=log2SA,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point(aes(size=Coef))+
  scale_size_continuous(breaks=seq(0,2,by=0.5))+
  geom_text(aes(label=Replicate),color='black',size=2)+
  labs(size='Weight')+
  scale_y_continuous(breaks=seq(-20,20,by=1) )+
  ylab('log 2 Standart abundance')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Strain,labeller = labeller(Strain = label_both),ncol = 5)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_log2SA_by_treatment.pdf",sep=''),
             width=40,height=10, useDingbats=FALSE)


#Check

ggplot(prot.clean,aes(x=Strain,y=log2SA,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point(aes(size=Coef))+
  scale_y_continuous(breaks=seq(-20,20,by=1))+
  geom_text(aes(label=Replicate),color='black',size=2)+
  ylab('log 2 Standart abundance')+
  scale_size_continuous(breaks=seq(0,2,by=0.5))+
  labs(size='Weight')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Gene,ncol = 5)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_log2SA_by_Gene.pdf",sep=''),
             width=15,height=80, useDingbats=FALSE)



#Summarise manually

#Use weighted means and SDs
summary<-ddply(prot.clean,.(Uniprot,Gene,Strain,Metformin_mM,Group),summarise,Mean=wtd.mean(log2SA,Weight,na.rm=TRUE,normwt=TRUE),SD=sqrt(wtd.var(log2SA,Weight,na.rm=TRUE,normwt = TRUE)) )
sum.m<-melt(summary,measure.vars = c('Mean','SD'),variable.name = 'Stat',value.name = 'Value')

sum.m$Index<-paste(sum.m$Group,sum.m$Gene,sep=' ')



sum.c<-dcast(sum.m,Uniprot+Gene+Strain+Metformin_mM+Group+Index~Stat,value.var = 'Value',drop = TRUE)

sum.c$VarPrc<-(2^(sum.c$SD)-1)*100

#Find outliers in the most variable groups
indxord<-sum.c[order(sum.c$SD),'Index']

sum.c$Index<-factor(sum.c$Index,levels=indxord,labels=indxord)


sum.c$Index[duplicated(sum.c$Index)]


ggplot(sum.c,aes(x=Index,y=SD))+
  geom_point()+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_log2SA_SD.pdf",sep=''),
             width=5,height=60, useDingbats=FALSE)


ggplot(sum.c,aes(x=Index,y=VarPrc))+
  geom_point()+
  scale_y_continuous(breaks = seq(0,150,by=25))+
  ylab('In-group variation, %')+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=5,height=60, useDingbats=FALSE)



#Linear modelling

lmshape<-dcast(prot.clean,Group+ReplicateUniq~Gene,value.var = 'log2SA',drop=TRUE)
lmshape.w<-dcast(prot.clean,Group+ReplicateUniq~Gene,value.var = 'Volume',drop=TRUE)


genes<-setdiff(colnames(lmshape),c('Sample','Group','ReplicateUniq'))

sel.groups<-as.character(unique(lmshape$Group))



lmshape.ws<-ddply(lmshape.w,.(Group),function(x) apply(x[,genes],2,sum,na.rm=TRUE))
lmshape.ws


lmshape.wn<-lmshape.w
lmshape.wn[,genes]<-lmshape.w[,genes]/lmshape.ws[match(lmshape.w$Group,lmshape.ws$Group),genes]



#Test if worked
ddply(lmshape.wn,.(Group),function(x) apply(x[,genes],2,sum,na.rm=TRUE))


contrasts<-read.contrasts('!Contrasts_Ecoli_proteomics.xlsx','Contrasts_values',sel.groups)
contrasts.table<-contrasts$Contrasts.table
contr.matrix<-contrasts$Contrasts.matrix

strainlist<-c('OP50','OP50-MR')


contrasts.table$Strain<-factor(contrasts.table$Strain,levels=strainlist,labels=strainlist) #
allresults<-hypothesise(lmshape,genes,contr.matrix,formula="0+Group",weights = lmshape.wn)

results<-merge(allresults$All,contrasts.table[,c('Contrast','Description','Contrast_type','Strain')],by='Contrast',all.x=TRUE)

results<-rename(results,c("Variable"="Gene"))

results.a<-merge(results,map,by.x='Gene',by.y='Gene',all.x=TRUE)
results.a<-rename(results.a,c("UniProtID"="Uniprot"))



egSYM<-toTable(org.EcK12.egSYMBOL2EG)
#Missing Transcript to EntrezGene
dim(subset(egSYM,is.na(gene_id) | is.na(symbol)))
mat<-match(results.a$Gene,egSYM$symbol)
results.a$EntrezGene<-egSYM$gene_id[mat]



egALIAS<-toTable(org.EcK12.egALIAS2EG)
dim(subset(egALIAS,is.na(gene_id) | is.na(alias_symbol)))
mat<-match(results.a$Gene,egALIAS$alias_symbol)
results.a$EntrezGene<-ifelse(is.na(results.a$EntrezGene),egALIAS$gene_id[mat],results.a$EntrezGene)


results.m<-melt(results.a,measure.vars = c('logFC','p.value','SE','t.value','PE','NE','FDR','logFDR'),variable.name = 'Stat',value.name = 'Value')

results.castfull<-dcast(results.m,Gene+Uniprot+EntrezGene~Contrast+Stat,value.var = 'Value')
results.cast<-dcast(subset(results.m,Stat %in% c('logFC','FDR')),Gene+Uniprot+EntrezGene~Contrast+Stat,value.var = 'Value')


results.exp<-results.a[,c('Contrast','Description','Contrast_type','Strain','Gene','Uniprot','logFC','FDR','p.value','SE','PE','NE','t.value')]


write.csv(results.exp,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)


write.table(results.cast,paste(odir,'/All_results_sidebyside.tsv',sep=''),
            row.names = FALSE,col.names=FALSE,sep='\t',quote = FALSE)

contrs<-unique(as.character(results$Contrast))
stats<-c('logFC','FDR','PE','NE')

contcombs<-apply(expand.grid(contrs, stats), 1, paste, collapse="_")
results.cm<-subset(results.a,Contrast=='C_T-C_C')
results.jT<-merge(results.cm,subset(results.a,Contrast!='C_T-C_C'),by=c('Gene','Uniprot'),suffixes = c('_C',''),all.y=TRUE)

#Plotting starts

erralpha<-0.6
errcolor<-'grey80'
lblsize<-2
amp<- 4
cbrks<-seq(-amp,amp,by=1)
#gradcols<-c('black','purple','purple')
maincomp<-'Interaction strength'


gradcols<-c('blue4','blue','gray80','red','red4')

ggplot(subset(results.jT,Contrast_type=='Treatment'),aes(x=logFC_C,y=logFC,color=logFC-logFC_C))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=NE_C,xmax=PE_C),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=NE,ymax=PE),alpha=erralpha,color=errcolor,width=0)+
  geom_point()+
  scale_x_continuous(breaks=seq(-10,10,by=0.5))+
  scale_y_continuous(breaks=seq(-10,10,by=0.5))+
  geom_text_repel(aes(label=Gene),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  xlab('Metformin effect on OP50')+
  ylab('Metformin effect on other strain')+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  facet_wrap(~Strain)+
  theme(panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_treatment.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)



#unique(results.jT$Contrast_type)

ggplot(subset(results.jT,Contrast_type=='Bacterial mutant'),aes(x=logFC_C,y=logFC,color=logFC-logFC_C))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=NE_C,xmax=PE_C),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=NE,ymax=PE),alpha=erralpha,color=errcolor,width=0)+
  geom_point()+
  scale_x_continuous(breaks=seq(-10,10,by=0.5))+
  scale_y_continuous(breaks=seq(-10,10,by=0.5))+
  geom_text_repel(aes(label=Gene),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  xlab('Metformin effect on OP50')+
  ylab('Strain difference')+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  facet_wrap(~Strain)+
  theme(panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_mutant.pdf',sep = ''),
             width=20,height=6,useDingbats=FALSE)




#Volcano plots
ggplot(subset(results,Contrast_type=="Treatment"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+
  #ylim(0,4)+
  ggtitle('Metformin treatment effect in different strains')+
  geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Gene),'')),size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_treatment.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)


ggplot(subset(results,Contrast_type %in% c("Bacterial mutant","Worm mutant","Bacterial strain") ),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+
  #ylim(0,15)+
  ggtitle('Mutant difference vs OP50')+
  geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Gene),'')),size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_mutant.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)




ggplot(subset(results,Contrast_type=="Interaction"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+
  #labs(size='Average metabolite\nconcentration, nmol')+
  #ylim(0,10)+
  ggtitle('Interaction between treatment and strain in comparison to OP50')+
  geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Gene),'')),size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_interaction.pdf',sep = ''),
             width=10,height=6,useDingbats=FALSE)


#Heatmap for summary

comparisons<-c("Treatment effect on OP50","Treatment effect on OP50-MR",
               "Mutant difference for OP50-MR","Interaction with treatment for OP50-MR")


TG<-read.table('Pathway_details.txt',sep='\t',header = TRUE)

heatsum<-dcast(results.a,Gene~Description,value.var = 'logFC')
rownames(heatsum)<-heatsum$Gene
heatsum$Gene<-NULL


amp<-2

minv<- -amp
maxv<- amp

nstep<-maxv-minv

nstep<-8


clrbrks<-seq(-amp,amp,by=1)

brks<-seq(minv,maxv,by=(maxv-minv)/(nstep))
bgg <- colorRampPalette(c("blue", "gray90", "red"))(n = nstep)

clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)


reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }

hm<-heatmap3(as.matrix(heatsum),key=TRUE,Colv=FALSE,trace='none',col=bgg,
             xlab='Comparison',Rowv=TRUE,breaks = brks,dendrogram="row",scale="none")

ordmet<-rownames(heatsum[hm$rowInd,])


if (length(ordmet)!=length(unique(ordmet))){
  print("Non unique metabolites!")
}



#By alphabet

results.a$Gene<-factor(results.a$Gene,levels=ordmet,labels=ordmet)
results.a$Description<-factor(results.a$Description,levels=comparisons,labels=comparisons)
results.a$FDRstars<-stars.pval(results.a$FDR)

pltsel<-results.a

ggplot(subset(pltsel,Description %in% comparisons),aes(x=Description,y=Gene ))+
  geom_tile(aes(fill=logFC))+
  #geom_point(aes(size=FDRc,colour=logFC),alpha=0.9)+
  theme_minimal()+
  geom_text(aes(label=as.character(FDRstars)))+
  #scale_size_discrete(range = c(2,4))+#,breaks=brks
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  #scale_fill_gradient2(low = "purple", mid = "gray", high = "red", midpoint = 0, breaks = clrbrks)+
  xlab("Comparison")+
  theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_heatmap_amp2.pdf',sep = ''),
             width=6,height=25,useDingbats=FALSE)



colnames(results.castfull)
longevity<-subset(results.castfull,`C_T-C_C_FDR`<0.05 & `MR_T-MR_C_FDR`>0.05 & `MR_C-C_C_FDR`>0.05)$Gene

amp<-3

ggplot(subset(results.a,Description %in% comparisons & Gene %in% longevity),aes(x=Description,y=Gene))+
  geom_tile(aes(fill=logFC))+
  #geom_point(aes(size=FDRc,colour=logFC),alpha=0.9)+
  theme_minimal()+
  geom_text(aes(label=as.character(FDRstars)))+
  #scale_size_discrete(range = c(2,4))+#,breaks=brks
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  #scale_fill_gradient2(low = "purple", mid = "gray", high = "red", midpoint = 0, breaks = clrbrks)+
  xlab("Comparison")+
  theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

dev.copy2pdf(device=cairo_pdf,file=paste0(odir,'/Comparison_heatmap_longevity.pdf'),
             width=6,height=12,useDingbats=FALSE)


#Ecocyc enrichment
ecyc<-read.table('../Annotations/Ecoli/EcoCyc_Patwhays.tsv',header=TRUE,sep='\t')
ecyc<-data.frame(splitstackshape::cSplit(ecyc,'Link',sep='='))
ecyc$ID<-ecyc$Link_3
ecyc<-ecyc[,-grep('Link',colnames(ecyc))]


ecyc.m<-data.frame(splitstackshape::cSplit(ecyc,'Genes',sep=';',direction = 'long'))
ecyc.m<-rename(ecyc.m,c('Genes'='Gene'))

ece<-merge(ecyc.m,prota.uc,by='Gene',all.x=TRUE)


idvariables<-c('Pathway','ID','Gene')
selstats<-c('logFC','p.value','FDR')

ece.m<-enrichment.melt(ece,idvariables,selstats)

ece.en<-enrichment(ece.m,terms = c('ID','Pathway'),IDs = 'Gene',comparisons = 'Contrast',change = 'logFC',sign = 'FDR')


write.csv(ece.en,paste0(odir,'/EcoCyc_enrichment.csv'))



#KEGG plotting

gp.KEGG<-subset(getGeneKEGGLinks('eco',convert = TRUE),!is.na(GeneID) & ! is.na(PathwayID) & PathwayID!='path:cel01100')


keggxml<-'~/Dropbox/Projects/2015-Metformin/Annotations/Ecoli/KEGG_pathways/'

gdata<-subset(results.cast,! is.na(EntrezGene))
nrow(gdata)
dupl<-duplicated(gdata[,'EntrezGene'])
print('Entrez Gene duplicates')
print(table(dupl))
gdataf<-gdata[!dupl,]
rownames(gdataf)<-gdataf[,'EntrezGene']


allpaths<-c('00020','00010','00030','01130','00330')



print(paste('Total KEGG pathways to plot:',length(allpaths)))


keggdir<-paste(odir,'/KEGG',sep='')
dir.create(keggdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

pathcomp<-c('C_T-C_C_logFC')#comp.list[[comp.gr]]
outsuffx<-'CT-CC'

setwd(keggdir)
pv.out <- pathview(gene.data = gdataf[,pathcomp,drop=FALSE],
                   pathway.id = allpaths,
                   species = "eco",
                   out.suffix = outsuffx,
                   kegg.dir = keggxml,
                   same.layer=FALSE,kegg.native = T,map.symbol=FALSE,
                   map.null=FALSE,new.signature=FALSE, plot.col.key=TRUE,
                   res=300,
                   min.nnodes = 0,
                   limit=list(gene=2,cpd=1),node.sum = "mean",
                   low = list(gene = "blue", cpd = "green"),
                   mid=list(gene = "gray", cpd = "gray"),
                   high = list(gene = "red", cpd ="yellow"))
setwd(cwd)




#DAVID enrichment

D.all<-read.table('Results/Enrichment/DAVID_CM-CC.txt',sep='\t',header=TRUE)
D.all$Type<-'All'
D.up<-read.table('Results/Enrichment/DAVID_CM-CC_UP.txt',sep='\t',header=TRUE)
D.up$Type<-'Up'
D.down<-read.table('Results/Enrichment/DAVID_CM-CC_DOWN.txt',sep='\t',header=TRUE)
D.down$Type<-'Down'

D.comb<-rbind(D.all,D.up,D.down)

#write.csv(D.comb,paste0(odir,'/DAVID_enrichment.csv'))



#Enrichment combined
#TF enrichment
TFsum<-read.csv('Results/TF_enrichment.csv')

#TF.en<-TFsum[TFsum$TvC_FDR<0.05,c('TF','TvC_p','TvC_FDR')]
TF.en<-TFsum[,c('TF','TvC_p','TvC_FDR')]
TF.en$Test<-'All'
TF.en$Category<-'Transcription factor'
TF.en<-rename(TF.en,c('TvC_p'='p','TvC_FDR'='FDR','TF'='Term'))
TF.en<-TF.en[order(TF.en$FDR),]


#KEGG enrichment

#FDR<0.05 &
KEGG.en<-subset(D.comb,
                  !Term %in% c('eco01130:Biosynthesis of antibiotics') &
                  !Category %in% c('GOTERM_CC_DIRECT','GOTERM_MF_DIRECT','UP_SEQ_FEATURE'))[,c('Type','Category','Term','PValue','FDR')]
KEGG.en<-rename(KEGG.en,c('Type'='Test','PValue'='p'))

#EcoCyc Enrichment



ece.en<-read.csv('Results/EcoCyc_enrichment.csv')
EC.en<-subset(ece.en,Contrast=='C_T-C_C')[,c('Test','Pathway','p','FDR')] #FDR<0.05 & 
EC.en$Category<-'EcoCyc'
EC.en<-rename(EC.en,c('Pathway'='Term'))


all.en<-rbind(TF.en,EC.en,KEGG.en)[,c('Category','Test','Term','p','FDR')]

write.csv(all.en,paste0(odir,'/Enrichment_all_complete.csv'))


all.en.all<-subset(all.en,Test=='All')
all.en.all$Test<-NULL

all.en.all$p <- format(all.en.all$p, scientific = TRUE,digits = 2)
all.en.all$FDR <- format(all.en.all$FDR, scientific = TRUE,digits = 2)

#acetylation<-as.character(as.data.frame(splitstackshape::cSplit(subset(D.comb,Term=='Acetylation' & Type=='All'),splitCols = 'Genes',sep=',',direction='long'))$Genes)
#subset(ece,ID=='NONOXIPENT-PWY')

library(gridExtra)
library(grid)


grid.table(all.en.all)
dev.copy2pdf(device=cairo_pdf,file=paste0(odir,'/Enrichment_table.pdf'),
             width=10,height=6,useDingbats=FALSE)

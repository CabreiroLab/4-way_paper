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



## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("gage")



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

# theme_update(panel.background = element_rect(colour = "black"),
#              axis.text = element_text(colour = "black"))

remove(enrichmentprep)

devtools::install_github("PNorvaisas/PFun")
library(PFun)


cwd<-"~/Dropbox/Projects/2015-Metformin/Proteomics/"
setwd(cwd)
odir<-'Results'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#load('.RData')

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
dim(annot)

annot$Uniprot<-as.character(annot$Uniprot)



annot$Uniprot<-ifelse(annot$Uniprot=='B1XGK9','P77581',annot$Uniprot)

annot


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
data.annot
spots<-setdiff(colnames(data),data.annot)
spots


data.m<-melt(data,measure.vars = spots,variable.name = 'Spot',value.name = 'log10SA',na.rm = FALSE)
data.m$log10SA<-ifelse(data.m$log10SA=='NaN',NA,data.m$log10SA)
head(data.m)




norm<-read.xlsx2('Filipe POI analysis.xlsx',sheetName='Volumes',
                    header=TRUE,endRow=19,check.names = FALSE,
                    colClasses=c(rep("character",5), rep("numeric", 313)))
norm[,1:10]

norm$Sample<-data$Sample

norm.m<-melt(norm,id.vars = c('Sample','Group'),measure.vars = spots,variable.name = 'Spot',value.name = 'Volume',na.rm = FALSE)
norm.m$Volume<-ifelse(norm.m$Volume=='NaN',NA,norm.m$Volume)

head(norm.m)

alldata<-merge(data.m,norm.m,by=c('Sample','Group','Spot'),all=TRUE)
head(alldata)



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

#alldata.as$log2SA_W<-alldata.as$log2SA*alldata.as$Weight

alldata.asc<-subset(alldata.as,!Uniprot %in% c(NA,'NA','/'))

table(alldata.asc$Uniprot)


#Fill data for PCA and heatmaps
data.uni<-ddply(alldata.asc,.(Group,Uniprot),summarise,log2SA_Uni=mean(log2SA,na.rm=TRUE))
data.spot<-ddply(alldata.as,.(Group,Spot),summarise,log2SA_Spot=mean(log2SA,na.rm=TRUE))


alldata.fu<-merge(alldata.as,data.uni,by=c('Group','Uniprot'),all.x=TRUE)

alldata.fs<-merge(alldata.fu,data.spot,by=c('Group','Spot'),all.x=TRUE)
head(alldata.fs)

alldata.fs$log2SA_F<-ifelse(is.na(alldata.fs$log2SA),alldata.fs$log2SA_Uni,alldata.fs$log2SA)
alldata.fs$log2SA_F<-ifelse(is.na(alldata.fs$log2SA_F),alldata.fs$log2SA_Spot,alldata.fs$log2SA_F)
head(alldata.fs)


dim(alldata.fs)
alldata.fsa<-merge(alldata.fs,map,by.x='Uniprot',by.y='UniProtID',all.x=TRUE)

dim(alldata.fsa)


write.csv(alldata.fsa,paste(odir,'/All_raw_data.csv',sep=''),row.names = FALSE)

head(map)

prot.c<-subset(alldata.fsa,!Uniprot %in% c(NA,'NA','/')) 


prot.c$ReplicateUniq<-makereplicates(prot.c[,c("Uniprot","Gene","Group")])
# 
# prot.c$ReplicateUniq<-1
# alluniq<-FALSE
# while(!alluniq){
#   print(max(prot.c[,'ReplicateUniq']))
#   duplics<-duplicated(prot.c[,c("Uniprot","Gene","Group","ReplicateUniq")])
#   alluniq<-all(duplics==FALSE)
#   if (!alluniq) {
#     prot.c[duplics,'ReplicateUniq']<-max(prot.c[duplics,'ReplicateUniq'])+1
#   }
# }

subset(prot.c,Gene=='astC')

dim(prot.c)

#Clean data


metslm.f<-dcast(prot.c,Sample+Group+Strain+Metformin_mM~Spot,value.var = 'log2SA_F')
rownames(metslm.f)<-metslm.f$Sample


spots.clean<-unique(as.character(prot.c$Spot))

metslm.f[,1:10]


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

#samples.clean<-samples
samples.missing


#samples.clean<-samples
#samples.clean<-setdiff(samples,c('C_T_1','C_T_4','R_T_3',samples.missing))
#samples.clean<-setdiff(samples,c('C_T_1','C_T_4','R_T_3',c('C_T_3','C_C_4')))#


samples.clean<-setdiff(samples,c('C_T_1','R_T_3'))#,'R_C_1','R_C_5'



ancols<-setdiff(colnames(metslm.f),spots)

metslm<-metslm.f[samples.clean,c(ancols,spots.clean)]

dim(metslm)



metslm[,1:10]


#Find compounds with missing values
miss.cols<-apply(metslm, 2, function(x) any(is.na(x)))
miss.rows<-apply(metslm, 1, function(x) any(is.na(x)))

missing.cols<-names(miss.cols[miss.cols==TRUE])
missing.rows<-rownames(metslm)[miss.rows==TRUE]

missing.cols
missing.rows

# 
# metslm[missing.rows,]
# cleancols<-setdiff(cols,missing.cols)
# cleanrows<-setdiff(rownames(metslm),missing.rows)


pca.group<-metslm$Sample

pca.dat<-metslm[,spots.clean]


ir.pca <- prcomp(pca.dat,
                 center = TRUE,
                 scale. = TRUE) 
plot(ir.pca,type='l')



pcadata<-data.frame(ir.pca$x)
pcadata[,c('Sample','Strain','Group','Metformin_mM')]<-metslm[,c('Sample','Strain','Group','Metformin_mM')]

#pcadata<-merge(data.frame(ir.pca$x),metslm[,c('Sample','Strain','Metformin_mM')],by = 0)

pcaresult<-summary(ir.pca)$importance
PC1prc<-round(pcaresult['Proportion of Variance',][[1]]*100,0)
PC2prc<-round(pcaresult['Proportion of Variance',][[2]]*100,0)

ellipses<-ddply(pcadata,.(Strain,Group,Metformin_mM), summarise, x=getellipse(PC1,PC2,1)$x,y=getellipse(PC1,PC2,0.75)$y ) 


ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=5,stroke=1,shape=21)+
  scale_linetype_manual("Drug, mM",values=c("0"=1,"50"=2))+
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
head(prot.clean)

heatshape<-dcast(prot.cl,Spot~Group+Replicate,value.var = 'log2SA_F')

rownames(heatshape)<-heatshape$Spot
heatshape$Spot<-NULL

groups<- substr(colnames(heatshape),1,nchar(colnames(heatshape))-2)
str.tr<- matrix(unlist(strsplit(groups,'_')),2,byrow=FALSE)
strains<-str.tr[1,]
treat<-str.tr[2,]


colnames(heatshape)
groups

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
#prot.clean<-subset(alldata.fs,Sample %in% samples.clean)

#prot.clean<-subset(prot.c,Sample %in% samples.clean)

prot.clean<-subset(prot.c,Sample %in% samples.clean) #,samples.missing


ggplot(prot.clean,aes(x=log2SA,fill=Group))+geom_histogram(position='identity',alpha=0.5)


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
#P36683 P37902

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

#,labeller = labeller(Metablite = label_both)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_log2SA_by_Gene.pdf",sep=''),
             width=15,height=80, useDingbats=FALSE)



#Summarise manually

#Use weighted means and SDs
summary<-ddply(prot.clean,.(Uniprot,Gene,Strain,Metformin_mM,Group),summarise,Mean=wtd.mean(log2SA,Weight,na.rm=TRUE,normwt=TRUE),SD=sqrt(wtd.var(log2SA,Weight,na.rm=TRUE,normwt = TRUE)) )
sum.m<-melt(summary,measure.vars = c('Mean','SD'),variable.name = 'Stat',value.name = 'Value')

head(sum.m)

sum.m$Index<-paste(sum.m$Group,sum.m$Gene,sep=' ')



sum.c<-dcast(sum.m,Uniprot+Gene+Strain+Metformin_mM+Group+Index~Stat,value.var = 'Value',drop = TRUE)

sum.c$VarPrc<-(2^(sum.c$SD)-1)*100

head(sum.c)

#Find outliers in most variable groups
indxord<-sum.c[order(sum.c$SD),'Index']

sum.c$Index<-factor(sum.c$Index,levels=indxord,labels=indxord)


sum.c$Index[duplicated(sum.c$Index)]

subset(prot.clean,Gene=='astC')


ggplot(sum.c,aes(x=Index,y=SD))+
  geom_point()+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_log2SA_SD.pdf",sep=''),
             width=5,height=60, useDingbats=FALSE)


max(sum.c$VarPrc,na.rm = TRUE)

ggplot(sum.c,aes(x=Index,y=VarPrc))+
  geom_point()+
  scale_y_continuous(breaks = seq(0,150,by=25))+
  ylab('In-group variation, %')+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=5,height=60, useDingbats=FALSE)


subset(prot.clean,Uniprot=='P76217' & Group=='R_T')




#Linear modelling

lmshape<-dcast(prot.clean,Group+ReplicateUniq~Gene,value.var = 'log2SA',drop=TRUE)
lmshape.w<-dcast(prot.clean,Group+ReplicateUniq~Gene,value.var = 'Volume',drop=TRUE)


genes<-setdiff(colnames(lmshape),c('Sample','Group','ReplicateUniq'))

sel.groups<-as.character(unique(lmshape$Group))



lmshape.ws<-ddply(lmshape.w,.(Group),function(x) apply(x[,genes],2,sum,na.rm=TRUE))
lmshape.ws
#lmshape.ws<-apply(lmshape.w[,genes],2,mean,na.rm=TRUE)


lmshape.wn<-lmshape.w
lmshape.wn[,genes]<-lmshape.w[,genes]/lmshape.ws[match(lmshape.w$Group,lmshape.ws$Group),genes]



#Test if worked
ddply(lmshape.wn,.(Group),function(x) apply(x[,genes],2,sum,na.rm=TRUE))


contrasts<-read.contrasts('!Contrasts_Ecoli_proteomics.xlsx','Contrasts_values',sel.groups)
contrasts.table<-contrasts$Contrasts.table
contr.matrix<-contrasts$Contrasts.matrix

strainlist<-c('OP50','OP50-MR')


contrasts.table$Strain<-factor(contrasts.table$Strain,levels=strainlist,labels=strainlist) #
contr.matrix

dim(lmshape)
dim(lmshape.wn)


allresults<-hypothesise(lmshape,genes,contr.matrix,formula="0+Group",weights = lmshape.wn)


results<-merge(allresults$All,contrasts.table[,c('Contrast','Description','Contrast_type','Strain')],by='Contrast',all.x=TRUE)

results<-rename(results,c("Variable"="Gene"))

head(results)


results.a<-merge(results,map,by.x='Gene',by.y='Gene',all.x=TRUE)
results.a<-rename(results.a,c("UniProtID"="Uniprot"))



egSYM<-toTable(org.EcK12.egSYMBOL2EG)
#Missing Transcript to EntrezGene
dim(subset(egSYM,is.na(gene_id) | is.na(symbol)))
mat<-match(results.a$Gene,egSYM$symbol)
results.a$EntrezGene<-egSYM$gene_id[mat]



egALIAS<-toTable(org.EcK12.egALIAS2EG)
head(egALIAS)
dim(subset(egALIAS,is.na(gene_id) | is.na(alias_symbol)))
mat<-match(results.a$Gene,egALIAS$alias_symbol)
results.a$EntrezGene<-ifelse(is.na(results.a$EntrezGene),egALIAS$gene_id[mat],results.a$EntrezGene)

subset(results.a,is.na(EntrezGene))

#Needs to be added manually
subset(egALIAS,grepl('tal',alias_symbol))





c('galF')
c('ECK2036')



head(results.a)

subset(results.a,Gene=='astC')


results.m<-melt(results.a,measure.vars = c('logFC','p.value','SE','t.value','PE','NE','FDR','logFDR'),variable.name = 'Stat',value.name = 'Value')


head(results.m)

results.castfull<-dcast(results.m,Gene+Uniprot+EntrezGene~Contrast+Stat,value.var = 'Value')
results.cast<-dcast(subset(results.m,Stat %in% c('logFC','FDR')),Gene+Uniprot+EntrezGene~Contrast+Stat,value.var = 'Value')



head(results.castfull)
head(results.cast)

#View(results.cast)


subset(results.cast,Gene=='astC')


subset(results.cast,is.na(EntrezGene))


results.exp<-results.a[,c('Contrast','Description','Contrast_type','Strain','Gene','Uniprot','logFC','FDR','p.value','SE','PE','NE','t.value')]

head(results.exp)


write.csv(results.exp,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)


write.table(results.cast,paste(odir,'/All_results_sidebyside.tsv',sep=''),
            row.names = FALSE,col.names=FALSE,sep='\t',quote = FALSE)


head(results.cast)


dim(subset(results.cast,`MR_C-C_C_FDR`<0.05))
dim(subset(results.cast,`C_T-C_C_FDR`<0.05))


contrs<-unique(as.character(results$Contrast))
stats<-c('logFC','FDR','PE','NE')

contcombs<-apply(expand.grid(contrs, stats), 1, paste, collapse="_")



results.cm<-subset(results.a,Contrast=='C_T-C_C')


results.jT<-merge(results.cm,subset(results.a,Contrast!='C_T-C_C'),by=c('Gene','Uniprot'),suffixes = c('_C',''),all.y=TRUE)




head(results.jT)
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


unique(as.character(results.a$Description))

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
#ordmet<-ordmet[order(ordmet)]


results.a$Gene<-factor(results.a$Gene,levels=ordmet,labels=ordmet)
results.a$Description<-factor(results.a$Description,levels=comparisons,labels=comparisons)
results.a$FDRstars<-stars.pval(results.a$FDR)

#pltsel<-subset(results.a,Gene %in% TG$Gene.name)
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

longevity

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


head(results.castfull)


#Enrichment


annotuni<-read.csv('../Annotations/Ecoli/UniPort_GeneID_20170719.csv',header = TRUE)
head(annotuni)

head(map)



prota.uc<-merge(results.castfull,annotuni[,c('b_ID','Gene_name')],by.x='Gene',by.y='Gene_name',all.x=TRUE)
prota.uc$TvC_Rank<-ifelse(prota.uc$`C_T-C_C_logFC`<0,'-1','3')
prota.uc$RvS_Rank<-ifelse(prota.uc$`MR_C-C_C_logFC`<0,'-1','3')


colnames(prota.uc)
stats<-setdiff(colnames(results.castfull),'Uniprot')



stats.sel<-stats[grepl('logFC',stats)|grepl('p.value',stats)|grepl('FDR',stats)]

prota.uc<-prota.uc[,c('Gene','Uniprot','b_ID','TvC_Rank','RvS_Rank',stats.sel)]

write.csv(prota.uc,paste0(odir,'/Ecoli_proteomics_FBA.csv'),quote = FALSE,row.names = FALSE)


TF<-read.csv('../Annotations/Ecoli/Transcription_Factors/network_tf_gene.txt',comment.char = '#',sep = '\t',header=FALSE)
colnames(TF)<-c('TF','Gene','Effect','Evidence','Support','NAs')
TF$NAs<-NULL

head(TF)


TFa<-merge(TF,prota.uc,by='Gene',all=TRUE)
TFa<-subset(TFa,!is.na(TF))

subset(TFa,Gene=='aceE')
Gene_counts<-data.frame(table(TFa$Gene))
Gene_counts.u<-subset(Gene_counts,Freq==1)

#Is gene under single TF
TFa$Unique<-TFa$Gene %in% Gene_counts.u$Var1
TFa<-TFa[,c('TF','Gene','Effect','Support','Unique',stats.sel)]
head(TFa)

write.csv(TFa,paste0(odir,'/TF_gene.csv'),quote = FALSE)



#TFsel<-subset(TFa,!is.na(`C_T-C_C_logFC`) &  Gene %in% TG$Gene.name)# TF %in% c('Cra','CRP','ArgR','NtrC') &
TFsel<-subset(TFa,!is.na(`C_T-C_C_logFC`) & TF %in% c('Cra','CRP','ArgR','NtrC'))#  &

TFsel
duplicated(TFsel$TF,TFsel$Gene)


TFsel$ReplicateUniq<-makereplicates(TFsel[,c("TF","Gene")])

TFselcast<-dcast(TFsel,Gene+ReplicateUniq~TF,value.var = 'Effect')

TFselcast

write.csv(TFselcast,paste0(odir,'/TF_CraCRPArgRNtrC.csv'),quote = FALSE)


#Enrichment analysis
#TFs

TFa.sel<-TFa
#TFa.sel<-subset(TFa,Unique)

length(TFa$Gene)
length(TFa.sel$Gene)

TFau<-TFa.sel[!duplicated(TFa.sel[,c('Gene','TF')]),]
allgenes<-length(unique(TFau$Gene))

head(TFau)


TFsum<-ddply(TFau,.(TF),summarise,
             TF_total=length(Gene),
             TvC_sig=sum(`C_T-C_C_p.value`<0.05,na.rm = TRUE),
             RvS_sig=sum(`MR_C-C_C_p.value`<0.05,na.rm = TRUE),
             Longevity_sig=sum(`C_T-C_C-(MR_T-MR_C)_p.value`<0.05,na.rm = TRUE))

head(TFsum)



#How many significant in total
TvC_sig_total<-length(unique(subset(TFau,`C_T-C_C_FDR`<0.05)$Gene))
RvS_sig_total<-length(unique(subset(TFau,`MR_C-C_C_FDR`<0.05)$Gene))
Longevity_sig_total<-length(unique(subset(TFau,`C_T-C_C-(MR_T-MR_C)_p.value`<0.05)$Gene))

TFsum$TvC_total<-TvC_sig_total
TFsum$RvS_total<-RvS_sig_total
TFsum$Longevity_total<-Longevity_sig_total

#allgenes
#TvC_sig_total
#RvS_sig_total

TFsum$TvC_p<-phyper(TFsum$TvC_sig-1,TFsum$TvC_total,allgenes-TFsum$TvC_total,TFsum$TF_total,lower.tail =FALSE)
TFsum$RvS_p<-phyper(TFsum$RvS_sig-1,TFsum$RvS_total,allgenes-TFsum$RvS_total,TFsum$TF_total,lower.tail = FALSE)
TFsum$Longevity_p<-phyper(TFsum$Longevity_sig-1,TFsum$Longevity_total,allgenes-TFsum$Longevity_total,TFsum$TF_total,lower.tail = FALSE)


TFsum$TvC_FDR<-p.adjust(TFsum$TvC_p,method='fdr')
TFsum$RvS_FDR<-p.adjust(TFsum$RvS_p,method='fdr')
TFsum$Longevity_FDR<-p.adjust(TFsum$Longevity_p,method='fdr')


TFsum.n<-enrichment(TFau.c,'TF','Gene','Contrast','logFC','FDR')

TFsum<-TFsum[order(TFsum$TvC_p),]


subset(TFsum,TvC_FDR<0.05)

subset(TFsum,TF=='FNR')
TvC_sig/allgenes
30/513



TFsum<-TFsum[,c('TF','TF_total',
                'TvC_total','RvS_total','Longevity_total',
                'TvC_sig','RvS_sig','Longevity_sig',
                'TvC_p','TvC_FDR',
                'RvS_p','RvS_FDR',
                'Longevity_p','Longevity_FDR')]
write.csv(TFsum,paste0(odir,'/TF_enrichment.csv'))





subset(TFsum,TvC_FDR<0.05)


head(TFa)

head(TFsum)


all<-merge(TFa,TFsum,by='TF',all.x=TRUE)
all$Evidence<-NULL

write.csv(all,paste0(odir,'/TF_gene_all.csv'),quote = FALSE)

TFsum[order(TFsum$TF_total,decreasing = TRUE),]




#Ecocyc enrichment
ecyc<-read.table('../Annotations/Ecoli/EcoCyc_Patwhays.tsv',header=TRUE,sep='\t')
ecyc<-data.frame(splitstackshape::cSplit(ecyc,'Link',sep='='))
ecyc$ID<-ecyc$Link_3
ecyc<-ecyc[,-grep('Link',colnames(ecyc))]


ecyc.m<-data.frame(splitstackshape::cSplit(ecyc,'Genes',sep=';',direction = 'long'))
ecyc.m<-rename(ecyc.m,c('Genes'='Gene'))

head(ecyc.m)

ece<-merge(ecyc.m,prota.uc,by='Gene',all.x=TRUE)

head(ece)




idvariables<-c('Pathway','ID','Gene')
selstats<-c('logFC','p.value','FDR')

ece.m<-enrichment.melt(ece,idvariables,selstats)

head(ece.m)

#remove(enrichmentmelt)

ece.en<-enrichment(ece.m,terms = c('ID','Pathway'),IDs = 'Gene',comparisons = 'Contrast',change = 'logFC',sign = 'FDR')


head(ece.en)
subset(ece.en,FDR<0.05)

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


head(gdataf)

allpaths<-c('00020','00010','00030','01130','00330')



print(paste('Total KEGG pathways to plot:',length(allpaths)))


max(gdataf$`C_T-C_C_logFC`)
min(gdataf$`C_T-C_C_logFC`)


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

write.csv(D.comb,paste0(odir,'/DAVID_enrichment.csv'))




#Enrichment combined
TF.en<-TFsum[TFsum$TvC_FDR<0.05,c('TF','TvC_p','TvC_FDR')]
TF.en$Test<-'All'
TF.en$Category<-'Transcription factor'
TF.en<-rename(TF.en,c('TvC_p'='p','TvC_FDR'='FDR','TF'='Term'))
TF.en<-TF.en[order(TF.en$FDR),]

KEGG.en<-subset(D.comb,FDR<0.05 &
                  !Term %in% c('eco01130:Biosynthesis of antibiotics') &
                  !Category %in% c('GOTERM_CC_DIRECT','GOTERM_MF_DIRECT','UP_SEQ_FEATURE'))[,c('Type','Category','Term','PValue','FDR')]
KEGG.en<-rename(KEGG.en,c('Type'='Test','PValue'='p'))

EC.en<-subset(ece.en,FDR<0.05 & Contrast=='C_T-C_C')[,c('Test','Pathway','p','FDR')]
EC.en$Category<-'EcoCyc'
EC.en<-rename(EC.en,c('Pathway'='Term'))


head(TF.en)

head(EC.en)

head(KEGG.en)

all.en<-rbind(TF.en,EC.en,KEGG.en)[,c('Category','Test','Term','p','FDR')]
all.en<-rename(all.en,c('FDR'='q'))


all.en.all<-subset(all.en,Test=='All')
all.en.all$Test<-NULL

all.en.all$p <- format(all.en.all$p, scientific = TRUE,digits = 2)
all.en.all$q <- format(all.en.all$q, scientific = TRUE,digits = 2)

#acetylation<-as.character(as.data.frame(splitstackshape::cSplit(subset(D.comb,Term=='Acetylation' & Type=='All'),splitCols = 'Genes',sep=',',direction='long'))$Genes)
#subset(ece,ID=='NONOXIPENT-PWY')

library(gridExtra)
library(grid)


grid.table(all.en.all)
dev.copy2pdf(device=cairo_pdf,file=paste0(odir,'/Enrichment_table.pdf'),
             width=10,height=6,useDingbats=FALSE)

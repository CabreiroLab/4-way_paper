#Data transformation and analysis
library(xlsx)
library(plyr)
library(reshape2)
library(multcomp)
library(contrast)
library(ellipse)
library(car)

#Plotting and visualisation
library(gplots)
library(heatmap3)
library(ggplot2)
library(ggbiplot)
library(ggrepel)
library(RColorBrewer)
library(plot3D)


library(gtools)


devtools::install_github("PNorvaisas/PFun")
library(PFun)



theme_set(theme_light())

theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))

setwd("~/Dropbox/Projects/2015-Metformin/Metabolomics/")





odir<-'Summary_Ecoli_nucleotides'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



getinfo<-function(cof) {
  df<-data.frame(cof)
  df$Comparisons<-rownames(df)
  return(df)
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

getellipse<-function(x,y,sc=1) {
  as.data.frame(ellipse::ellipse( cor(x, y),
                         scale=c(sd(x)*sc,sd(y)*sc),
                         centre=c( mean(x),mean(y)) ))
}



#Needs to be manual

nuc.raw<-read.xlsx2('Ecoli_nucleotides.xlsx',sheetName = 'Clean',header = TRUE,
                    colClasses = c(rep('character',6),rep('numeric',32)),check.names=FALSE)

metabolites<-colnames(nuc.raw)[7:dim(nuc.raw)[2]]

annot<-read.xlsx2('Nucleotide_annotation.xlsx','Annotation')



norms<-apply(nuc.raw[,metabolites],1,mean)
normalise<-norms/mean(norms)


nuc.norm<-nuc.raw

nuc.norm[,metabolites]<-nuc.norm[,metabolites]/normalise

#Test normalisation
apply(nuc.norm[,metabolites],1,mean)



nuc<-melt(nuc.norm,measure.vars = metabolites,variable.name = 'Metabolite',value.name = 'Conc')


#nuc<-merge(nuc.m,annot,by.x = 'Metabolite',by.y='Compound',all.x = TRUE)
nuc<-subset(nuc,Metabolite!="")

nuc$Conc<-ifelse(nuc$Conc==0,NA,nuc$Conc)
nuc$Metformin_mM<-as.factor(nuc$Metformin)
strains<-c('OP50','OP50-MR')

nuc$Strain<-factor(nuc$Strain,levels=strains,labels=strains)



nuc$Metabolite<-as.factor(nuc$Metabolite)
nuc$Conc<-as.numeric(nuc$Conc) 
nuc$logConc<-log(nuc$Conc,2)

head(nuc)



metslm<-dcast(nuc,Sample+Group+Strain+Metformin_mM+Replicate~Metabolite,value.var = c('logConc'),fill = as.numeric(NA),drop=TRUE)
cols<-setdiff(colnames(metslm),c('Sample','Group','Strain','Glucose','Metformin_mM','Replicate'))
rownames(metslm)<-metslm$Sample

head(metslm)


#Find compounds with missing values
miss.cols<-apply(metslm, 2, function(x) any(is.na(x)))
miss.rows<-apply(metslm, 1, function(x) any(is.na(x)))

missing.cols<-names(miss.cols[miss.cols==TRUE])
missing.rows<-rownames(metslm)[miss.rows==TRUE]

missing.cols
missing.rows

metslm[missing.rows,]

cleancols<-setdiff(cols,missing.cols)
cleanrows<-setdiff(rownames(metslm),missing.rows)



pca.group<-metslm$Sample

pca.dat<-metslm[,metabolites]



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



ellipses<-ddply(pcadata,.(Strain,Group,Metformin_mM), summarise, x=getellipse(PC1,PC2,0.75)$x,y=getellipse(PC1,PC2,0.75)$y ) 


ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=5,stroke=1,shape=21)+
  scale_linetype_manual("Drug, mM",values=c("0"=1,"50"=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  #geom_text(aes(label=Sample))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA.pdf",sep=''),
             width=12,height=9)


allmets<-unique(as.character(nuc$Metabolite))
cleanmets<-setdiff(allmets,missing.cols)





#Heatmap
heatshape<-dcast(nuc,Metabolite~Sample,value.var = 'logConc')

rownames(heatshape)<-heatshape$Metabolite
heatshape$Metabolite<-NULL



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
colnames(clab)<-c('Metformin','Strain')


hmap<-heatmap3(as.matrix(heatshape),
               scale = 'row',
               #hclustfun = function(x) hclust(dist(x,method="manhattan"),method="complete"),
               #distfun = function(x) dist(x,method="manhattan"),
               #labRow=FALSE,
               #key=TRUE,
               #ColSideLabs = 'Group',
               ColSideColors = clab)
legend('topleft',legend=legtxt,fill=legcol, border=FALSE, bty="n",title='Strain')
legend('left',legend=c('Treatment','Control'),fill=c('black','white'), border=TRUE, bty="n",title='Metformin')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap.pdf",sep=''),
             width=12,height=9, useDingbats=FALSE)








ggplot(nuc,aes(x=Metabolite,y=logConc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point()+
  geom_text(aes(label=Replicate),color='black',size=2)+
  scale_y_continuous(breaks=seq(-20,20,by=1) )+
  ylab('log 2 Concentration, nmol')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Strain,labeller = labeller(Strain = label_both),ncol = 5)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_treatment.pdf",sep=''),
             width=25,height=10, useDingbats=FALSE)




ggplot(nuc,aes(x=Strain,y=logConc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point()+
  scale_y_continuous(breaks=seq(-20,20,by=1))+
  geom_text(aes(label=Replicate),color='black',size=2)+
  ylab('log 2 Concentration, nmol')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Metabolite,ncol = 5)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_metabolite.pdf",sep=''),
             width=10,height=20, useDingbats=FALSE)






unique(as.character(met$Group))

#Summarise manually

met.mes<-melt(nuc,measure.vars = c('Conc','logConc'),variable.name = 'Measure',value.name = 'Value')
head(met.mes)

#Concentrations by groups
summary<-ddply(met.mes,.(Metabolite,Strain,Metformin_mM,Group,Measure),summarise,Mean=mean(Value),SD=sd(Value))
sum.m<-melt(summary,measure.vars = c('Mean','SD'),variable.name = 'Stat',value.name = 'Value')


sum.m$Index<-paste(sum.m$Group,sum.m$Metabolite,sep=' ')


sum.c<-dcast(sum.m,Metabolite+Strain+Metformin_mM+Group+Index~Measure+Stat,value.var = 'Value',drop = TRUE)
sum.c$VarPrc<-(2^(sum.c$logConc_SD)-1)*100


#Find outliers in most variable groups
indxord<-sum.c[order(sum.c$logConc_SD),'Index']
sum.c$Index<-factor(sum.c$Index,levels=indxord,labels=indxord)

ggplot(sum.c,aes(x=Index,y=logConc_SD))+
  geom_point()+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_logConc_SD.pdf",sep=''),
             width=5,height=20, useDingbats=FALSE)

ggplot(sum.c,aes(x=Index,y=VarPrc))+
  geom_point()+
  scale_y_continuous(breaks = seq(0,500,by=50))+
  ylab('In-group variation, %')+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=5,height=20, useDingbats=FALSE)




#Absolute concentrations over all samples
met.abs<-ddply(met.mes,.(Metabolite,Measure),summarise,Mean=mean(Value),SD=sd(Value))
met.am<-melt(met.abs,measure.vars = c('Mean','SD'),variable.name = 'Stat',value.name = 'Value')
#Average concentration over all samples
met.c<-dcast(met.am,Metabolite~Measure+Stat,value.var = 'Value',drop = TRUE)

head(met.am)

head(met.c)

sel.groups<-unique(nuc$Group)
sel.groups








#Linear modelling

lmshape<-dcast(nuc,Sample+Group+Replicate~Metabolite,value.var = 'logConc',drop=TRUE)
rownames(lmshape)<-lmshape$Sample
head(lmshape)

lmshape<-subset(lmshape,Group %in% sel.groups)

unique(lmshape$Group)



contrasts<-read.contrasts('!Contrasts_Ecoli_nucleotide_metabolomics.xlsx','Contrasts_values',sel.groups)
contrasts.table<-contrasts$Contrasts.table
contr.matrix<-contrasts$Contrasts.matrix

strainlist<-c('OP50','OP50-MR')


contrasts.table$Strain<-factor(contrasts.table$Strain,levels=strainlist,labels=strainlist) #
contr.matrix


allresults<-hypothesise(lmshape,metabolites,contr.matrix,formula="0+Group")


results<-merge(allresults$All,contrasts.table[,c('Contrast','Description','Contrast_type','Strain')],by='Contrast',all.x=TRUE)

results<-rename(results,c("Variable"="Metabolite"))

results.a<-merge(results,met.c,by='Metabolite')

head(results.a)



results.m<-melt(results.a,measure.vars = c('logFC','p.value','SE','t.value','PE','NE','FDR','logFDR'),variable.name = 'Stat',value.name = 'Value')


head(results.m)
results.castfull<-dcast(results.m,Metabolite+Conc_Mean+Conc_SD+logConc_Mean+logConc_SD~Contrast+Stat,value.var = 'Value')
results.cast<-dcast(subset(results.m,Stat %in% c('logFC','FDR')),Metabolite+Conc_Mean+Conc_SD+logConc_Mean+logConc_SD~Contrast+Stat,value.var = 'Value')



head(results.castfull)
head(results.cast)



results.exp<-results.a[,c('Contrast','Description','Contrast_type','Strain','Metabolite','logFC','FDR','p.value','SE','PE','NE','t.value')]
head(results.exp)


write.csv(results.exp,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)




contrs<-unique(as.character(results$Contrast))
stats<-c('logFC','FDR','PE','NE')

contcombs<-apply(expand.grid(contrs, stats), 1, paste, collapse="_")



results.cm<-subset(results.a,Contrast=='C_T-C_C')


results.jT<-merge(results.cm,subset(results,Contrast!='C_T-C_C'),by='Metabolite',suffixes = c('_C',''),all.y=TRUE)




head(results.jT)
#Plotting starts

erralpha<-0.6
errcolor<-'grey80'
lblsize<-2
amp<- 5
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
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  geom_text_repel(aes(label=Metabolite),
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
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  geom_text_repel(aes(label=Metabolite),
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
             width=12,height=6,useDingbats=FALSE)




#Volcano plots
ggplot(subset(results.a,Contrast_type=="Treatment"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  scale_y_continuous(breaks = seq(0,20,by=1),limits = c(0,6))+
  labs(size='Average metabolite\nconcentration, pmol/ug')+
  ggtitle('Metformin treatment effect in different strains')+
  geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_treatment.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)


ggplot(subset(results.a,Contrast_type %in% c("Bacterial mutant","Worm mutant","Bacterial strain") ),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  scale_y_continuous(breaks = seq(0,20,by=1),limits = c(0,6))+
  labs(size='Average metabolite\nconcentration, nmol')+
  ggtitle('Mutant difference vs OP50 (N2)')+
  #geom_text_repel(aes(label=ifelse(FDR<0.05,as.character(Metabolite),'')),size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_mutant.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)




ggplot(subset(results.a,Contrast_type=="Interaction"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  labs(size='Average metabolite\nconcentration, pmol/ug')+
  scale_y_continuous(breaks = seq(0,20,by=1),limits = c(0,3))+
  ggtitle('Interaction between treatment and strain in comparison to OP50 (N2)')+
  geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_interaction.pdf',sep = ''),
             width=5,height=6,useDingbats=FALSE)





#Heatmap for summary


unique(as.character(results.a$Description))

comparisons<-c("Treatment effect on OP50","Treatment effect on OP50-MR",
               "Mutant difference for OP50-MR","Interaction with treatment for OP50-MR")



heatsum<-dcast(results.a,Metabolite~Description,value.var = 'logFC')
rownames(heatsum)<-heatsum$Metabolite
heatsum$Metabolite<-NULL


amp<-5

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


results.a$Metabolite<-factor(results.a$Metabolite,levels=ordmet,labels=ordmet)
results.a$Description<-factor(results.a$Description,levels=comparisons,labels=comparisons)
results.a$FDRstars<-stars.pval(results.a$FDR)

max(results.a$logFC)
min(results.a$logFC)

amp<-5

ggplot(subset(results.a,Description %in% comparisons),aes(x=Description,y=Metabolite))+
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

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_heatmap.pdf',sep = ''),
             width=6,height=10,useDingbats=FALSE)




longevity<-subset(results.castfull,`C_T-C_C_FDR`<0.05 & `MR_T-MR_C_FDR`>0.05 & `MR_C-C_C_FDR`>0.05)$Metabolite


amp<-5

ggplot(subset(results.a,Description %in% comparisons & Metabolite %in% longevity),aes(x=Description,y=Metabolite))+
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

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_heatmap_longevity.pdf',sep = ''),
             width=6,height=6,useDingbats=FALSE)



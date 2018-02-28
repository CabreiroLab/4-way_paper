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


library(gtools)


devtools::install_github("PNorvaisas/PFun")
library(PFun)



theme_set(theme_light())

theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))

setwd("~/Dropbox/Projects/2015-Metformin/Metabolomics/")





odir<-'Summary_Ecoli_nucleotides_new'
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



#Needs to be manual

nuc.raw<-read.and.clean('Ecoli_nucleotides_New_extended.xlsx','Clean')
annot<-read.and.clean('Nucleotide_annotation.xlsx','Annotation')

nuc.m<-melt(nuc.raw,id.vars = c('Sample','Replicate'),variable.name = 'Metabolite',value.name = 'Conc')


nuc<-merge(nuc.m,annot,by.x = 'Metabolite',by.y='Compound',all.x = TRUE)
nuc<-subset(nuc,Metabolite!="")

nuc$Conc<-ifelse(nuc$Conc==0,NA,nuc$Conc)

nuc$Metformin<-ifelse(grepl("Metf",nuc$Sample),50,0)
nuc$Metformin<-as.factor(nuc$Metformin)
nuc$Strain<-gsub(" Metf","",nuc$Sample)
nuc$Glucose<-ifelse(grepl("Glu",nuc$Strain),20,0)
nuc$Glucose<-as.factor(nuc$Glucose)
nuc$Strain<-ifelse(grepl("Glu",nuc$Strain),"OP50-C",nuc$Strain)

strainlevels<-c('OP50-C','MR','CRP','cyaA')
strainlist<-c('OP50C','OP50MR','crp','cyaA')
nuc$Strain<-factor(nuc$Strain,levels=strainlevels,labels=strainlist)



nuc$Metabolite<-trimws(nuc$Metabolite)
nuc$Metabolite<-as.factor(nuc$Metabolite)
nuc$Conc<-as.numeric(nuc$Conc) 
nuc$logConc<-log(nuc$Conc,2)


nuc$Sample<-paste(nuc$Strain,nuc$Glucose,nuc$Metformin,sep='_')


nuc

unique(nuc$Sample)
# nuc$Replicate<-as.numeric(gsub("\\D", "", nuc$Sample))
# nuc$Type<-as.character(gsub("[[:digit:]]", "", nuc$Sample))
# nuc$Class<-as.factor(nuc$Class)
#nuc$Type<-factor(nuc$Type,levels=c('C','CM','R','RM'),labels=c('C','CM','R','RM'))




metslm<-dcast(nuc,Sample+Strain+Metformin+Glucose+Replicate~Metabolite,value.var = c('logConc'),fill = as.numeric(NA),drop=TRUE)
cols<-setdiff(colnames(metslm),c('Sample','Strain','Glucose','Metformin','Replicate'))

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

pca.dat<-metslm[cleanrows,cleancols]

ir.pca <- prcomp(pca.dat,
                 center = TRUE,
                 scale. = TRUE) 
plot(ir.pca,type='l')



pcadata<-data.frame(ir.pca$x)
pcadata[,c('Sample','Strain','Metformin','Glucose')]<-metslm[,c('Sample','Strain','Metformin','Glucose')]

pcadata<-merge(data.frame(ir.pca$x),metslm[,c('Sample','Strain','Metformin','Glucose')],by = 0)

pcaresult<-summary(ir.pca)$importance
PC1prc<-round(pcaresult['Proportion of Variance',][[1]]*100,0)
PC2prc<-round(pcaresult['Proportion of Variance',][[2]]*100,0)



ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  stat_ellipse(aes(group=interaction(Sample),linetype=Metformin,fill=factor( ifelse(Metformin==0,Strain,NA) ) ) ,
               level=0.68, type = "norm",geom="polygon",size=1,alpha=0.2)+
  geom_point(aes(shape=Glucose,fill=factor( ifelse(Metformin==0,Strain, NA ) ) ),size=5,stroke=2)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_shape_manual("Glucose, mM",values=c(21,23)) +
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=3,colour='black',linetype=c(1,3),fill=c(1,NA))))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA.pdf",sep=''),
             width=12,height=9)


allmets<-unique(as.character(nuc$Metabolite))
cleanmets<-setdiff(allmets,missing.cols)




#Heatmap
heatshape<-dcast(nuc,Metabolite~Sample+Replicate,value.var = 'logConc')

rownames(heatshape)<-heatshape$Metabolite
heatshape$Metabolite<-NULL

#heatshape<-heatshape[apply(heatshape,1,function(x) all(!is.na(x))) ,]

#Remove numbers from sample IDs
#groups<-gsub('_[[1-6]]?', '', colnames(heatshape))

groups<- substr(colnames(heatshape),1,nchar(colnames(heatshape))-2)
colnames(heatshape)
groups

annData<-factor(groups)
#ColSideColors<-rainbow(length(unique(annData)))[as.numeric(annData)]

ColSideColors<-rainbow(length(unique(annData)))[as.numeric(annData)]
legtxt<-as.character(unique(annData))
legcol<-rainbow(length(unique(annData)))[as.numeric(unique(annData))]

hmap<-heatmap3(as.matrix(heatshape),
               scale = 'row',
               #hclustfun = function(x) hclust(dist(x,method="manhattan"),method="complete"),
               #distfun = function(x) dist(x,method="manhattan"),
               #labRow=FALSE,
               ColSideLabs = 'Group',
               ColSideColors = ColSideColors,
               legendfun=function() showLegend(legend=legtxt,col=legcol,cex=1))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap.pdf",sep=''),
             width=9,height=9, useDingbats=FALSE)





ggplot(nuc,aes(x=Metabolite,y=logConc,color=Metformin))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  scale_y_continuous(breaks=seq(-20,20,by=1) )+
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  facet_wrap(Glucose~Strain,labeller = labeller(Glucose = label_both, Strain = label_both),ncol = 5) 

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_treatment.pdf",sep=''),
             width=25,height=10, useDingbats=FALSE)


ggplot(subset(nuc,Glucose==0),aes(x=Strain,y=logConc,color=Metformin))+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  geom_point(position = "dodge")+
  ggtitle('Comparison between strains - No Glucose. Boxplot: +/-SD, Min/Max')+
  scale_y_continuous(breaks=seq(-20,20,by=1) )+
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  facet_wrap(~Metabolite,nrow=3)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_metabolite_noGlu.pdf",sep=''),
             width=25,height=10, useDingbats=FALSE)


ggplot(subset(nuc,Strain=='OP50C'),aes(x=Glucose,y=logConc,color=Metformin))+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  geom_point(position = "dodge")+
  scale_y_continuous(breaks=seq(-20,20,by=1) )+
  ggtitle('Comparison in OP50C - Glucose. Boxplot: +/-SD, Min/Max')+
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  facet_wrap(~Metabolite,nrow=3)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_metabolite_Glu.pdf",sep=''),
             width=25,height=10, useDingbats=FALSE)







#Linear modelling


summary<-ddply(nuc,.(Metabolite,Sample),summarise,Mean=mean(logConc),sd=sd(logConc))
subset(summary,Metabolite=='cAMP')
#Select samples to work with

#Check order of contrasts
#sel.samples<-c("OP50C_0_0","OP50C_0_50")

sel.samples<-unique(nuc$Sample)

lmshape<-dcast(nuc,Sample+Replicate~Metabolite,value.var = 'logConc',drop=FALSE)

lmshape<-subset(lmshape,Sample %in% sel.samples)

unique(lmshape$Sample)

sample.names<-unique(as.character(lmshape$Sample))
metabolites<-setdiff(colnames(lmshape),c('Sample','Replicate'))




descriptors<-c('Strain')
contrasts<-read.contrasts('!Contrasts_Ecoli_nucleotide_metabolomics.xlsx','Contrasts_values',sel.samples)
contrasts.table<-contrasts$Contrasts.table
contr.matrix<-contrasts$Contrasts.matrix

strainlist<-c('OP50C','OP50MR','crp','cyaA')
#strainlist<-c('OP50C') #,'OP50MR','crp','cyaA'

contrasts.table$Strain<-factor(contrasts.table$Strain,levels=strainlist,labels=strainlist) #

contr.matrix


#Arrangement of sampels must correspond to the contrast table!!!
#lmshape$Sample<-factor(lmshape$Sample,levels=colnames(contr.matrix),labels=colnames(contr.matrix))




#lm2<-subset(lmshape,Sample %in% c("OP50C_0_0","OP50C_0_50") )

allresults<-hypothesise(lmshape,metabolites,contr.matrix)



results<-merge(allresults$All,contrasts.table[,c('Contrast','Description','Contrast_type','Strain')],by='Contrast',all.x=TRUE)
results<-rename(results,c("Variable"="Metabolite"))

results.castfull<-allresults$CastFull
results.cast<-allresults$Cast

subset(results,Metabolite=='cAMP')

results.exp<-results[,c('Contrast','Description','Contrast_type','Strain','Metabolite','logFC','SE','PE','NE','t.value','p.value','FDR')]
head(results.exp)
write.csv(results.exp,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)



write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)




#Plotting starts

erralpha<-1
errcolor<-'grey80'

brks<-seq(0,8,by=1)
gradcols<-c('black','purple','purple')


ggplot(subset(results,Contrast_type=="Treatment"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+
  ggtitle('Metformin treatment effect in different strains')+
  geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_treatment.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)


ggplot(subset(results,Contrast_type=="Strain"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+
  ylim(0,10)+
  ggtitle('Strain differeces in control')+
  geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_strain.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)



ggplot(subset(results,Contrast_type=="Interaction treatment-strain"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+

  geom_point()+
  #ylim(0,10)+
  ggtitle('Interaction between treatment and strain')+
  geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_interaction_treatment-strain.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)


ggplot(subset(results,Contrast_type=="Interaction treatment-glucose"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  scale_y_continuous(breaks=seq(-10,10,by=1) )+
  geom_point()+
  #ylim(0,10)+
  ggtitle('Interaction between treatment and glucose')+
  geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
  facet_wrap(~Description,ncol = 4)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_interaction_treatment-glucose.pdf',sep = ''),
             width=9,height=12,useDingbats=FALSE)


subset(results.exp,Metabolite=='cAMP' )




#Heatmap for summary

heatsum<-dcast(results,Metabolite~Description,value.var = 'logFC')
rownames(heatsum)<-heatsum$Metabolite
heatsum$Metabolite<-NULL


amp<-14

minv<- -amp
maxv<- amp

nstep<-maxv-minv


clrbrks<-seq(-amp,amp,by=2)

brks<-seq(minv,maxv,by=(maxv-minv)/(nstep))
bgg <- colorRampPalette(c("blue", "gray90", "red"))(n = nstep)

clrscale <- colorRampPalette(c("blue","blue", "gray90", "red","red"))(n = nstep)


reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }

hm<-heatmap3(as.matrix(heatsum),key=TRUE,Colv=FALSE,trace='none',col=bgg,
              xlab='Comparison',Rowv=TRUE,breaks = brks,dendrogram="row",scale="none")

ordmet<-rownames(heatsum[hm$rowInd,])


if (length(ordmet)!=length(unique(ordmet))){
  print("Non unique metabolites!")
}


#  "Treatment interaction for Glucose"   "Treatment interaction for crp"      
# [10] "Treatment interaction for cyaA"                          
# [13] "Treatment interaction for OP50MR

comparisons<-c("Treatment effect on OP50C","Treatment effect on OP50C + Glucose","Treatment effect on OP50MR",
"Treatment effect on crp","Treatment effect on cyaA","Glucose effect on OP50C","Strain effect in OP50MR",
"Strain effect in crp","Strain effect in cyaA")


results$Metabolite<-factor(results$Metabolite,levels=ordmet,labels=ordmet)
results$Description<-factor(results$Description,levels=comparisons,labels=comparisons)
results$FDRstars<-stars.pval(results$FDR)



ggplot(subset(results,Description %in% comparisons),aes(x=Description,y=Metabolite))+
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
             width=9,height=9,useDingbats=FALSE)






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

#devtools::install_github("PNorvaisas/PFun")
library(PFun)




theme_set(theme_light())

theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))

setwd("~/Dropbox/Projects/2015-Metformin/Metabolomics/")



odir<-'Summary_Ecoli_100met'
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


mycolsplit<-function(column,cols,sep=':'){
  elems<-unlist(strsplit(column,paste("([\\",sep,"])",sep=''), perl=TRUE))
  return(as.data.frame(matrix(elems,ncol=cols,byrow=TRUE)))
}

MinMeanSDMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}




#Needs to be manual

met.raw<-read.and.clean('Ecoli_100met_data.xlsx','Export')
norm<-read.and.clean('Ecoli_100met_normalisation.xlsx','Normalisation')

dim(met.raw)

miss.cols<-apply(met.raw, 2, function(x) any(x==""))
miss.rows<-apply(met.raw, 1, function(x) any(x==""))

missing.cols<-names(miss.cols[miss.cols==TRUE])
missing.rows<-rownames(met.raw)[miss.rows==TRUE]

missing.cols
missing.rows

met<-met.raw[,-match(missing.cols,colnames(met.raw))]
dim(met)


metm<-melt(met,id.vars = c('Sample','Replicate'),variable.name = 'Metabolite',value.name = 'Conc')

metm[,c('Number','Sample')]<-mycolsplit(metm$Sample,2,sep='_')
#metm[,c('Sample','Color')]<-mycolsplit(metm$Sample,2,sep='\w')

metm$Sample



metn<-merge(metm,norm[,c('Sample','Code','Protein (ug/ml)')],by.x='Sample',by.y='Sample',all.x=TRUE)
#metn$Replicate<-as.numeric(gsub("[^\\d]+", "", metn$Sample, perl=TRUE))

metn$Sample<-metn$Code
metn$Code<-paste(metn$Sample,metn$Replicate,sep='_')
head(metn)

#metn$Sample<-gsub('[[:digit:]]+','',metn$Sample)
#metn$Sample<-ifelse(metn$Sample %in% c('CM','RM'),metn$Sample,paste(metn$Sample,'0',sep=''))
metn$Metabolite<-trimws(metn$Metabolite)
metn$Metabolite<-as.factor(metn$Metabolite)


samples<-c('C0','CM','R0','RM','P0')
strains<-c('OP50-C','OP50-MR','OP50-P5')

metn$Sample<-factor(metn$Sample,levels=samples,labels=samples)

metn$Strain<-ifelse(metn$Sample %in% c('C0','CM'),'OP50-C',NA)
metn$Strain<-ifelse(metn$Sample %in% c('R0','RM'),'OP50-MR',metn$Strain)
metn$Strain<-ifelse(metn$Sample %in% c('P0','PM'),'OP50-P5',metn$Strain)
metn$Strain<-factor(metn$Strain,levels=strains,labels=strains)
metn$Metformin<-as.factor(ifelse(metn$Sample %in% c('CM','RM'),50,0))


metn$Conc<-as.numeric(metn$Conc)
metn$`Protein (ug/ml)`<-as.numeric(metn$`Protein (ug/ml)`)
metn$Concnorm<-metn$Conc / metn$`Protein (ug/ml)`
metn$logConc<-log(metn$Concnorm)



unique(metn$Sample)
unique(metn$Strain)


subset(metn,Code %in% c('RM_2','R0_3'))



metslm<-dcast(metn,Code+Sample+Strain+Metformin+Replicate~Metabolite,value.var = c('logConc'),fill = as.numeric(NA),drop=TRUE)
rownames(metslm)<-metslm$Code
cols<-setdiff(colnames(metslm),c('Code','Sample','Strain','Metformin','Replicate'))

head(metslm)


#Find compounds with missing values
miss.cols<-apply(metslm, 2, function(x) any(is.na(x)))
miss.rows<-apply(metslm, 1, function(x) any(is.na(x)))

missing.cols<-names(miss.cols[miss.cols==TRUE])
missing.rows<-rownames(metslm)[miss.rows==TRUE]

missing.cols
missing.rows


cleancols<-setdiff(cols,missing.cols)
cleanrows<-setdiff(rownames(metslm),c('RM_2','R0_3'))


pca.dat<-metslm[cleanrows,cleancols]

ir.pca <- prcomp(pca.dat,
                 center = TRUE,
                 scale. = TRUE) 
plot(ir.pca,type='l')



pcadata<-data.frame(ir.pca$x)
#pcadata[,c('Sample','Strain','Metformin')]<-metslm[,c('Sample','Strain','Metformin')]

pcadata<-merge(data.frame(ir.pca$x),metslm[,c('Code','Sample','Strain','Metformin','Replicate')],by = 0)

pcaresult<-summary(ir.pca)$importance
PC1prc<-round(pcaresult['Proportion of Variance',][[1]]*100,0)
PC2prc<-round(pcaresult['Proportion of Variance',][[2]]*100,0)



ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  stat_ellipse(aes(group=interaction(Sample),linetype=Metformin,fill=factor( ifelse(Metformin==0,Strain,NA) ) ) ,
               level=0.68, type = "norm",geom="polygon",size=1,alpha=0.2)+
  geom_point(aes(fill=factor( ifelse(Metformin==0,Strain, NA ) ) ),size=5,stroke=2,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  geom_text(aes(label=Code),nudge_y=1)+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=3,colour='black',linetype=c(1,3),fill=c(1,NA))))+
  guides(colour = guide_legend(override.aes = list(shape=c(1),size=3)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

  

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA.pdf",sep=''),
             width=12,height=9)


allmets<-unique(as.character(metn$Metabolite))
#cleanmets<-setdiff(allmets,missing.cols)




#Heatmap
heatshape<-dcast(subset(metn,!Code %in% c('RM_2','R0_3')),Metabolite~Code,value.var = 'logConc')

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





ggplot(metn,aes(x=Strain,y=logConc,color=Metformin))+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  geom_point(position = "dodge")+
  ggtitle('Comparison between strains. Boxplot: +/-SD, Min/Max')+
  scale_y_continuous(breaks=seq(-20,20,by=1) )+
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  facet_wrap(~Metabolite,ncol=6)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_metabolite.pdf",sep=''),
             width=10,height=30, useDingbats=FALSE)




#Linear modelling

lmshape<-dcast(subset(metn,!Code %in% c("RM_2",'R0_3')),Sample+Replicate~Metabolite,value.var = 'logConc',drop=TRUE)


sample.names<-unique(as.character(lmshape$Sample))
metabolites<-setdiff(colnames(lmshape),c('Sample','Replicate'))


contrasts<-read.contrasts('!Contrasts_Ecoli_100met_metabolomics.xlsx','Contrasts_values',sample.names)
contrasts.table<-contrasts$Contrasts.table
cont.matrix<-contrasts$Contrasts.matrix

contrasts.table

cont.matrix

contrasts.table$Strain<-factor(contrasts.table$Strain,levels=strains,labels=strains)


#Arrangement of sampels must correspond to the contrast table!!!

allresults<-hypothesise(lmshape,metabolites,cont.matrix)

results<-merge(allresults$All,contrasts.table[,c('Contrast','Description','Contrast_type','Strain')],by='Contrast',all.x=TRUE)
results<-rename(results,c("Variable"="Metabolite"))

results.castfull<-allresults$CastFull
results.cast<-allresults$Cast


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


ggplot(results,aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+
  ggtitle('Results from different comparisons')+
  geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_all.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)
# 
# 
# ggplot(subset(results,Contrast_type=="Strain"),aes(x=logFC,y=logFDR,color=Strain))+
#   geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
#   geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
#   geom_point()+
#   ylim(0,10)+
#   ggtitle('Strain differences in control')+
#   geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
#   facet_wrap(~Description)
# 
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_strain.pdf',sep = ''),
#              width=12,height=9,useDingbats=FALSE)
# 
# 
# 
# ggplot(subset(results,Contrast_type=="Interaction treatment-strain"),aes(x=logFC,y=logFDR,color=Strain))+
#   geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
#   geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
#   geom_point()+
#   #ylim(0,10)+
#   ggtitle('Interaction between treatment and strain')+
#   geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
#   facet_wrap(~Description)
# 
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_interaction_treatment-strain.pdf',sep = ''),
#              width=12,height=9,useDingbats=FALSE)
# 
# 
# ggplot(subset(results,Contrast_type=="Interaction treatment-glucose"),aes(x=logFC,y=logFDR,color=Strain))+
#   geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
#   geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
#   geom_point()+
#   #ylim(0,10)+
#   ggtitle('Interaction between treatment and glucose')+
#   geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
#   facet_wrap(~Description)
# 
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_interaction_treatment-glucose.pdf',sep = ''),
#              width=9,height=9,useDingbats=FALSE)



#Heatmap for summary

heatsum<-dcast(results,Metabolite~Description,value.var = 'logFC')
rownames(heatsum)<-heatsum$Metabolite
heatsum$Metabolite<-NULL



amp<-6

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

comparisons<-c("Treatment effect on OP50-C",
               "Treatment effect on OP50-MR",
               "Stain effect in OP50-MR",
               "Strain effect in OP50-P5",
               "Treatment interaction for OP50MR")


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
             width=6,height=20,useDingbats=FALSE)



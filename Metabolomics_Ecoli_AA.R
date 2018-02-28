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





odir<-'Summary_Ecoli_AA'
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

met.raw<-read.xlsx2('Ecoli_old_AA.xlsx',sheetName='AA',
                    header=TRUE,endRow=17,
                    colClasses=c("character", rep("numeric", 24)))
met.raw<-subset(met.raw,Code!='CM1')
AAs<-setdiff(colnames(met.raw),'Code')


rowmeans<-apply(met.raw[,AAs],1,mean)
allmean<-mean(rowmeans,na.rm = TRUE)
norms<-rowmeans/allmean

met.norm<-met.raw
met.norm[,AAs]<-met.raw[,AAs]/norms



met.m<-melt(met.norm,id.vars = c('Code'),variable.name = 'AA',value.name = 'Conc')


#write.csv(metabolites,'AA.csv')

annot<-read.and.clean('AA_annotation.xlsx','Annotation')
met<-merge(met.m,annot,by.x = 'AA',by.y='AA',all.x = TRUE)


met<-subset(met,Metabolite!="")

met$Conc<-ifelse(met$Conc==0,NA,met$Conc)

met$Sample<-gsub('[[:digit:]]', '', met$Code)
met$Sample<-ifelse(met$Sample %in% c('C','R'),paste(met$Sample,'0',sep=''),met$Sample)
met$Replicate<-gsub('[[:alpha:]]', '', met$Code)
met$Metformin_mM<-ifelse(grepl("M",met$Sample),50,0)
met$Metformin_mM<-as.factor(met$Metformin_mM)
met$Strain<-ifelse(grepl("R",met$Sample),'OP50-MR','OP50-C')


strains<-c('OP50-C','OP50-MR')
met$Strain<-factor(met$Strain,levels=strains,labels=strains)
met


met$Metabolite<-trimws(met$Metabolite)
met$Metabolite<-as.factor(met$Metabolite)
met$Conc<-as.numeric(met$Conc) 
met$logConc<-log(met$Conc,2)

unique(met$Code)
metabolites

metabolites<-as.character(unique(met$Metabolite))

# nuc$Replicate<-as.numeric(gsub("\\D", "", nuc$Sample))
# nuc$Type<-as.character(gsub("[[:digit:]]", "", nuc$Sample))
# nuc$Class<-as.factor(nuc$Class)
#nuc$Type<-factor(nuc$Type,levels=c('C','CM','R','RM'),labels=c('C','CM','R','RM'))

metslm<-dcast(met,Code+Sample+Strain+Metformin_mM+Replicate~Metabolite,value.var = c('logConc'),fill = as.numeric(NA),drop=TRUE)
rownames(metslm)<-metslm$Code
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
pcadata[,c('Sample','Strain','Metformin_mM')]<-metslm[,c('Sample','Strain','Metformin_mM')]

#pcadata<-merge(data.frame(ir.pca$x),metslm[,c('Sample','Strain','Metformin_mM')],by = 0)

pcaresult<-summary(ir.pca)$importance
PC1prc<-round(pcaresult['Proportion of Variance',][[1]]*100,0)
PC2prc<-round(pcaresult['Proportion of Variance',][[2]]*100,0)



ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  stat_ellipse(aes(group=interaction(Sample),linetype=Metformin_mM,fill=factor( ifelse(Metformin_mM==0,Strain,NA) ) ) ,
               level=0.68, type = "norm",geom="polygon",size=1,alpha=0.2)+
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=5,stroke=2,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(size=3,colour='black',linetype=c(1,3),fill=c(1,NA),shape=c(21,21))))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA.pdf",sep=''),
             width=12,height=9)


#cleanmets<-setdiff(allmets,missing.cols)




#Heatmap
heatshape<-dcast(met,Metabolite~Sample+Replicate,value.var = 'logConc')

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




ggplot(met,aes(x=Metabolite,y=logConc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  scale_y_continuous(breaks=seq(-3,8,by=1) )+
  ylab('log 2 Concentration, nmol')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Strain,labeller = labeller(Strain = label_both),ncol = 5)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_treatment.pdf",sep=''),
             width=25,height=10, useDingbats=FALSE)

ggplot(met,aes(x=Metabolite,y=Conc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  ylab('Concentration, nmol')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Strain,labeller = labeller(Strain = label_both),ncol = 5)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_Conc_by_treatment.pdf",sep=''),
             width=25,height=10, useDingbats=FALSE)



#Linear modelling

met.mes<-melt(met,measure.vars = c('Conc','logConc'),variable.name = 'Measure',value.name = 'Value')
summary<-ddply(met.mes,.(Metabolite,Strain,Metformin_mM,Sample,Measure),summarise,Mean=mean(Value),SD=sd(Value))

met.abs<-ddply(met.mes,.(Metabolite,Measure),summarise,Mean=mean(Value),SD=sd(Value))

sum.m<-melt(summary,measure.vars = c('Mean','SD'),variable.name = 'Stat',value.name = 'Value')

met.am<-melt(met.abs,measure.vars = c('Mean','SD'),variable.name = 'Stat',value.name = 'Value')
sum.c<-dcast(sum.m,Metabolite+Strain+Metformin_mM+Sample~Measure+Stat,value.var = 'Value',drop = TRUE)

met.c<-dcast(met.am,Metabolite~Measure+Stat,value.var = 'Value',drop = TRUE)


subset(summary,Metabolite=='L-Glutamate')
#Select samples to work with

#Check order of contrasts
#sel.samples<-c("OP50C_0_0","OP50C_0_50")

sel.samples<-unique(met$Sample)

lmshape<-dcast(met,Sample+Replicate~Metabolite,value.var = 'logConc',drop=FALSE)

lmshape<-subset(lmshape,Sample %in% sel.samples)

unique(lmshape$Sample)

sample.names<-unique(as.character(lmshape$Sample))
metabolites<-setdiff(colnames(lmshape),c('Sample','Replicate'))

contrasts<-read.contrasts('!Contrasts_Ecoli_AA_metabolomics.xlsx','Contrasts_values',sel.samples)
contrasts.table<-contrasts$Contrasts.table
contr.matrix<-contrasts$Contrasts.matrix

strainlist<-c('OP50-C','OP50-MR')
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
results.castfull<-rename(results.castfull,c("Variable"="Metabolite"))
results.cast<-allresults$Cast
results.cast<-rename(results.cast,c("Variable"="Metabolite"))


results.joined<-merge(results,met.c,by='Metabolite')



subset(results,Metabolite=='L-Glutamate')

results.exp<-results[,c('Contrast','Description','Contrast_type','Strain','Metabolite','logFC','SE','PE','NE','t.value','p.value','FDR')]
head(results.exp)


write.csv(results.exp,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)




#Plotting starts

erralpha<-0.6
errcolor<-'grey80'
lblsize<-2
amp<- 3
cbrks<-seq(0,amp,by=1)
gradcols<-c('black','purple','purple')
intname<-'Interaction strength'




results.castfull.met<-merge(results.castfull,met.c,by='Metabolite')


ggplot(results.castfull.met,aes(x=`CM-C_logFC`,y=`RM-R_logFC`,color=abs(`RM-R-(CM-C)_logFC`)))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=`CM-C_NE`,xmax=`CM-C_PE`),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=`RM-R_NE`,ymax=`RM-R_PE`),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=abs(`RM-R-(CM-C)_logFC`)))+
  ggtitle('Interaction between treatment and strain')+
  #labs(color='Interaction strength')+
  geom_text_repel(aes(label=Metabolite),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(0,amp),name=intname)+
  scale_size(range = c(0.25, 5),name=intname)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_treatment.pdf',sep = ''),
             width=9,height=5,useDingbats=FALSE)



aamp<- 200
acbrks<-seq(0,aamp,by=50)
abscols<-c('black','blue','purple')

aintname<-'Average metabolite\nconc., nmol'


ggplot(results.castfull.met,aes(x=`CM-C_logFC`,y=`RM-R_logFC`,color=abs(`RM-R-(CM-C)_logFC`)))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=`CM-C_NE`,xmax=`CM-C_PE`),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=`RM-R_NE`,ymax=`RM-R_PE`),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  ggtitle('Interaction between treatment and strain')+
  #labs(color='Interaction strength')+
  geom_text_repel(aes(label=Metabolite),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(0,amp),name=intname)+
  scale_size(range = c(0.25, 5),name=aintname)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_treatment_Conc.pdf',sep = ''),
             width=9,height=5,useDingbats=FALSE)





ggplot(subset(results.joined,Contrast_type=="Treatment"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  labs(size='Average metabolite\nconcentration, nmol')+
  ggtitle('Metformin treatment effect in different strains')+
  geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.3,size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_treatment.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)




ggplot(subset(results.joined,Contrast_type=="Strain"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  labs(size='Average metabolite\nconcentration, nmol')+
  ylim(0,10)+
  ggtitle('Strain differences in control')+
  geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_strain.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)



ggplot(subset(results.joined,Contrast_type=="Interaction treatment-strain"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  labs(size='Average metabolite\nconcentration, nmol')+
  #ylim(0,10)+
  ggtitle('Interaction between treatment and strain')+
  geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_interaction_treatment-strain.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)


subset(results.exp,Metabolite=='L-Glutamate' )




ggplot(subset(results.joined,Contrast_type=="Treatment"),aes(x=Conc_Mean,y=logFC,color=Strain))+
  #geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=Conc_Mean-Conc_SD,xmax=Conc_Mean+Conc_SD),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=NE,ymax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)) )+
  scale_size_continuous(breaks=seq(0,200,by=50))+
  xlab('Average metabolite\nconcentration, nmol')+
  labs(size='Average metabolite\nconcentration, nmol')+
  ggtitle('Correlation between absolute concentration and fold changes')+
  geom_text_repel(aes(label=Metabolite),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  annotation_logticks(sides="b")+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Correlation_between_conc_and_logFC.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)




#Heatmap for summary

heatsum<-dcast(results,Metabolite~Description,value.var = 'logFC')
rownames(heatsum)<-heatsum$Metabolite
heatsum$Metabolite<-NULL


amp<-4

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

comparisons<-c("Treatment effect on OP50-C","Treatment effect on OP50-MR",
              "Strain effect in OP50-MR","Treatment interaction for OP50MR")


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
             width=6,height=9,useDingbats=FALSE)






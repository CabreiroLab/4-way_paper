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

library(ellipse)

library(Unicode)

library(ggthemes)


devtools::install_github("PNorvaisas/PFun")
library(PFun)



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

theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))

setwd("~/Dropbox/Projects/2015-Metformin/Metabolomics/")





odir<-'Summary_Celegans_AA_old'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



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



#Needs to be manual

met.raw<-read.xlsx2('Celegans_old_AA.xls',sheetName='AA_worms_clean',
                    header=TRUE,endRow=17,
                    colClasses=c(rep("character",6), rep("numeric", 19)))








head(met.raw)

AAs<-colnames(met.raw)[7:length(colnames(met.raw))]


rowmeans<-apply(met.raw[,AAs],1,mean)
allmean<-mean(rowmeans,na.rm = TRUE)
norms<-rowmeans/allmean

met.norm<-met.raw

#No normalisation?
#met.norm[,AAs]<-met.raw[,AAs]/norms



met<-melt(met.norm,measure.vars = AAs,variable.name = 'Metabolite',value.name = 'Conc')

met$Metformin_mM<-as.factor(met$Metformin_mM)
met$Group<-as.factor(met$Group)



strains<-c('OP50','OP50-MR')
met$Strain<-factor(met$Strain,levels=strains,labels=strains)
met


met$Metabolite<-trimws(met$Metabolite)
met$Metabolite<-as.factor(met$Metabolite)

met$Conc<-as.numeric(met$Conc) 
met$logConc<-log(met$Conc,2)


#met<-subset(met,Metabolite!='Citrulline')


# 
# met<-subset(met,!( Group %in% c('C_C','C_T','crp_C','crp_T') & Replicate %in% c(1,4,6) |
#                     Group %in% c('MR_C','MR_T') & Replicate %in% c(1,2) )    )
# 


unique(met$Group)
metabolites<-as.character(unique(met$Metabolite))
metabolites

# 
# removals<-read.xlsx2('Metformin AA worm data.xlsx',sheetName='Removals',
#                      header=TRUE,
#                      colClasses=c(rep("character",2)))
# 
# 
# removals


met.clean<-met
# for (i in 1:nrow(removals)){
#   met.clean<-subset(met.clean,!(Sample==as.character(removals[i,'Sample']) & Metabolite==as.character(removals[i,'Metabolite'])) )
# }

metslm<-dcast(met,Sample+Group+Strain+Metformin_mM+Replicate~Metabolite,value.var = c('logConc'),drop=TRUE)#,fill = as.numeric(NA)
rownames(metslm)<-metslm$Sample
head(metslm)

#View(metslm)


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



ellipses<-ddply(pcadata,.(Strain,Group,Metformin_mM), summarise, x=getellipse(PC1,PC2,1)$x,y=getellipse(PC1,PC2,0.75)$y ) 


ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=5,stroke=1,shape=21)+
  scale_linetype_manual("Drug, mM",values=c("0"=1,"50"=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA.pdf",sep=''),
             width=12,height=9)


#cleanmets<-setdiff(allmets,missing.cols)





#Heatmap
heatshape<-dcast(met.clean,Metabolite~Group+Replicate,value.var = 'logConc')

rownames(heatshape)<-heatshape$Metabolite
heatshape$Metabolite<-NULL

#heatshape<-heatshape[apply(heatshape,1,function(x) all(!is.na(x))) ,]

#Remove numbers from sample IDs
#groups<-gsub('_[[1-6]]?', '', colnames(heatshape))





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
               #labRow=FALSE,
               #key=TRUE,
               #ColSideLabs = 'Group',
               
               ColSideColors = clab)
legend('topleft',legend=legtxt,fill=legcol, border=FALSE, bty="n",title='Strain/Mutant')
legend('left',legend=c('Treatment','Control'),fill=c('black','white'), border=TRUE, bty="n",title='Drug')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap.pdf",sep=''),
             width=12,height=9, useDingbats=FALSE)






ggplot(met.clean,aes(x=Metabolite,y=logConc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point()+
  geom_text(aes(label=Replicate),color='black',size=2)+
  scale_y_continuous(breaks=seq(-3,8,by=1) )+
  ylab('log 2 Concentration, nmol')+

  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Strain,labeller = labeller(Strain = label_both),ncol = 5)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_treatment.pdf",sep=''),
             width=25,height=10, useDingbats=FALSE)

ggplot(met.clean,aes(x=Metabolite,y=Conc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point()+
  geom_text(aes(label=Replicate),color='black',size=2)+
  ylab('Concentration, nmol')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Strain,labeller = labeller(Strain = label_both),ncol = 5)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_Conc_by_treatment.pdf",sep=''),
             width=25,height=10, useDingbats=FALSE)




ggplot(met.clean,aes(x=Strain,y=logConc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point()+
  scale_y_continuous(breaks=seq(-5,15,by=1))+
  geom_text(aes(label=Replicate),color='black',size=2)+
  ylab('log 2 Concentration, nmol')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Metabolite,ncol = 5)

#,labeller = labeller(Metablite = label_both)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_metabolite.pdf",sep=''),
             width=20,height=20, useDingbats=FALSE)







head(met)


unique(as.character(met$Group))

#Summarise manually

met.mes<-melt(met.clean,measure.vars = c('Conc','logConc'),variable.name = 'Measure',value.name = 'Value')
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
  scale_y_continuous(breaks = seq(0,100,by=5))+
  ylab('In-group variation, %')+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=5,height=20, useDingbats=FALSE)



#Find outliers in most variable groups
# indxord<-sum.c[order(sum.c$Conc_SD),'Index']
# sum.c$Index<-factor(sum.c$Index,levels=indxord,labels=indxord)
# 
# ggplot(sum.c,aes(x=Index,y=Conc_SD))+
#   geom_point()+
#   coord_flip()
# 
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Ingroup_variation_Conc_SD.pdf",sep=''),
#              width=5,height=20, useDingbats=FALSE)










#Absolute concentrations over all samples
met.abs<-ddply(met.mes,.(Metabolite,Measure),summarise,Mean=mean(Value),SD=sd(Value))
met.am<-melt(met.abs,measure.vars = c('Mean','SD'),variable.name = 'Stat',value.name = 'Value')
#Average concentration over all samples
met.c<-dcast(met.am,Metabolite~Measure+Stat,value.var = 'Value',drop = TRUE)

head(met.am)

head(met.c)

sel.groups<-unique(met$Group)
sel.groups








#Linear modelling

lmshape<-dcast(met.clean,Group+Replicate~Metabolite,value.var = 'logConc',drop=TRUE)
#lmshape$Sample<-lmshape$Group
head(lmshape)

lmshape<-subset(lmshape,Group %in% sel.groups)

unique(lmshape$Group)


metabolites<-setdiff(colnames(lmshape),c('Sample','Group','Replicate'))




contrasts<-read.contrasts('!Contrasts_Celegans_AA_old_metabolomics.xlsx','Contrasts_values',sel.groups)
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


results.jT<-merge(results.cm,subset(results,Contrast_type=='Treatment' & Contrast!='C_T-C_C'),by='Metabolite',suffixes = c('_C',''),all.y=TRUE)




head(results.jT)
#Plotting starts

erralpha<-0.6
errcolor<-'grey80'
lblsize<-2
amp<- 2
cbrks<-seq(-amp,amp,by=1)
#gradcols<-c('black','purple','purple')
maincomp<-'Interaction strength'


gradcols<-c('blue4','blue','gray80','red','red4')

ggplot(results.jT,aes(x=logFC_C,y=logFC,color=logFC-logFC_C))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=NE_C,xmax=PE_C),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=NE,ymax=PE),alpha=erralpha,color=errcolor,width=0)+
  geom_point()+
  geom_text_repel(aes(label=Metabolite),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  xlab('Metformin effect on OP50 (N2)')+
  ylab('Metformin effect on other strain')+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  facet_wrap(~Strain)+
  theme(panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_treatment.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)



ggplot(subset(results.a,Contrast_type=="Treatment"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  ylim(0,15)+
  labs(size='Average metabolite\nconcentration, nmol')+
  ggtitle('Metformin treatment effect in different strains')+
  geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_treatment.pdf',sep = ''),
             width=12,height=12,useDingbats=FALSE)


ggplot(subset(results.a,Contrast_type %in% c("Bacterial mutant","Worm mutant","Bacterial strain") ),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  ylim(0,15)+
  labs(size='Average metabolite\nconcentration, nmol')+
  ggtitle('Mutant difference vs OP50 (N2)')+
  #geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_mutant.pdf',sep = ''),
             width=12,height=12,useDingbats=FALSE)




ggplot(subset(results.a,Contrast_type=="Interaction"),aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=Conc_Mean))+
  labs(size='Average metabolite\nconcentration, nmol')+
  #ylim(0,10)+
  ggtitle('Interaction between treatment and strain in comparison to OP50 (N2)')+
  #geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),size=2)+
  facet_wrap(~Description)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_interaction_treatment-strain.pdf',sep = ''),
             width=12,height=12,useDingbats=FALSE)



# ggplot(subset(results.jT,Contrast_type=="Treatment"),aes(x=Conc_Mean,y=logFC,color=Strain))+
#   #geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
#   geom_errorbarh(aes(xmin=Conc_Mean-Conc_SD,xmax=Conc_Mean+Conc_SD),alpha=erralpha,color=errcolor,height=0)+
#   geom_errorbar(aes(ymin=NE,ymax=PE),alpha=erralpha,color=errcolor,height=0)+
#   geom_point(aes(size=Conc_Mean))+
#   scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x)) )+
#   scale_size_continuous(breaks=seq(0,200,by=50))+
#   xlab('Average metabolite\nconcentration, nmol')+
#   labs(size='Average metabolite\nconcentration, nmol')+
#   ggtitle('Correlation between absolute concentration and fold changes')+
#   geom_text_repel(aes(label=Metabolite),
#                   size=lblsize,
#                   force=2,
#                   segment.colour=errcolor,
#                   segment.alpha =erralpha)+
#   annotation_logticks(sides="b")+
#   facet_wrap(~Description)
# 
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Correlation_between_conc_and_logFC.pdf',sep = ''),
#              width=12,height=6,useDingbats=FALSE)
# 



#Heatmap for summary


unique(as.character(results.a$Description))

comparisons<-c("Treatment effect on OP50","Treatment effect on OP50-MR",
               "Mutant difference for OP50-MR","Interaction with treatment for OP50-MR")



heatsum<-dcast(results.a,Metabolite~Description,value.var = 'logFC')
rownames(heatsum)<-heatsum$Metabolite
heatsum$Metabolite<-NULL


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


results.a$Metabolite<-factor(results.a$Metabolite,levels=ordmet,labels=ordmet)
results.a$Description<-factor(results.a$Description,levels=comparisons,labels=comparisons)
results.a$FDRstars<-stars.pval(results.a$FDR)



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
             width=6,height=9,useDingbats=FALSE)






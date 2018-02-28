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






setwd("~/Dropbox/Projects/2015-Metformin/Metabolomics/")





odir<-'Summary_Metabolomics_comparisons'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



ecoli.AA<-read.csv('Summary_Ecoli_AA/All_results.csv',header = TRUE)
ecoli.AA$Dataset<-'AA'
ecoli.nuc<-read.csv('Summary_Ecoli_nucleotides/All_results.csv',header = TRUE)
ecoli.nuc$Dataset<-'Nuc'
ecoli.100<-read.csv('Summary_Ecoli_100met/All_results.csv',header = TRUE)
ecoli.100$Dataset<-'100'


head(ecoli.100)


ecoli.HMT<-read.csv('HMT_MPFT/Results.csv',header = TRUE,check.names = FALSE)#,colClasses = c(rep('character',5),rep('numeric',2))
ecoli.HMT$Dataset<-'HMT'


head(ecoli.HMT)

ecoli.HMT<-subset(ecoli.HMT,grepl('[[:digit:]]',`OP50-C`) & grepl('[[:digit:]]',`OP50-D`) & `HMDB ID`!='No ID')


ecoli.HMT$`OP50-C`<-as.numeric(as.character(ecoli.HMT$`OP50-C`))
ecoli.HMT$`OP50-D`<-as.numeric(as.character(ecoli.HMT$`OP50-D`))


ecoli.HMT$Metf_logFC<-log2(ecoli.HMT$`OP50-D`/ecoli.HMT$`OP50-C`)



mets.HMT<-unique(as.character(ecoli.HMT$`Compound name`))
mets.100<-unique(as.character(ecoli.100$Metabolite))

write.csv(mets.100,paste(odir,'/Ecoli_100_metabolites.csv',sep = ''))


intersect(mets.HMT,mets.100)



HMDB<-read.csv('Summary_Ecoli_100met/HMDB_matches.csv',header=TRUE,check.names = TRUE)
ecoli.100.HMDB<-merge(subset(ecoli.100,Contrast=='C_T-C_C'),HMDB[,c('Query','HMDB')],by.x='Metabolite',by.y='Query',all.x=TRUE)



ecoli.100.HMT<-merge(ecoli.100.HMDB,ecoli.HMT[,c('HMDB ID','Metf_logFC')],by.x='HMDB',by.y='HMDB ID',all=TRUE)



View(ecoli.100.HMT)


erralpha<-0.6
errcolor<-'grey80'
lblsize<-2



ggplot(subset(ecoli.100.HMT,!is.na(logFC)),aes(x=logFC,y=Metf_logFC,color=FDR<0.05))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  xlab('Ecoli 100 metabolite logFC')+
  ylab('Ecoli HMT logFC')+
  labs(color='E. coli 100\nFDR<0.05')+
  geom_point()+
  geom_text_repel(aes(label=Metabolite))


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_Ecoli_100met_HMT_comparison.pdf',sep = ''),
             width=10,height=6,useDingbats=FALSE)




cel.AA<-read.csv('Summary_Celegans_AA/All_results.csv',header = TRUE)
cel.AA$Dataset<-'AA'
cel.AA$Common.name<-cel.AA$Metabolite
cel.FA<-read.csv('Summary_Celegans_FA/All_results.csv',header = TRUE)
cel.FA$Dataset<-'FA'


head(ecoli.AA)
head(ecoli.nuc)
head(ecoli.100)

ecoli.comb<-rbind(ecoli.100,ecoli.AA,ecoli.nuc)
ecoli.comb$Organism<-'Eco'
ecoli.comb$Common.name<-ecoli.comb$Metabolite

head(cel.FA)
head(cel.AA)

cel.comb<-rbind(cel.FA,cel.AA)
cel.comb$Organism<-'Cel'

all<-rbind(cel.comb,ecoli.comb)



head(all)


stats<-c('logFC','FDR','p.value','SE','PE','NE','t.value')

all.m<-melt(all,measure.vars = stats,variable.name = 'Stat',value.name = 'Value')


head(all)
write.csv(all,paste(odir,'/All_results.csv',sep = ''))




head(all)

subset(all,Metabolite=='cAMP')

maincontrast<-c("C_T-C_C","MR_T-MR_C","MR_C-C_C","MR_T-MR_C-(C_T-C_C)")


unique(as.character(all[duplicated(all[,c('Contrast','Organism','Metabolite')]),'Metabolite']))



all.m$Organism_indx<-makereplicates(all.m[,c('Contrast','Organism','Metabolite','Stat')])






all.c<-dcast(all.m,Strain+Contrast+Metabolite+Common.name~Organism+Dataset+Stat,value.var = 'Value')
all.co<-dcast(all.m,Strain+Contrast+Metabolite+Common.name+Organism_indx~Organism+Stat,value.var = 'Value')


sharedecocel<-as.character(subset(all.co,!is.na(Cel_logFC) & !is.na(Eco_logFC) )[,'Metabolite'])



eco.shared<-subset(all.m,Metabolite %in% sharedecocel & Contrast %in% maincontrast & Organism=='Eco')
cel.shared<-subset(all.m,Metabolite %in% sharedecocel & Contrast %in% maincontrast & Organism=='Cel')


cel.shared.c<-dcast(cel.shared,Contrast+Metabolite~Stat,value.var = 'Value')
eco.shared.c<-dcast(eco.shared,Contrast+Metabolite+Dataset~Stat,value.var = 'Value')
  
ec.shared<-merge(eco.shared.c,cel.shared.c,by=c('Contrast','Metabolite'),all=TRUE,suffixes = c('_eco','_cel'))




erralpha<-0.6
errcolor<-'grey80'
lblsize<-2

#gradcols<-c('black','purple','purple')
maincomp<-'Difference'


gradcols<-c('blue4','blue','gray80','red','red4')

amp<-4
cbrks<-seq(-amp,amp,by=1)
camp<-5
ggplot(subset(all.c,Contrast %in% maincontrast ),
       aes(x=Eco_100_logFC,y=Eco_Nuc_logFC,color=Eco_Nuc_logFC-Eco_100_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=Eco_100_NE,xmax=Eco_100_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=Eco_Nuc_NE,ymax=Eco_Nuc_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_point()+
  scale_x_continuous(breaks=seq(-10,10,by=1),limits = c(-camp,camp))+
  scale_y_continuous(breaks=seq(-10,10,by=1),limits = c(-camp,camp))+
  ggtitle("Comparison of E. coli metabolomics datasets")+
  geom_text_repel(aes(label=as.character(Metabolite) ),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  xlab('Changes in 100met dataset')+
  ylab('Changes in nucleotide dataset')+
  scale_colour_gradientn(colours = gradcols,breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  facet_wrap(~Contrast)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_100met_Nuc_comparison.pdf',sep = ''),
             width=12,height=12,useDingbats=FALSE)








amp<-4
cbrks<-seq(-amp,amp,by=1)
camp<-4
ggplot(subset(all.c,Contrast %in% c("C_T-C_C","MR_T-MR_C","MR_C-C_C","MR_T-MR_C-(C_T-C_C)") ),
       aes(x=Eco_100_logFC,y=Eco_AA_logFC,color=Eco_AA_logFC-Eco_100_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=Eco_100_NE,xmax=Eco_100_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=Eco_AA_NE,ymax=Eco_AA_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_point()+
  scale_x_continuous(breaks=seq(-10,10,by=1),limits = c(-camp,camp))+
  scale_y_continuous(breaks=seq(-10,10,by=1),limits = c(-camp,camp))+
  ggtitle("Comparison of E. coli metabolomics datasets")+
  geom_text_repel(aes(label=as.character(Metabolite) ),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  xlab('Changes in 100met dataset')+
  ylab('Changes in AA dataset')+
  scale_colour_gradientn(colours = gradcols,breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  facet_wrap(~Contrast)
  



amp<-4
cbrks<-seq(-amp,amp,by=1)
camp<-4
ggplot(subset(ec.shared,Contrast %in% maincontrast ),
       aes(x=logFC_eco,y=logFC_cel,color=logFC_cel-logFC_eco))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=NE_eco,xmax=PE_eco),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=NE_cel,ymax=PE_cel),alpha=erralpha,color=errcolor,width=0)+
  geom_point()+
  scale_x_continuous(breaks=seq(-10,10,by=1),limits = c(-camp,camp))+
  scale_y_continuous(breaks=seq(-10,10,by=1),limits = c(-camp,camp))+
  ggtitle("Comparison of E. coli and C. elegans metabolomics datasets")+
  geom_text_repel(aes(label=as.character(Metabolite) ),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  xlab('Changes in E. coli')+
  ylab('Changes in C. elegans')+
  scale_colour_gradientn(colours = gradcols,breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  facet_grid(Contrast~Dataset,labeller=labeller(Dataset=label_both))

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_Eco_Cel_comparison.pdf',sep = ''),
             width=12,height=20,useDingbats=FALSE)











#Venn diagram
library('Vennerable')

EC.AA<-as.character(subset(all,Organism=='Eco' & Dataset=='AA')$Metabolite)
EC.100<-as.character(subset(all,Organism=='Eco' & Dataset=='100')$Metabolite)
Cel.AA<-as.character(subset(all,Organism=='Cel' & Dataset=='AA')$Metabolite)

intersect(EC.AA,EC.100)


setdiff(Cel.AA,EC.100)


setdiff(EC.AA,EC.100)


setdiff(EC.100,EC.AA)

ecoli.sets<-list('E. coli AA'=EC.AA,
             'E. coli 100 metabolites'=EC.100)


ecocel.sets<-list('C. elegans AA'=Cel.AA,
                  'E. coli 100 metabolites'=EC.100)



plot(Vennerable::Venn(ecoli.sets),show = list(Faces = FALSE),doWeights = FALSE)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Ecoli_AA_100met_overlap.pdf",sep=''),
             width=8,height=8, useDingbats=FALSE)




plot(Vennerable::Venn(ecocel.sets),show = list(Faces = FALSE),doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Celegans_Ecoli_100met_overlap.pdf",sep=''),
             width=8,height=8, useDingbats=FALSE)


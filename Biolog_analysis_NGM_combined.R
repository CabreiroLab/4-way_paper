library(tidyverse)
library(ggrepel)

#Uses older version of PFun
#devtools::install_github("PNorvaisas/PFun@v0.2.0")
#devtools::install_github("PNorvaisas/PFun")
library(PFun)


pole<-function(x,y) {
  r<-sqrt(x^2+y^2)
  o<- -2*atan(y/(x))/(pi/2)
  return(o) 
}

pole2<-function(x,y) {
  o<- (atan(y/(x)))/(2*pi)
  return(o) 
}

#theme_set(theme_Publication())
#theme_set(theme_light())


#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau



cwd<-"~/Dropbox/Projects/2015-Metformin/Biolog_Met_NGM/"
setwd(cwd)


odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#load('Biolog_combined.RData')
#save.image('Biolog_combined.RData')



info<-read_csv('../Biolog/Biolog_metabolites_EcoCyc_Unique_PM1-PM5.csv') %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A'))


PM1PM2<-subset(info,Plate %in% c('PM1','PM2A'))

#Pick data to show

logFDRbreaks<-c(-1,1.3,2,3,14)
logFDRbins<-c('N.S.','p<0.05','p<0.01','p<0.001')
logFDRbinsl<-c('N.S.','<0.05','<0.01','<0.001')

logFDRcols<-c("gray40", "red4", "red3",'red')


#Carboxylic acids to remove
coxy<-c('Itaconic Acid','Caproic Acid','Capric Acid','4-Hydroxy Benzoic Acid','2-Hydroxy Benzoic Acid')


selmets<-info %>%
  filter(!(Plate %in% c('PM3B','PM4A') &
             Metabolite %in% PM1PM2$Metabolite ) &
           !Metabolite %in% c('Negative Control',coxy))

#Do analysis on the spot
# data<-read_csv('Data/Summary.csv') %>%
#   filter(Plate %in% c('PM1','PM2A','PM3B','PM4A')) %>%
#   mutate(Sample=paste(Strain,ifelse(Type=='Control','C','T'),sep='_'),
#          SampleID=paste(Sample,Replicate,sep='_'),
#          Sample=factor(Sample,
#                        levels=c("OP50Sens_C","OP50Sens_T"),
#                        labels=c("OP50Sens_C","OP50Sens_T") ) ) %>%
#   left_join(info[,c('Index','MetaboliteU')],by='Index') %>%
#   select(File:Metformin_mM,Replicate,Well,Index,Sample:MetaboliteU,Metabolite=Name,EcoCycID:Group,G=Int_750nm_log,GR=a_log) %>%
#   gather(Measure,Value,G,GR) %>%
#   group_by(Plate,Type,Replicate,Group,Measure) %>%
#   mutate(Value_ref=Value[Metabolite=='Negative Control'],
#          Value_norm=ifelse(Metabolite=="Negative Control",Value,Value-Value_ref)) %>%
#   ungroup %>%
#   mutate_at(c('SampleID','Sample','Strain','Metformin_mM'),as.factor)
# 
# 
# data.nc<-data %>%
#   filter(Plate !='PM5')
# 
# 
# 
# 
# contrasts<-read.contrasts2('!Contrasts.xlsx')
# 
# contrasts.desc<-contrasts$Contrasts.table %>%
#   select(Description:Strain)
# 
# contr.matrix<-contrasts$Contrasts.matrix
# contr.matrix
# 
# 
# results.all<-data.nc %>%
#   #filter(Metabolite !="Negative Control") %>%
#   #filter(!MetaboliteU %in% selmets$MetaboliteU) %>%
#   group_by(Measure,Group,Plate,Well,Index,Metabolite,MetaboliteU,EcoCycID,KEGG_ID) %>%
#   do(hypothesise2(.,'Value_norm~0+Sample',contr.matrix)) %>%
#   getresults(contrasts.desc,c("Measure"))
# 
# 
# #Separate results
# results<-results.all$results %>%
#   filter(Contrast!="T-C_none" ) %>% #
#   group_by(Measure,MetaboliteU) %>%
#   mutate(FDR=p.adjust(p.value,method = "fdr")) %>%
#   ungroup
#   
#   #filter(MetaboliteU %in% selmets$MetaboliteU)
# 
# results %>%
#   group_by(Measure,Contrast) %>%
#   summarise(Total=n(),
#             Ant=sum(logFC>0 & FDR<=0.05,na.rm=TRUE),
#             Syn=sum(logFC<0 & FDR<=0.05,na.rm=TRUE),
#             All=sum(FDR<=0.05,na.rm=TRUE))
# 
# 
# 
# results.eco<-results




results.eco<-read_csv('Summary/Ecoli_results.csv') 

metorder<-results.eco %>%
  filter(MetaboliteU %in% selmets$MetaboliteU & Contrast=="T-C" & Measure=="G") %>%
  arrange(logFC) %>%
  mutate(MetaboliteU=factor(MetaboliteU,levels=MetaboliteU)) %>%
  pull(MetaboliteU) %>%
  as.character


results.cel<-read_csv('Celegans/Summary/Celegans_results_withNGM.csv') %>%
  filter(! Index %in% c('Controls-A1','Controls-B1')) %>%
  select(Index,logFC:logFDR_bin) %>%
  left_join(info) %>%
  mutate(Organism='C. elegans',
         Measure="F",
         Description="C. elegans rescue via GFP::Pacs-2 fluorescence reduction",
         Contrast_type="Treatment",
         Contrast="T-C",
         logFCrev=-logFC,
         PErev=-NE,
         NErev=-PE) %>%
  select(Organism,Measure,Contrast,Contrast_type,Description,Plate,Well,Index,MetaboliteU,EcoCycID,KEGG_ID,logFC=logFCrev,SE,PE=PErev,NE=NErev,FDR,logFDR)


results.ecocel<-results.eco %>%
  mutate(Organism='E. coli') %>%
  select(Organism,Measure,Contrast,Contrast_type,Description,Plate,Well,Index,MetaboliteU,EcoCycID,KEGG_ID,logFC,SE,PE,NE,FDR,logFDR) %>%
  rbind(results.cel) %>%
  filter(MetaboliteU %in% selmets$MetaboliteU) %>%
  mutate(MetaboliteU=factor(MetaboliteU,levels=metorder,labels=metorder),
         FDRStars=pStars(FDR),
         logFDR_bin=cut(logFDR, breaks=logFDRbreaks,labels=logFDRbinsl),
         Organism=factor(Organism,levels=c('E. coli','C. elegans') ),
         Organism_short=str_replace(Organism,". ","") %>% str_sub(1,2)) %>%
  unite(OMC,Organism_short,Measure,Contrast,remove = FALSE) %>%
  group_by(Organism,Measure,Contrast) %>%
  arrange(logFC) %>%
  ungroup




ecocelmulti2<-multiplex(results.ecocel,c("Plate","Well","Index","MetaboliteU","EcoCycID","KEGG_ID"),2)
ecocelmulti<-multiplex(results.ecocel,c("Plate","Well","Index","MetaboliteU","EcoCycID","KEGG_ID"),3)
ecocelmulti4<-multiplex(results.ecocel,c("Plate","Well","Index","MetaboliteU","EcoCycID","KEGG_ID"),4)




results.ecocel %>%
  filter(Measure %in% c('G',"F") & Contrast=="T-C" & Organism=="E. coli") %>%
  arrange(desc(logFC))


results.ecocel %>%
  group_by(Organism,Measure,Contrast) %>%
  summarise(Total=n(),
            Up=sum(logFC>0 & FDR<0.05,na.rm=TRUE),
            Down=sum(logFC<0 & FDR<0.05,na.rm=TRUE),
            All=sum(FDR<0.05,na.rm=TRUE))


TallPlot<-function(data){
  data %>%
  ggplot(aes(y=MetaboliteU,x=logFC,color=logFDR_bin))+
    geom_vline(xintercept = 0,color='red',alpha=0.5)+
    geom_errorbarh(aes(xmin=NE,xmax=PE),height=0,size=0.5)+
    geom_point(size=0.5)+
    ylab('Metabolites')+
    labs(color='FDR')+
    scale_colour_manual(values = logFDRcols)+
    scale_x_continuous(breaks=seq(-10,10,by=2) )+
    theme(axis.ticks.y=element_blank(),
          axis.line = element_line(colour = NA),
          axis.line.x = element_line(colour = NA),
          axis.line.y = element_line(colour = NA),
          #panel.border=element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())
}


results.ecocel %>%
  filter(OMC %in% c("Ec_G_T-C","Ce_F_T-C") ) %>%
  TallPlot+
  facet_grid(~Organism)

dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/EcoliCelegans_logFC_tall_Treatment.pdf",sep=''),
             width=9,height=50)


results.ecocel %>%
  filter(OMC %in% c("Ec_G_T-C","Ce_F_T-C")) %>%
  TallPlot+
  facet_grid(~Organism)+
  theme( axis.text.y=element_blank())

ggsave(file=paste0(odir,"/EcoliCelegans_logFC_tall_Treatment_tiny.pdf"),
             width=60,height=60,units='mm',scale=2,device=cairo_pdf,family="Arial")




results.ecocel %>%
  filter(OMC %in% c("Ec_G_T-C")) %>%
  TallPlot+
  xlab('Metabolite - metformin interaction as growth logFC vs NGM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ecoli_logFC_tall_Treatment.pdf",sep=''),
             width=9,height=50)



results.ecocel %>%
  filter(OMC %in% c("Ce_F_T-C")) %>%
  arrange(logFC)%>%
  mutate(MetaboliteU=factor(MetaboliteU,levels=MetaboliteU)) %>%
  TallPlot+
  scale_x_continuous(breaks=-10:10)+
  xlab('Metabolite - metformin interaction as growth logFC vs NGM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Celegans_logFC_tall_Fluorescence.pdf",sep=''),
             width=9,height=50)





head(selectcast)
#Scatter plots

lblsize<-4

erralpha<-1
errcolor<-'grey80'
segalpha=0.5


grp<-'Complete'

fit<-ecocelmulti %>%
  filter(x_OMC=="Ec_G_C",
         y_OMC=="Ec_G_T",
         z_OMC=="Ce_F_T-C") %>%
  lm(y_logFC~x_logFC,data=.)


a<-fit$coefficients[[2]]
b<-fit$coefficients[[1]]

maincomp<-'C. elegans\nphenotype rescue\n(acs-2 GFP)'


amp<-2
cbrks<-seq(-amp,amp,by=1)
gradcols<-c('blue4','blue','gray80','red','red4')


thres<-0.05

gradcols<-c('blue4','blue','gray80','red','red4')



showmets<-c("D-Arabinose","L-Arabinose","Acetoacetic Acid","D-Ribose","Glycerol","alpha-D-Glucose","Phosphono Acetic Acid","")

ecocelmulti %>%
  filter(x_OMC=="Ec_G_C",
         y_OMC=="Ec_G_T",
         z_OMC=="Ce_F_T-C") %>%
  ggplot(aes(y=y_logFC,x=x_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(intercept=1,slope=1),color='blue',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  #ggtitle('Scatterplot of metformin and metabolite supplementation effects')+
  scale_x_continuous(breaks=-10:10) +
  scale_y_continuous(breaks=-10:10) +
  geom_point(aes(color=z_logFC),size=3)+
  coord_cartesian(xlim=c(-2,3),ylim = c(-2,3.5))+
  xlab('Growth logFC vs NGM - Control')+
  ylab('Growth logFC vs NGM - +50mM Metformin')+
  labs(color='C. elegans\nphenotype rescue\n(Pacs-2::GFP)')+
  geom_text_repel(aes(label=ifelse(MetaboliteU %in% showmets, as.character(MetaboliteU),"" ),color=z_logFC)  )+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(file=paste0(odir,"/Scatter_Control-Treatment_Complete_Selmets.pdf"),
             width=70,height=70,units='mm',scale=2,device=cairo_pdf,family="Arial")


#C elegans, E coli comparison - pole 

xmin<- -3
xmax<-3
ymin<- -2
ymax<- 2

errcolor<-"grey90"
sbrks<-seq(0,3,by=1)

cbrks<-seq(-1,1,by=1)
gradcols<-c('gray80','blue3','gray80','red3','gray80')

grp<-'Complete'
ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T-C", 
         y_OMC=="Ce_F_T-C") %>%
  ggplot(aes(x=x_logFC,y=y_logFC,color=pole(-y_logFC,x_logFC)))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_errorbar(aes(ymin=y_PE,ymax=y_NE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=sqrt( y_logFC^2 + x_logFC^2)))+
  coord_cartesian(xlim=c(xmin,xmax),ylim = c(ymin,ymax))+
  scale_x_continuous(breaks=seq(xmin,xmax,by=1))+
  scale_y_continuous(breaks=seq(ymin,ymax,by=1))+
  xlab('Growth rescue in E. coli (metabolite effect normalised), logFC')+
  ylab('Phenotype rescue in C. elegans (acs-2 GFP), logFC')+
  geom_text_repel(aes(label=ifelse(y_FDR<0.05 &  x_FDR<0.05, as.character(MetaboliteU),'')),
                  size=lblsize,nudge_y = 0.3,force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-2,2))+
  scale_size(range = c(0.25, 7),breaks=sbrks,name='Strength of effect')+
  labs(color='Consistency of effect')+
  ggtitle(paste('C. elegans and E. coli phenotype rescue in metformin treatment by metabolites: ',grp,sep='') )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Celegans_Ecoli_effect_correlation_Complete_BothSignificant.pdf",sep=''),
             width=13,height=9, useDingbats=FALSE)




#Driving effects for interaction


maincomp<-'C. elegans\nphenotype rescue\n(acs-2 GFP)'


xmin<- -3
xmax<- 4
ymin<- -3
ymax<- 4


lblsize<-4

erralpha<-1
errcolor<-'grey80'
segalpha=0.5

amp<-2
cbrks<-seq(-amp,amp,by=1)
gradcols<-c('blue4','blue','gray80','red','red4')

#Needs error checking
ecocelmulti4 %>%
  filter(x_OMC=="Ec_G_C", 
         y_OMC=="Ec_G_T",
         z_OMC=="Ec_G_T-C",
         w_OMC=="Ce_F_T-C") %>%
  mutate(Pole=ifelse(x_logFC>0,
                   pole2(x_logFC,y_logFC),
                   pole2(x_logFC,y_logFC)+0.5),
       Pole=ifelse(Pole>0.625,Pole-1,Pole),
       Pole360=Pole*360) %>%
  ggplot(aes(y=y_logFC,x=Pole360,color=w_logFC))+
  geom_point()+
  geom_point(aes(size=abs(w_logFC)))+
  geom_text_repel(aes(label=ifelse(w_FDR<0.05 & z_FDR<0.05, as.character(MetaboliteU),'')),
                  size=lblsize,nudge_y = 0.3,
                  force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)+
  ggtitle('Breakdown of driving components for metabolite-metformin interaction')+
  #scale_x_discrete(labels=c("-90"="-T","-45"="-TC","0" = "C", "45" = "TC","90" = "T","135"="T-C","180"="-C","225"="-T-C","270"="-T"))+
  scale_x_continuous(breaks=seq(-135,225,by=45),labels=c("-T-C","-T","-TC","C","TC","T","T-C","-C","-T-C")  )+
  scale_y_continuous(breaks=seq(-4,4,by=1),limits = c(-3,3) )+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp))+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  xlab('Driving component')+
  ylab('Metabolite - metformin interaction as growth logFC vs NGM')+
  scale_size(range = c(0.25, 7),name=maincomp)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Celegans_Ecoli_Driving_effects_for_interaction_CelSignificant.pdf",sep=''),
             width=16,height=12, useDingbats=FALSE)

#Volcano plots


ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T-C",
         y_OMC=="Ce_F_T-C") %>%
  ggplot(aes(x=x_logFC,y=x_logFDR,color=y_logFC))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=abs(y_logFC) ))+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  scale_size(range = c(0.25, 7),name=maincomp)+
  ylab('-logFDR')+
  xlab('Metabolite - metformin interaction as growth logFC vs NGM')+
  labs(color='Significance (FDR)')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Volcano_Ecoli_Growth.pdf",sep=''),
             width=12,height=9)



ecocelmulti2 %>%
  filter(x_OMC=="Ec_GR_T-C",
         y_OMC=="Ce_F_T-C") %>%
  ggplot(aes(x=x_logFC,y=x_logFDR,color=y_logFC))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=abs(y_logFC) ))+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  scale_size(range = c(0.25, 7),name=maincomp)+
  ylab('-logFDR')+
  xlab('Doubling rate difference (Treatment-Control), log2(OD)/h')+
  labs(color='Significance (FDR)')


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Volcano_Ecoli_GrowthRate.pdf",sep=''),
             width=12,height=9)




#Growth rate analysis
fit<-ecocelmulti2 %>%
  filter(x_OMC=="Ec_GR_C",
         y_OMC=="Ec_GR_T") %>%
  lm(y_logFC~x_logFC,data=.)


a<-fit$coefficients[[2]]
b<-fit$coefficients[[1]]

thres<-0.05

grp<-'All'

ecocelmulti %>%
  filter(x_OMC=="Ec_GR_C",
         y_OMC=="Ec_GR_T",
         z_OMC=="Ce_F_T-C") %>%
  ggplot(aes(y=y_logFC,x=x_logFC,color=z_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(intercept=b,slope=a),alpha=1,color='red')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  ggtitle(paste('Scatterplot of metformin and metabolite supplementation effects ',grp,sep=''),
          subtitle = paste('Metabolites with FDR<',thres,' are marked',sep='') )+
  geom_point(size=3)+#aes(size=abs(Cel_logFC))
  #coord_cartesian(xlim=c(xmin,xmax),ylim = c(ymin,ymax))+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  
  xlab('Doubling rate - Control, log2(OD)/h')+
  ylab('Doubling rate - +50mM Metformin, log2(OD)/h')+
  #eval(parse(text = intfdr)) < thres & abs( eval(parse(text = intvar)) ) > 0.75 
  # geom_text_repel(aes(label=ifelse(Cel_FDR<thres, as.character(Metabolite),'')),
  #                 size=lblsize,nudge_y = 0.3,
  #                 force=1,
  #                 segment.colour=errcolor,
  #                 segment.alpha =segalpha)+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp))+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  #scale_size(range = c(0.25, 7),name=maincomp)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Scatter_Control-Treatment_Complete_GrowthRate_Celegans.pdf",sep=''),
             width=13,height=9, useDingbats=FALSE)


#Growth rate and growth
ecocelmulti %>%
  filter(x_OMC=="Ec_G_T-C",
         y_OMC=="Ec_GR_T-C",
         z_OMC=="Ce_F_T-C") %>%
  ggplot(aes(x=x_logFC,y=y_logFC,color=z_logFC) )+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  #geom_abline(aes(intercept=b,slope=a),alpha=1,color='red')+
  #geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=abs(z_logFC)))+
  geom_text_repel(aes(label=ifelse(z_FDR<thres, as.character(MetaboliteU),'')),
                  size=lblsize,nudge_y = 0.3,
                  force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)+
  ylab('Doubling rate difference (Treatment-Control), log2(OD)/h')+
  xlab('Metabolite - metformin interaction as growth logFC vs NGM')+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp))+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  scale_size(range = c(0.25, 7),name=maincomp)+
  coord_cartesian(xlim=c(-3,3),ylim = c(-2,2))


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Scatter_GrowthRate_vs_Growth_vs_Celegans.pdf",sep=''),
             width=13,height=9, useDingbats=FALSE)




#For summaries
#Time series
ectimem<-read_csv('Data/Data.csv') %>%
  filter(Data=='750nm_f') %>%
  gather(Time_sec,OD,ends_with(".0")) %>%
  mutate(Time_sec=as.numeric(as.character(Time_sec)),
         Time_min=Time_sec/60,
         Time_h=Time_min/60)

tsum<-ectimem %>%
  group_by(Index,Type,Time_sec,Time_h) %>%
  summarise(Mean=mean(OD),SD=sd(OD),SE=SD/sqrt(length(OD))) %>%
  left_join(info)


results.in<-read_csv('Celegans/Summary/Celegans_results_withNGM.csv')
data.Isum<-readRDS('Celegans/Celegans_fluorescence_distributions.rds')


indxsel<-c('PM1-A1','PM1-A2','PM4A-E6','PM1-G7') #'Controls-A1',

indxsel<-c('Controls-A1','PM1-A1','PM1-A2','PM1-G7','PM4A-E6')

mlevels<-c('Negative Control','L-Arabinose','Acetoacetic Acid','Phosphono Acetic Acid')
mlabels<-c('NGM','L-Arabinose','Acetoacetic Acid','Phosphono Acetic Acid')




#Distributions

Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- c("Control","Treatment")
Metlab<-'Metformin, mM'


tsum %>%
  filter(Index %in% indxsel) %>%
  mutate( Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels)) %>%
  ggplot(aes(x=Time_h,y=Mean,fill=Type,color=Type))+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),color=NA,alpha=0.2)+
  scale_colour_manual(name = Metlab,values =Metcols)+
  scale_fill_manual(name = Metlab,values =Metcols)+
  geom_line()+
  scale_x_continuous(breaks=seq(0,24,by=6))+
  ylab("OD")+
  xlab("Time, h")+
  labs(fill="Type")+
  facet_grid(~Metabolite)+
  theme(legend.position = "none")

ggsave(file=paste0(odir,"/Summary_Ecoli_growth_curves.pdf"),
       width=80,height=26,units='mm',scale=2,device=cairo_pdf,family="Arial")



data.Isel<-data.Isum %>%
  filter(Index %in% indxsel & Measure=='Absolute') %>%
  mutate(Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels))

res.sel<-results.in %>%
  filter(Index %in% indxsel) %>%
  mutate(Type=as.factor(ifelse(Index=='Controls-A1','Control','Treatment')),
         Metabolite=as.character(Metabolite),
         Metabolite=ifelse(Index=='Controls-A1','Negative Control',Metabolite),
         Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels),
         Metformin_mM=factor(Type,levels=c("Control","Treatment"),labels=c("0","50")) )


ggplot(res.sel,aes(color=Type))+
  geom_vline(aes(xintercept=res.sel[res.sel$Index=='PM1-A1','Raw_Q90_Mean']),color='grey80',linetype='longdash')+
  geom_rect(aes(xmin=Raw_Q90_Mean-Raw_Q90_SE,xmax=Raw_Q90_Mean+Raw_Q90_SE,fill=Type),ymin=-Inf,ymax=Inf,alpha=0.2,color=NA)+
  geom_vline(aes(xintercept=Raw_Q90_Mean,color=Type))+
  geom_ribbon(data=data.Isel,aes(x=LogBrightness_num,ymin=(Mean-SD)*100,ymax=(Mean+SD)*100,fill=Type),alpha=0.5,color=NA)+
  geom_line(data=data.Isel,aes(x=LogBrightness_num,y=Mean*100,color=Type))+
  scale_colour_manual(name = Metlab,values =Metcols)+
  scale_fill_manual(name = Metlab,values =Metcols)+
  ylab('Probability density, %')+
  xlab('log2 Brightness')+
  scale_x_continuous(breaks=seq(-20,20,by=1))+
  coord_cartesian(ylim=c(0,15))+
  facet_grid(.~Metabolite)+
  theme(legend.position = "none")

ggsave(file=paste(odir,"/Summary_Celegans_brightness_distribution.pdf",sep=''),
             width=80,height=26,units='mm',scale=2,device=cairo_pdf,family="Arial")


#Comparisons

res.Ecsel<-read_csv(paste0(odir,"/Ecoli_results_for_examples_with_NC.csv")) %>%
  filter(Contrast %in% c("C","T")  & Measure=="G")  %>%
  group_by(Contrast,Group) %>%
  mutate(logFCnew=ifelse(Metabolite!="Negative Control",logFC+logFC[Metabolite=="Negative Control"],logFC),
         NEnew=ifelse(Metabolite!="Negative Control",NE+logFC[Metabolite=="Negative Control"],NE),
         PEnew=ifelse(Metabolite!="Negative Control",PE+logFC[Metabolite=="Negative Control"],PE) ) %>%
  ungroup %>%
  filter(Index %in% indxsel) %>%
  mutate(Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels),
         Metformin_mM=factor(Contrast,levels=c("C","T"),labels=c("0","50")),
         Organism="E. coli") %>%
  select(Metabolite,Metformin_mM,logFC=logFCnew,PE=PEnew,NE=NEnew,SE=SE,Organism) 
  



res.Ecsel %>%
  ggplot(aes(x=Metabolite,y=logFCnew,color=Metformin_mM))+
  geom_hline(aes(yintercept=res.Ecsel[res.Ecsel$Index=='PM1-A1' & res.Ecsel$Contrast=='C','logFC']),color='grey80',linetype='longdash')+
  #scale_y_continuous(breaks=seq(-20,20,by=1),limits=c(-4,0))+
  geom_errorbar(aes(ymin=NEnew,ymax=PEnew),position = position_dodge(width = dwidth),width=0.25)+
  geom_point(position = position_dodge(width = dwidth),size=2)+
  labs(color="Metformin, mM")+
  ylab('log2 AUC, OD*h')+
  theme(axis.text.x = element_text(angle = 90,hjust=1))

dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Summary_Ecoli_growth_comparison.pdf",sep=''),
             width=4,height=4)


dwidth<-0.5
res.sel %>%
  ggplot(aes(x=Metabolite,y=Raw_Q90_Mean,color=Metformin_mM))+
  geom_hline(aes(yintercept=res.sel[res.sel$Index=='PM1-A1','Raw_Q90_Mean']),color='grey80',linetype='longdash')+
  scale_y_continuous(breaks=seq(-20,20,by=1),limits=c(-4,0))+
  geom_errorbar(aes(ymin=Raw_Q90_Mean-Raw_Q90_SE,ymax=Raw_Q90_Mean+Raw_Q90_SE),position = position_dodge(width = dwidth),width=0.25)+
  geom_point(position = position_dodge(width = dwidth),size=2)+
  ylab('log2 brightness Q90')+
  labs(color="Metformin, mM")+
  theme(axis.text.x = element_text(angle = 90,hjust=1))



dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Summary_Celegans_brightness_comparison.pdf",sep=''),
             width=4,height=4)




reflines<-data.frame(Organism=c("E. coli","C. elegans"),
                     Ref=c(2.03,-0.989))

res.comb<-res.sel %>%
  select(Metformin_mM,Metabolite,logFC=Raw_Q90_Mean,SE=Raw_Q90_SE) %>%
  mutate(NE=logFC-SE,PE=logFC+SE,
         Organism="C. elegans") %>%
  rbind(res.Ecsel) %>%
  mutate(Organism=factor(Organism,levels=c("E. coli","C. elegans")))


Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- c("0","50")
Metlab<-'Metformin, mM'


res.comb %>%
  ggplot(aes(x=Metabolite,y=logFC,color=Metformin_mM))+
  geom_hline(data=reflines,aes(yintercept=Ref),color='grey80',linetype='longdash')+
  scale_y_continuous(breaks=seq(-20,20,by=1))+
  geom_errorbar(aes(ymin=NE,ymax=PE),position = position_dodge(width = dwidth),width=0.25)+
  scale_colour_manual(name = Metlab,values =Metcols)+
  geom_point(position = position_dodge(width = dwidth),size=1)+
  labs(color="Metformin, mM")+
  #ylab('log2 AUC, OD*h')+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  facet_grid(Organism~.,scale="free_y")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(file=paste0(odir,"/Summary_CelegansEcoli_comparison.pdf"),
             width=30,height=60,units='mm',scale=2,device=cairo_pdf,family="Arial")






#Venn diagram
results.ecocel %>%
  filter(MetaboliteU=="D-Glucosaminic Acid" & Contrast=="T-C")


venn.pass<-results.ecocel %>%
  mutate(Ant=logFC>0 & FDR<0.05,
         Syn=logFC<0 & FDR<0.05,
         All=FDR<0.05) %>%
  #gather different thresholds
  gather(Type,Pass,Ant,Syn,All) %>%
  filter(Pass) %>%
  group_by(Organism,Organism_short,OMC,Measure,Contrast,Type) %>%
  do(List=c(as.character(.$MetaboliteU)))


rescue<-venn.pass %>%
  filter(OMC %in% c("Ec_G_T-C","Ce_F_T-C") & Type!="All") %>%
  unite(OT,Organism_short,Type) %>%
  select(OT,List) %>%
  spread(OT,List) %>%
  as.list %>%
  unlist(recursive=FALSE)



names(rescue)


intersect(rescue[["Ce_Ant"]],rescue[["Ec_Syn"]])

intersect(rescue[["Ce_Syn"]],rescue[["Ec_Ant"]])

intersect(rescue[["Ce_Syn"]],rescue[["Ec_Syn"]])


setdiff(Cel.rescue,EC.ant)





vcols<-c("red","red4","blue","blue4")
grid::grid.draw(VennDiagram::venn.diagram(rescue, col=vcols,cat.col=vcols ,NULL))


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Celegans_Ecoli_updated.pdf",sep=''),
             width=3,height=3, useDingbats=FALSE)
dev.off()#





#EcoCyc class enrichment

ecli.co<-read_csv('../Biolog/Ecoli_All_class_instances_non_redundant.csv') %>%
  select(-X1)

ecli.co %>%
  group_by(Class) %>%
  summarise()


#EcoCyc to Biolog links
pmcs<-read_csv('../Biolog/EColi_net/All_Media_Mix_clean_explicit.csv')

enrbrks<-c(0,-log(0.05,10),2,3,4,100)
enrlbls<-c('N.S.','<0.05','<0.01','<0.001','<0.0001')

EcoCyc.enrichment<-results.ecocel %>%
  left_join(pmcs %>% select(Plate,Well,Instance)) %>%
  left_join(ecli.co %>% filter(NoInstances>2)) %>%
  filter(!is.na(FDR) & !is.na(Class)) %>%
  enrichment(terms = c('Class','NoInstances','subclasses','supclasses','Hierarchy'),
             IDs = 'Instance',
             comparisons = c('Organism','Organism_short','Measure','Contrast','OMC'),
             change = 'logFC',sign = 'FDR') %>%
  mutate(logFDR=ifelse(-log10(FDR)<0,0,-log10(FDR)),
         logFDRbin=cut(logFDR,breaks=enrbrks,labels=enrlbls,right=FALSE),
         Test=factor(Test,levels=c("All","Up","Down")),
         Class=str_replace_all(as.character(Class),"\\|","")) 
  

EcoCyc.enrichment %>%
  filter(Test!="All" & OMC %in% c("Ec_G_T-C","Ce_F_T-C") & !Class %in% EChm$Class)%>%
  select(-c(NoInstances:Comparison),-c(Organism_short,Measure,Contrast,OMC) ) %>%
  write_csv(paste0(odir,'/EcoCyc_enrichment.csv'))


ec.heatmap<-EcoCyc.enrichment %>%
  filter(OMC %in% c("Ec_G_T-C","Ce_F_T-C") & Test!="All") %>%
  group_by(Class) %>%
  filter( any(FDR<0.05)) %>%
  ungroup %>%
  unite(OT,Organism,Test) %>%
  select(Class,OT,FDR) %>%
  spread(OT,FDR) 

ec.heatmap %>%
  write_csv(paste0(odir,'/EcoCyc_enrichment_heatmap_all_classes.csv'))
  
  

#Some manual ordering
EChm<-readxl::read_xlsx(paste0(odir,'/EcoCyc_enrichment_heatmap.xlsx'),sheet = 'Heatmap') %>%
  filter(is.na(Remove)) %>%
  mutate(Class=str_replace_all(as.character(Class),"\\|",""))


# classorder<-c("Carbohydrates","Aldehydes-Or-Ketones",
#               "Carboxylates","Amino-Sugars","Peptides","Alcohols",
#               "All-Nucleosides","Ribonucleosides","Purines",
#               "Alpha-Amino-Acids","Amino-Acids","L-Amino-Acids","dicarboxylate",
#               "Esters","Lactones","Glycosides")

EcoCyc.enrichment %>%
  filter(Test!="All" & Class %in% EChm$Class & OMC %in% c("Ec_G_T-C","Ce_F_T-C") ) %>% 
  group_by(Class) %>%
  filter( any(FDR<0.05)) %>%
  ungroup %>%
  clustorder('Class',c("Organism","Test"),'logFDR',descending=FALSE) %>% 
  #mutate(Class=factor(Class,levels=rev(classorder) )) %>%
  PlotEnrichment("Test","Class",ncols = length(enrbrks))+
  facet_grid(~Organism)+
  ylab("EcoCyc metabolite class")

# dev.copy2pdf(device=cairo_pdf,file=paste0(odir,'/Enrichment_heatmap_EcoCyc_UpDown.pdf'),
#              width=5,height=5,useDingbats=FALSE)




# groups<-c("Organism","Organism_short","Measure","Contrast","OMC")
# features<-c("Class","Hierarchy","NoInstances","subclasses","supclasses")
# 
# results.ecocel %>%
#   left_join(pmcs %>% select(Plate,Well,Instance)) %>%
#   right_join(ecli.co %>% filter(NoInstances>2)) %>%
#   filter(!is.na(FDR) & !is.na(Class) ) %>%
#   enrichment(groups,features) %>%
#   filter(Type!="All" & OMC %in% c("Ec_G_T-C","Ce_F_T-C") ) %>%
#   group_by(Class) %>%
#   filter( any(FDR<0.05)) %>%
#   ungroup %>%
#   clustorder('Class',c("Organism","Type"),'logFDR',descending=TRUE) %>%
#   #mutate(Class=factor(Class,levels=rev(classorder) )) %>%
#   PlotEnrichment("Type","Class",ncols=4)+
#   facet_grid(~Organism)+
#   ylab("EcoCyc metabolite class")

#Biolog class enrichment
#Classes preparation
classes<-info %>%
  select('Index','Plate','Well','MetaboliteU') %>%
  left_join(readxl::read_xlsx('../Biolog/Biolog_metabolites_EcoCyc.xlsx',sheet = 'Classes') %>% select('Index','Class1','Class2')) %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A')) %>%
  gather(Redundancy,Class,contains("Class")) %>%
  filter(Class!="") 



#phyper(q, m, n, k)
# pop size : 5260
# sample size : 131
# Number of items in the pop that are classified as successes : 1998
# Number of items in the sample that are classified as successes : 62
#phyper(62-1, 1998, 5260-1998, 131, lower.tail=FALSE)




# groups<-c("Organism","Organism_short","Measure","Contrast","OMC")
# feature<-"Class"

class.enrichment<-results.ecocel%>%
  left_join(classes) %>%
  filter(!is.na(FDR)) %>%
  #enrichment(groups,feature)
  enrichment(terms = c('Class'),
             IDs = 'MetaboliteU',
             comparisons = c('Organism','Organism_short','Measure','Contrast','OMC'),
             change = 'logFC',sign = 'FDR') %>%

  
  

#select thresholds

class.enrichment %>%
  filter(Contrast=="T-C" & Measure!="GR" & p.value<0.05) %>%
  View()


class.enrichment %>% 
  write_csv(paste0(odir,'/Metabolite_class_enrichment_nocarboxy_complete.csv'))

class.enrichment %>%
  filter(OMC %in% c("Ec_G_T-C","Ce_F_T-C") & Type!="All") %>%
  group_by(Class) %>%
  filter(any(FDR<0.05)) %>%
  ungroup %>%
  PlotEnrichment("Type","Class","logFDRbin")+
  facet_grid(~Organism)+
  xlab("Direction of change")+
  ylab('Biolog metbolite class')



#EcoCyc pathway enrichment

#EcoCyc to Biolog links
#Rocver possible links via classes
pmcs2<-read_csv('../Biolog/EColi_net/All_Media_Mix_clean_explicit.csv') %>%
  select(Plate,Well,EcoCycID=Instance)


ecopats<-read_csv("~/Dropbox/Projects/2015-Metformin/Annotations/Ecoli/EcoCyc_pathways_metabolites_2018-04-04.csv") %>%
  select(-X1,PID=`Pathway-ID`,Pathway=`Pathway-Name`,EcoCycID=Compounds,-Comment) %>%
  separate_rows(EcoCycID,sep = ";")


groups<-c("Organism","Measure","Contrast","Contrast_type","Description")


results.ecocel %>%
  filter(Contrast=="T-C") %>%
  mutate(Up=logFC>0 & FDR <0.05,
         Down=logFC<0 & FDR <0.05,
         All=FDR<0.05) %>%
  #gather different thresholds
  gather(Type,Pass,Up,Down,All) %>%
  group_by(Organism,OMC,Measure,Contrast,Contrast_type,Description,Type) %>%
  summarise(Unique_size=n(),
            Unique_pass=sum(Pass,na.rm = TRUE) ) 


EcoCycenr<-results.ecocel %>%
  filter(Contrast=="T-C" & Measure!="GR") %>%
  select(-EcoCycID) %>%
  left_join(pmcs2) %>%
  left_join(ecopats) %>%
  filter(!is.na(PID) & ! is.na(FDR)) %>%
  group_by(Organism,OMC,Measure,Contrast,Contrast_type,Description) %>%
  enrichment(c("PID","Pathway"),featureid="EcoCycID",enrtype="regular")


#phyper(Class_pass-1,Total_pass,Total_size-Total_pass,Class_size,lower.tail=FALSE)

EcoCycenr %>%
  group_by(Organism,Measure,Contrast,Description,Type) %>%
  summarise(Count=n())


EcoCycenr %>%
  filter(Measure %in% c("G","F") & Contrast=="T-C" & p.value<0.05) %>%
  View

# EcoCycenr %>%
#   filter(Measure %in% c("G","F") & Contrast=="T-C" ) %>%
#   write_csv(paste0(odir,"/EcoCyc_pathway_enrichment.csv"))



#KEGG pathway enrichment
#KEGG mapping

#Only unique values with preference for PM1 PM2A

#KEGG info
library(limma)
library(org.EcK12.eg.db)

#allpaths<-c('01110')


path2pathde<-getKEGGPathwayNames('eco',remove.qualifier = TRUE)
path2pathde$PathwayID<-gsub('path:','',path2pathde$PathwayID)

allpaths<-gsub('eco','',path2pathde$PathwayID)

write.csv(path2pathde,paste(odir,'/KEGG_pathways.csv',sep=''))



keggprep<-ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T-C" & y_OMC=="Ce_F_T-C") %>%
  filter(!KEGG_ID %in% c("") &
           !is.na(KEGG_ID)&  !( Plate=='PM4A' & KEGG_ID=='C00097')) %>%
  data.frame



mapkey<-'KEGG_ID'
nrow(keggprep)
gdata<-keggprep
#gdata<-all.results.rcp
nrow(gdata)

dupl<-duplicated(gdata[,"KEGG_ID"])
print('KEGG ID duplicates')
print(table(dupl))
gdataf<-gdata[!dupl,]

rownames(gdataf)<-gdataf[,"KEGG_ID"]




# allpaths<-read_csv(paste0(odir,'/KEGG_pathways.csv')) %>%
#   pull(PathwayID) %>%
#   unique

print(paste('Total KEGG pathways to plot:',length(allpaths)))


#comp.gr<-'Main'

library(pathview)



pathcomp<-c('x_logFC',"y_logFC")
keggxml<-'~/Dropbox/Projects/2015-Metformin/Annotations/Ecoli/KEGG_pathways'


keggdir<-paste(cwd,odir,'/KEGG',sep='')
dir.create(keggdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
setwd(keggdir)

limitcpd<- 2

pv.out <- pathview::pathview(cpd.data = gdataf[,pathcomp,drop=FALSE],
                   pathway.id = allpaths,
                   cpd.idtype = "kegg",
                   species = "eco",
                   out.suffix = 'Ecoli_Celegans',
                   kegg.dir = keggxml,
                   same.layer=FALSE,kegg.native = T,map.symbol=FALSE,
                   map.null=FALSE,new.signature=FALSE, plot.col.key=TRUE,
                   res=300,
                   min.nnodes = 0,
                   limit=list(gene=2,cpd=limitcpd),node.sum = "mean",
                   low = list(gene = "blue", cpd = "blue"),
                   mid=list(gene = "gray", cpd = "gray"),
                   high = list(gene = "red", cpd ="red"))
setwd(cwd)

names(pv.out)

detach("package:pathview", unload=TRUE)
detach("package:limma", unload=TRUE)
detach("package:org.EcK12.eg.db", unload=TRUE)



extractpat<-function(pat,pv.out){
  frame<-data.frame(pv.out[[pat]][['plot.data.cpd']]) 
  if ("all.mapped" %in% colnames(frame)) {
    frame<-frame%>%
      unnest(all.mapped) %>%
      unnest(type) %>%
      unnest(kegg.names) %>%
      filter(!all.mapped %in% c("",NA) ) 
    return(frame)
  } else {
    return(data.frame(NA))
  }
}

all.kegg.mappings<-path2pathde %>%
  left_join( data.frame(PathwayID=names(pv.out)) ) %>%
  group_by(PathwayID,Description) %>%
  do( extractpat(as.character(.$PathwayID),pv.out)  )
  

write.csv(all.kegg.mappings,paste(odir,'/All_KEGG_mappings_Complete.csv',sep = ''),row.names = FALSE)

#all.kegg.mappings<-read_csv(paste0(odir,"/All_KEGG_mappings_Complete.csv"))
  
KEGGmets<-all.kegg.mappings %>%
  filter(all.mapped!="") %>%
  group_by(PathwayID,Description,all.mapped) %>%
  summarise %>%
  rename(KEGG_ID=all.mapped,Pathway=Description)


KEGGgroups<-c("PathwayID","Pathway")

KEGGenrich<-results.ecocel %>%
  separate_rows(KEGG_ID,sep=",") %>%
  left_join(KEGGmets) %>%
  filter(!is.na(KEGG_ID) & !is.na(PathwayID) & !is.na(FDR) & !PathwayID %in% c("eco01502","eco00521")) %>%
  group_by(Organism,OMC,Measure,Contrast,Contrast_type,Description) %>%
  enrichment(KEGGgroups,"KEGG_ID",enrtype="regular")


KEGGenrich %>%
  filter(Contrast=="T-C" & Measure %in% c("G","F")) %>%
  filter(FDR<0.05) %>%
  View

KEGGenrich %>%
  filter(Contrast=="T-C" & Type %in% c("Up","Down") & Measure %in% c("G","F")) %>%
  group_by(PathwayID) %>%
  filter( any(FDR<0.05)) %>%
  ungroup %>%
  clustorder('Pathway',c("Organism","Type"),'logFDR',descending=TRUE) %>%
  PlotEnrichment("Type","Pathway",fillval = "logFDRbin")+
  facet_grid(.~Organism)+
  ylab("KEGG pathway")+
  labs(fill="FDR")+
  theme(axis.text.x = element_text(angle=45))


ggsave(file=paste(odir,'/Enrichment_heatmap_KEGG.pdf',sep = ''),
       width=85,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")




KEGGenrich %>%
  filter(Measure %in% c("G","F") & Contrast=="T-C" ) %>%
  write_csv(paste0(odir,"/KEGG_pathway_enrichment.csv"))






#Chebi enrichment

# source("https://bioconductor.org/biocLite.R")
# biocLite('ontoCAT')
# 
# install.packages('ontologyIndex')
# install.packages('ontoCAT')
# 
# library(ontologyIndex)
# 
# library(ontoCAT)
# efo <- getEFO()
# 
# chebi<-getOntology('../ChEbi/chebi.obo')



pmcs.c<-subset(pmcs,ChEBI_ID!='')[,c('Plate','Well','Instance','ChEBI_ID')]
pmcs.c$fChEBI_ID<-paste0('CHEBI:',pmcs.c$ChEBI_ID)

head(pmcs.c)

pmcs.counts<-ddply(pmcs.c,.(Plate,Well),summarise,Count=length(ChEBI_ID))

pmcs.co<-merge(pmcs.c,pmcs.counts,by=c('Plate','Well'))

head(pmcs.co)


head(encast)

chedata<-merge(pmcs.co,encast,by=c('Plate','Well'))


head(chedata)


chedata$Ecoli_All<-ifelse(chedata$Ecoli_FDR<0.05,1/chedata$Count,0)
chedata$Ecoli_Up<-ifelse(chedata$Ecoli_FDR<0.05 & chedata$Ecoli_logFC>0, 1/chedata$Count, 0)
chedata$Ecoli_Down<-ifelse(chedata$Ecoli_FDR<0.05 & chedata$Ecoli_logFC<0, 1/chedata$Count, 0)


chedata$Celegans_All<-ifelse(chedata$Celegans_FDR<0.05,1/chedata$Count,0)
chedata$Celegans_Up<-ifelse(chedata$Celegans_FDR<0.05 & chedata$Celegans_logFC>0, 1/chedata$Count, 0)
chedata$Celegans_Down<-ifelse(chedata$Celegans_FDR<0.05 & chedata$Celegans_logFC<0, 1/chedata$Count, 0)



head(chedata)

View(chedata)

write.csv(chedata,paste0(odir,'/ChEBI_annotation.csv'))




chebi<-data.frame()

chfls<-list.files('Summary/ChEBI_enrichment/')

for (fl in chfls) {
  print(fl)
  ed<-gsub('.txt','',fl,fixed = TRUE)
  ed<-gsub('ChEBI_','',ed)
  ce<-read.table(paste0('Summary/ChEBI_enrichment/',fl),sep='\t',header = TRUE)
  ce$Contrast<-strsplit(ed,'_')[[1]][1]
  ce$Type<-strsplit(ed,'_')[[1]][2]
  chebi<-rbind(chebi,ce)
}


View(chebi)

chebi.cast<-dcast(chebi,ChEBI_Name~Contrast+Type,value.var='Corr.PValue')








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


lm_eqn = function(m) {
  fres<-summary(m)
  l <- list(a = format(abs(coef(m)[[2]]), digits = 2),
            b = format(coef(m)[[1]], digits = 2),
            r2 = format(summary(m)$r.squared, digits = 2),
            p2 = format(pf(fres$fstatistic[[1]], fres$fstatistic[[2]], fres$fstatistic[[3]],lower.tail = FALSE)[[1]], digits = 2,scientific=TRUE));
  
  if (coef(m)[[2]] >= 0)  {
    cof <- substitute(italic(y) == b + a %.% italic(x),l)
    full <- substitute(italic(y) == b + a %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  } else {
    cof <- substitute(italic(y) == b - a %.% italic(x),l) 
    full <- substitute(italic(y) == b - a %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  }
  
  stat<-substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  return(list('Coef'=as.character(as.expression(cof)),
              'Stat'=as.character(as.expression(stat)),
              'Full'=as.character(as.expression(full)),
              'Atop'=as.character(as.expression(paste(cof,'\n',stat,sep='') ))))                 
}





info<-read_csv('~/Scripts/Growth_analysis/Biolog_metabolites.csv') %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A'))


PM1PM2<-subset(info,Plate %in% c('PM1','PM2A'))


#Pick data to show

logFDRbreaks<-c(-1,1.3,2,3,14)
logFDRbins<-c('N.S.','p<0.05','p<0.01','p<0.001')
logFDRbinsl<-c('N.S.','<0.05','<0.01','<0.001')

logFDRcols<-c("gray40", "red4", "red3",'red')


#Carboxylic acids to remove
coxy<-c('Caproic Acid','Capric Acid','4-Hydroxy Benzoic Acid','2-Hydroxy Benzoic Acid')


selmets<-info %>%
  filter(!(Plate %in% c('PM3B','PM4A') &
             Metabolite %in% PM1PM2$Metabolite ) &
           !Metabolite %in% c('Negative Control',coxy))


results.eco<-read_csv('Summary/Ecoli_results.csv')  %>%
  select(Description:Contrast_type,Measure,MetaboliteU,Index,logFC:NE,-c(pStars,FDRStars,t.value,p.value))


results.eco %>%
  filter(Index=='PM2A-E11')

metorder<-results.eco %>%
  filter(MetaboliteU %in% selmets$MetaboliteU & Contrast=="T-C" & Measure=="G") %>%
  arrange(logFC) %>%
  mutate(MetaboliteU=factor(MetaboliteU,levels=MetaboliteU)) %>%
  pull(MetaboliteU) %>%
  as.character


results.cel<-read_csv('Celegans/Summary/Celegans_results_withNGM.csv') %>%
  filter(! Index %in% c('Controls-A1','Controls-B1')) %>%
  select(Index,logFC:logFDR_bin) %>%
  mutate(Organism='C. elegans',
         Measure="F",
         Description="C. elegans rescue via GFP::Pacs-2 fluorescence reduction",
         Contrast_type="Treatment",
         Contrast="T-C",
         logFCrev=-logFC,
         PErev=-NE,
         NErev=-PE) %>%
  select(Organism,Measure,Contrast,Contrast_type,Description,Index,logFC=logFCrev,SE,PE=PErev,NE=NErev,FDR) # Reverse logFC for rescue



results.ecocel<-results.eco %>%
  select(-MetaboliteU) %>%
  mutate(Organism='E. coli') %>%
  rbind(results.cel) %>%
  left_join(info) %>%
  select(Organism,Measure,Contrast,Contrast_type,Description,Plate,Well,Index,MetaboliteU,EcoCycID,KEGG_ID,logFC,SE,PE,NE,FDR,everything()) %>%
  filter(MetaboliteU %in% selmets$MetaboliteU) %>%
  mutate(logFDR=-log10(FDR),
         MetaboliteU=factor(MetaboliteU,levels=metorder,labels=metorder),
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


kegg.mappings<-read_csv(paste0(odir,"/All_KEGG_mappings_Complete.csv")) %>%
  separate_rows(all.mapped) %>%
  filter(!is.na(all.mapped)) %>%
  group_by(all.mapped) %>%
  summarise(Pathway_IDs=paste(unique(PathwayID),collapse = ';'),
            Pathways=paste(unique(Description),collapse = ';'))



all.results.KEGG<-results.ecocel %>%
  left_join(kegg.mappings,by = c('KEGG_ID'='all.mapped')) %>%
  left_join(EC_mappings) %>%
  select(Organism,Organism_short,OMC:KEGG_ID,Pathway_IDs,Pathways,EcoCyc_Instances,EcoCyc_Classes,everything())

all.results.KEGG %>%
  write_csv(paste0(odir,"/All_results_withKEGGEcoCycClass.csv")) 




all.results.sbs<-results.ecocel %>%
  gather(Stat,Value,logFC,SE,FDR) %>%
  filter(Measure %in% c('G',"F") & Contrast %in% c("T-C","T","C") ) %>%
  mutate(OMC2=paste(Organism_short,Contrast,Stat,sep='_')) %>%
  select(OMC2,Plate:KEGG_ID,Value) %>%
  spread(OMC2,Value) %>%
  left_join(kegg.mappings,by = c('KEGG_ID'='all.mapped')) %>%
  left_join(EC_mappings) %>%
  select(Plate:KEGG_ID,Pathway_IDs,Pathways,EcoCyc_Instances,EcoCyc_Classes,`Ec_T-C_logFC`,`Ec_T-C_SE`,`Ec_T-C_FDR`,`Ce_T-C_logFC`,`Ce_T-C_SE`,`Ce_T-C_FDR`,Ec_C_logFC,Ec_C_SE,Ec_C_FDR,Ec_T_logFC,Ec_T_SE,Ec_T_FDR)
  

all.results.sbs %>%
  write_csv(paste0(odir,"/All_results_side_by_side_withKEGGEcoCycClass.csv")) 




#For pathway analysis
amp<-2
cbrks<-seq(-amp,amp,by=1)
gradcols<-c('blue4','blue','gray80','red','red4')

  
PTSmets<-c("N-Acetyl-D-Glucosamine",
           "D-Mannose-6-Phosphate",
           "D-Glucosamine-6-Phosphate")

all.results.KEGG %>%
  mutate(MetKEGG=paste(MetaboliteU,KEGG_ID)) %>%
  filter(Measure %in% c('G',"F") & Contrast=="T-C" ) %>%
  filter(str_detect(Pathway_IDs,'eco02060') | MetaboliteU %in% PTSmets ) %>% 
  ggplot(aes(x=Organism,y=MetKEGG,fill=logFC))+
  geom_tile()+
  scale_fill_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp))+
  theme(axis.text.x = element_text(angle=90))

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PTS_Celegans_Ecoli_logFC.pdf",sep=''),
             width=5,height=10)



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
  xlab("E. coli growth rescue/C. elegans phenotype rescue, logFC")+
  theme( axis.text.y=element_blank(),
         legend.position = "top")

ggsave(file=paste0(odir,"/EcoliCelegans_logFC_tall_Treatment_tiny.pdf"),
             width=55,height=80,units='mm',scale=2,device=cairo_pdf,family="Arial")

ggsave(file=paste0(odir,"/EcoliCelegans_logFC_tall_Treatment_tiny_7cm.pdf"),
       width=60,height=70,units='mm',scale=2,device=cairo_pdf,family="Arial")




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



showmets<-c('alpha-D-Glucose',
  'D-Ribose',
  'D-Arabinose',
  'Glycerol',
  'L-Serine',
  'Adenosine',
  'Acetoacetic Acid',
  'Itaconic Acid')

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

ggsave(file=paste0(odir,"/Scatter_Control-Treatment_Complete_Selmets_large.pdf"),
             width=110,height=110,units='mm',scale=2,device=cairo_pdf,family="Arial")






amp<-2.2
cbrks<-c(-2,-1,0,1,2)
gradcols<-c('blue4','blue','gray80','red','red4')

eccol<-'E. coli\ngrowth\nrescue,\nlogFC'
cecol<-'C. elegans\nphenotype\nrescue,\nlogFC'

thres<-0.05

gradcols<-c('blue4','blue','gray80','red','red4')


fit_TC<-lm(y_logFC~x_logFC,
           data=filter(ecocelmulti2,
                       x_OMC=="Ec_G_C",
                       y_OMC=="Ec_G_T"))
result_TC<-summary(fit_TC)


gsum_TC<-ecocelmulti2 %>%
  filter(x_OMC %in% c("Ec_G_C"),
         y_OMC=="Ec_G_T") %>%
  group_by(x_OMC) %>%
  do(broom::tidy(lm(y_logFC~x_logFC,data=.)))

rsum_TC<-ecocelmulti2 %>%
  filter(x_OMC %in% c("Ec_G_C"),
         y_OMC=="Ec_G_T") %>%
  group_by(x_OMC) %>%
  do(broom::glance(lm(y_logFC~x_logFC,data=.)))



TC_r<-result_TC$r.squared
TC_i<-gsum_TC[gsum_TC$term=='(Intercept)',]$estimate
TC_s<-gsum_TC[gsum_TC$term=='x_logFC',]$estimate


ecocelmulti %>%
  filter(x_OMC=="Ec_G_C",
         y_OMC=="Ec_G_T",
         z_OMC=="Ec_G_T-C") %>%
  ggplot(aes(y=y_logFC,x=x_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(intercept=1.125,slope=1),color='green4',linetype='longdash')+
  geom_abline(aes(intercept=TC_i,slope=TC_s),color='red')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  #ggtitle('Scatterplot of metformin and metabolite supplementation effects')+
  scale_x_continuous(breaks=-10:10) +
  scale_y_continuous(breaks=-10:10) +
  geom_point(aes(color=z_logFC))+
  coord_cartesian(xlim=c(-2,3),ylim = c(-2,3.5))+
  xlab('Growth logFC vs NGM - Control')+
  ylab('Growth logFC vs NGM - +50mM Metformin')+
  labs(color=eccol)+
  geom_text_repel(aes(label=ifelse(MetaboliteU %in% showmets, as.character(MetaboliteU),"" )), color="black"  )+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=eccol)+
  annotate('text',x = 0.5, y =-2, label = lm_eqn(fit_TC)$Full, parse = TRUE,color ='red',size=4)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file=paste0(odir,"/Scatter_Control-Treatment_Complete_Selmets_clean_small.pdf"),
       width=55,height=42,units='mm',scale=2,device=cairo_pdf,family="Arial")



# Update simple cor
ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_C",
         y_OMC=="Ec_G_T") %>%
  ggplot(aes(x=x_logFC,y=y_logFC) )+
  #geom_rug(size=1,alpha=0.2,color="gray20")+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1, intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_abline(aes(intercept=1.125,slope=1),color='green4',linetype='longdash')+
  geom_abline(aes(slope=TC_s,
                  intercept=TC_i),
              color='red',size=0.5)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(color=ifelse(MetaboliteU %in% showmets,"red","black")),show.legend = FALSE)+
  scale_color_manual(values = c("black","red"))+
  xlab('Growth logFC vs NGM - Control')+
  ylab('Growth logFC vs NGM - +50mM Metformin')+
  scale_x_continuous(breaks=seq(-5,5,by=1))+
  scale_y_continuous(breaks=seq(-5,5,by=1))+
  geom_text_repel(aes(label=ifelse(MetaboliteU %in% showmets, as.character(MetaboliteU),"" ),color=ifelse(MetaboliteU %in% showmets,"red","black")), show.legend = FALSE)+
  annotate('text',x = 0, y =-2, label = lm_eqn(fit_TC)$Full, parse = TRUE,color ='red',size=4)+
  coord_cartesian(xlim=c(-2,3),ylim = c(-2,3.5))


ggsave(file=paste0(odir,"/Correlations_EcT_EcC.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



# Updated based on nutrient groups
metgroups <- c('Sign. in E. coli',
               'Sign. in C. elegans',
               'Concordant',
               'Discordant',
               'N.S.')
metcolours <- c('blue3','red3','purple3','orange2','gray70')


ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T-C", 
         y_OMC=="Ce_F_T-C") %>%
  mutate(Metgroup = case_when(x_FDR<0.05 & y_FDR<0.05 & sign(x_logFC) == sign(y_logFC) ~ 'Concordant',
                              x_FDR<0.05 & y_FDR<0.05 & sign(x_logFC) != sign(y_logFC) ~ 'Discordant',
                              x_FDR<0.05 & y_FDR>=0.05 ~ 'Sign. in E. coli',
                              x_FDR>=0.05 & y_FDR<0.05 ~ 'Sign. in C. elegans',
                              TRUE ~ 'N.S.'),
         Metgroup = factor(Metgroup, levels = metgroups)) %>%
  ggplot(aes(x=x_logFC,y=y_logFC,colour=Metgroup))+
  geom_vline(xintercept = 0,alpha=0.5)+#color='gray70',
  geom_hline(yintercept = 0,alpha=0.5)+#color='gray70',
  geom_abline(aes(slope=1,intercept=0),color='gray70',alpha=0.5,linetype='longdash')+
  geom_vline(xintercept = 1.125,color='green4', linetype='longdash')+
  geom_abline(aes(slope=TCs,
                  intercept=TCi),
              color='red',size=0.5)+
  geom_errorbar(aes(ymin=y_PE,ymax=y_NE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  #geom_point(size=3,alpha=0.6,shape=20)+
  geom_point(size=3,shape=20)+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  xlab('Growth rescue in E. coli, logFC')+
  ylab('Phenotype rescue in C. elegans, logFC')+
  geom_text_repel(aes(label=ifelse(MetaboliteU %in% metselection, as.character(MetaboliteU),'')),
                  size=lblsize,nudge_y = 0.3,force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha,
                  show.legend=FALSE)+
  scale_color_manual(values = metcolours)+
  # scale_colour_gradientn(colours = gradcols,
  #                        breaks=cbrks,limits=c(-2,2))+
  labs(color='Nutrient effect')+
  coord_cartesian(xlim=c(-2.5,2.5),ylim = c(-1.5,2))+
  #ggtitle(paste('C. elegans and E. coli phenotype rescue in metformin treatment by metabolites: ',grp,sep='') )+
  annotate('text',x = 0, y =-1.5, label = lm_eqn(fitTC)$Full, parse = TRUE,color ='red',size=4)+
  guides(colour = guide_legend(nrow=2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top")


ggsave(file=paste0(odir,"/Celegans_Ecoli_effects_solid.pdf"),
       width=74,height=70,units='mm',scale=2,device=cairo_pdf,family="Arial")



metgroups2 <- c('Concordant',
               'Discordant',
               'Significant',
               'N.S.')
metcolours2 <- c('red3','orange2','blue3','gray70')

ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T-C", 
         y_OMC=="Ce_F_T-C") %>%
  mutate(Metgroup = case_when(x_FDR<0.05 & y_FDR<0.05 & sign(x_logFC) == sign(y_logFC) ~ 'Concordant',
                              x_FDR<0.05 & y_FDR<0.05 & sign(x_logFC) != sign(y_logFC) ~ 'Discordant',
                              x_FDR<0.05 | y_FDR<0.05 ~ 'Significant',
                              TRUE ~ 'N.S.'),
         Metgroup = factor(Metgroup, levels = metgroups2)) %>%
  ggplot(aes(x=x_logFC,y=y_logFC,colour=Metgroup))+
  geom_vline(xintercept = 0,alpha=0.5)+#color='gray70',
  geom_hline(yintercept = 0,alpha=0.5)+#color='gray70',
  geom_abline(aes(slope=1,intercept=0),color='gray70',alpha=0.5,linetype='longdash')+
  geom_vline(xintercept = 1.125,color='green4', linetype='longdash')+
  geom_abline(aes(slope=TCs,
                  intercept=TCi),
              color='red',size=0.5)+
  geom_errorbar(aes(ymin=y_PE,ymax=y_NE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  #geom_point(size=3,alpha=0.6,shape=20)+
  geom_point(size=3,shape=20)+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  xlab('Growth rescue in E. coli, logFC')+
  ylab('Phenotype rescue in C. elegans, logFC')+
  geom_text_repel(aes(label=ifelse(MetaboliteU %in% metselection, as.character(MetaboliteU),'')),
                  size=lblsize,nudge_y = 0.3,force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha,
                  show.legend=FALSE)+
  scale_color_manual(values = metcolours2)+
  # scale_colour_gradientn(colours = gradcols,
  #                        breaks=cbrks,limits=c(-2,2))+
  labs(color='Nutrient effect')+
  coord_cartesian(xlim=c(-2.5,2.5),ylim = c(-1.5,2))+
  #ggtitle(paste('C. elegans and E. coli phenotype rescue in metformin treatment by metabolites: ',grp,sep='') )+
  annotate('text',x = 0, y =-1.5, label = lm_eqn(fitTC)$Full, parse = TRUE,color ='red',size=4)+
  guides(colour = guide_legend(nrow=2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top")


ggsave(file=paste0(odir,"/Celegans_Ecoli_effects_only_conc_disc_solid.pdf"),
       width=74,height=70,units='mm',scale=2,device=cairo_pdf,family="Arial")


metgroups3 <- c('Concordant rescue',
                'Concordant aggravate',
                'Discordant',
                'Significant',
                'N.S.')
metcolours3 <- c('red3','purple2','orange3','blue3','gray70')

ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T-C", 
         y_OMC=="Ce_F_T-C") %>%
  mutate(Metgroup = case_when(x_FDR<0.05 & y_FDR<0.05 & x_logFC>0 & y_logFC>0 ~ 'Concordant rescue',
                              x_FDR<0.05 & y_FDR<0.05 & x_logFC<0 & y_logFC<0 ~ 'Concordant aggravate',
                              x_FDR<0.05 & y_FDR<0.05 & sign(x_logFC) != sign(y_logFC) ~ 'Discordant',
                              x_FDR<0.05 | y_FDR<0.05 ~ 'Significant',
                              TRUE ~ 'N.S.'),
         Metgroup = factor(Metgroup, levels = metgroups3)) %>%
  ggplot(aes(x=x_logFC,y=y_logFC,colour=Metgroup))+
  geom_vline(xintercept = 0,alpha=0.5)+#color='gray70',
  geom_hline(yintercept = 0,alpha=0.5)+#color='gray70',
  geom_abline(aes(slope=1,intercept=0),color='gray70',alpha=0.5,linetype='longdash')+
  geom_vline(xintercept = 1.125,color='green4', linetype='longdash')+
  geom_abline(aes(slope=TCs,
                  intercept=TCi),
              color='red',size=0.5)+
  geom_errorbar(aes(ymin=y_PE,ymax=y_NE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  #geom_point(size=3,alpha=0.6,shape=20)+
  geom_point(size=3,shape=20)+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  xlab('Growth rescue in E. coli, logFC')+
  ylab('Phenotype rescue in C. elegans, logFC')+
  geom_text_repel(aes(label=ifelse(MetaboliteU %in% metselection, as.character(MetaboliteU),'')),
                  size=lblsize,nudge_y = 0.3,force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha,
                  show.legend=FALSE)+
  scale_color_manual(values = metcolours3)+
  # scale_colour_gradientn(colours = gradcols,
  #                        breaks=cbrks,limits=c(-2,2))+
  labs(color='Nutrient effect')+
  coord_cartesian(xlim=c(-2.5,2.5),ylim = c(-1.5,2))+
  #ggtitle(paste('C. elegans and E. coli phenotype rescue in metformin treatment by metabolites: ',grp,sep='') )+
  annotate('text',x = 0, y =-1.5, label = lm_eqn(fitTC)$Full, parse = TRUE,color ='red',size=4)+
  guides(colour = guide_legend(nrow=2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top")


ggsave(file=paste0(odir,"/Celegans_Ecoli_effects_only_conc_disc_with_direction_solid.pdf"),
       width=74,height=70,units='mm',scale=2,device=cairo_pdf,family="Arial")







ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T-C", 
         y_OMC=="Ce_F_T-C") %>%
  filter(x_logFC< -2) %>%
  View


# Update Gradient
ecocelmulti %>%
  filter(x_OMC=="Ec_G_T-C",
         y_OMC=="Ce_F_T-C",
         z_OMC=="Ce_F_T-C") %>%
  ggplot(aes(y=y_logFC,x=x_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(intercept=TCi,slope=TCs),color='red')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_vline(xintercept = 1.125,color='green4', linetype='longdash')+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  #ggtitle('Scatterplot of metformin and metabolite supplementation effects')+
  scale_x_continuous(breaks=-10:10) +
  scale_y_continuous(breaks=-10:10) +
  geom_point(aes(color=z_logFC))+
  coord_cartesian(xlim=c(-2.5,2.5),ylim = c(-1.5,2))+
  xlab('Growth rescue in E. coli, logFC')+
  ylab('Phenotype rescue in C. elegans, logFC')+
  labs(color=cecol)+
  geom_text_repel(aes(label=ifelse(MetaboliteU %in% showmets, as.character(MetaboliteU),"" ),color=z_logFC), show.legend = FALSE )+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=cecol)+
  annotate('text',x = 0, y =-1.5, label = lm_eqn(fitTC)$Full, parse = TRUE,color ='red',size=4)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file=paste0(odir,"/Celegans_Ecoli_Ce_gradient.pdf"),
       width=84,height=70,units='mm',scale=2,device=cairo_pdf,family="Arial")





#C elegans, E coli comparison - pole 

xmin<- -3
xmax<-3
ymin<- -2
ymax<- 2

errcolor<-"grey90"
sbrks<-seq(0,3,by=1)

cbrks<-seq(-1,1,by=1)
gradcols<-c('gray80','blue3','gray80','red3','gray80')



metselection<-c('alpha-D-Glucose',
           'D-Ribose',
           'D-Arabinose',
           'L-Serine',
           'Adenosine',
           'Acetoacetic Acid',
           'Itaconic Acid',
           'Glycerol')




compdata<-ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T-C", 
         y_OMC=="Ce_F_T-C")


fitTC<-lm(y_logFC~x_logFC,data=compdata)
resultTC<-summary(fitTC)



gsum<-ecocelmulti2 %>%
  filter(x_OMC %in% c("Ec_G_C","Ec_G_T",'Ec_G_T-C'),
         y_OMC=="Ce_F_T-C") %>%
  group_by(x_OMC) %>%
  do(broom::tidy(lm(y_logFC~x_logFC,data=.)))

rsum<-ecocelmulti2 %>%
  filter(x_OMC %in% c("Ec_G_C","Ec_G_T",'Ec_G_T-C'),
         y_OMC=="Ce_F_T-C") %>%
  group_by(x_OMC) %>%
  do(broom::glance(lm(y_logFC~x_logFC,data=.)))



TCr<-resultTC$r.squared
TCi<-gsum[gsum$x_OMC=='Ec_G_T-C' & gsum$term=='(Intercept)',]$estimate
TCs<-gsum[gsum$x_OMC=='Ec_G_T-C' & gsum$term=='x_logFC',]$estimate




grp<-'Complete'
ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T-C", 
         y_OMC=="Ce_F_T-C") %>%
  ggplot(aes(x=x_logFC,y=y_logFC,color=pole(-y_logFC,x_logFC)))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5)+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5)+
  geom_vline(xintercept = 1.125,color='green4',linetype='longdash')+
  geom_abline(aes(slope=TCs,
                  intercept=TCi),
              color='red',size=0.5)+
  geom_errorbar(aes(ymin=y_PE,ymax=y_NE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=sqrt( y_logFC^2 + x_logFC^2)))+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  xlab('Growth rescue in E. coli, logFC')+
  ylab('Phenotype rescue in C. elegans, logFC')+
  geom_text_repel(aes(label=ifelse(MetaboliteU %in% metselection, as.character(MetaboliteU),'')),
                  size=lblsize,nudge_y = 0.3,force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-2,2))+
  scale_size(range = c(0.25, 7),breaks=sbrks,name='Strength\nof effect')+
  labs(color='Consistency\nof effect')+
  coord_cartesian(xlim=c(-2.5,2.5),ylim = c(-1.5,2))+
  #ggtitle(paste('C. elegans and E. coli phenotype rescue in metformin treatment by metabolites: ',grp,sep='') )+
  #annotate('text',x = -1, y =-1.5, label = lm_eqn(fitTC)$Full, parse = TRUE,color ='red',size=4)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



ggsave(file=paste0(odir,"/Celegans_Ecoli_effect_correlation_Complete_Selection.pdf"),
       width=110,height=82,units='mm',scale=2,device=cairo_pdf,family="Arial")






#Distribution

ecocelmulti2 %>%
  filter(x_OMC=='Ec_G_T-C' & y_OMC=='Ce_F_T-C' & y_logFC>0 & y_FDR<0.05) %>%
  ggplot(aes(x=x_logFC))+
  geom_histogram(binwidth = 0.2)+
  xlab('E. coli growth rescue, logFC')+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  coord_cartesian(xlim=c(-2,2))

ggsave(file=paste0(odir,"/Ecoli_growth_T-C_distribution_Ce-filter.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")


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



gcomp<-ecocelmulti %>%
  filter(x_OMC=="Ec_G_C",
         y_OMC=="Ec_G_T",
         z_OMC=="Ce_F_T-C")




fitC<-lm(z_logFC~x_logFC,data=gcomp)
resultC<-summary(fitC)
resultC
resultC$r.squared


fitT<-lm(z_logFC~y_logFC,data=gcomp)
resultT<-summary(fitT)
resultT
resultT$r.squared


gsum<-ecocelmulti2 %>%
  filter(x_OMC %in% c("Ec_G_C","Ec_G_T",'Ec_G_T-C'),
         y_OMC=="Ce_F_T-C") %>%
  group_by(x_OMC) %>%
  do(broom::tidy(lm(y_logFC~x_logFC,data=.)))

rsum<-ecocelmulti2 %>%
  filter(x_OMC %in% c("Ec_G_C","Ec_G_T",'Ec_G_T-C'),
         y_OMC=="Ce_F_T-C") %>%
  group_by(x_OMC) %>%
  do(broom::glance(lm(y_logFC~x_logFC,data=.)))



gsum
rsum

Cr<-resultC$r.squared
Ci<-gsum[gsum$x_OMC=='Ec_G_C' & gsum$term=='(Intercept)',]$estimate
Cs<-gsum[gsum$x_OMC=='Ec_G_C' & gsum$term=='x_logFC',]$estimate

Ti<-gsum[gsum$x_OMC=='Ec_G_T' & gsum$term=='(Intercept)',]$estimate
Ts<-gsum[gsum$x_OMC=='Ec_G_T' & gsum$term=='x_logFC',]$estimate
Tr<-resultT$r.squared




ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_C",
         y_OMC=="Ce_F_T-C") %>%
  ggplot(aes(x=x_logFC,y=y_logFC) )+
  geom_rug(size=1,alpha=0.2,color="gray20")+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=Cs,
                  intercept=Ci),
              color='red',size=0.5)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+
  ylab('C. elegans phenotype rescue, logFC')+
  xlab('E. coli growth vs NGM, logFC')+
  scale_x_continuous(breaks=seq(-5,5,by=1))+
  scale_y_continuous(breaks=seq(-5,5,by=1))+
  #annotate('text',x = -1, y =-2, label = lm_eqn(fitC)$Full, parse = TRUE,color ='red',size=4)+
  coord_cartesian(xlim=c(-3,4),ylim = c(-2,2))


ggsave(file=paste0(odir,"/Correlations_CeF_EcC_rug.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")





ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T",
         y_OMC=="Ce_F_T-C") %>%
  ggplot(aes(x=x_logFC,y=y_logFC) )+
  geom_rug(size=1,alpha=0.2,color="gray20")+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=Ts,
                  intercept=Ti),
              color='red',size=0.5)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+

  ylab('C. elegans phenotype rescue, logFC')+
  xlab('E. coli growth vs NGM +50 mM Metformin, logFC')+
  scale_x_continuous(breaks=seq(-5,5,by=1))+
  scale_y_continuous(breaks=seq(-5,5,by=1))+
  #annotate('text',x = -1, y =-2, label = lm_eqn(fitT)$Full, parse = TRUE,color ='red',size=4)+
  coord_cartesian(xlim=c(-3,4),ylim = c(-2,2))

ggsave(file=paste0(odir,"/Correlations_CeF_EcT_rug.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



ecocelmulti2 %>%
  filter(x_OMC=="Ec_G_T-C",
         y_OMC=="Ce_F_T-C") %>%
  ggplot(aes(x=x_logFC,y=y_logFC) )+
  #geom_rug(size=1,alpha=0.2,color="gray20")+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_vline(xintercept = 1.125,color='blue',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=TCs,
                  intercept=TCi),
              color='red',size=0.5)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+
  ylab('C. elegans phenotype rescue, logFC')+
  xlab('E. coli growth rescue, logFC')+
  scale_x_continuous(breaks=seq(-5,5,by=1))+
  scale_y_continuous(breaks=seq(-5,5,by=1))+
  #annotate('text',x = -1, y =-2, label = lm_eqn(fitTC)$Full, parse = TRUE,color ='red',size=4)+
  coord_cartesian(xlim=c(-2.5,2.5),ylim = c(-1.5,2))

ggsave(file=paste0(odir,"/Correlations_CeF_EcT-C.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")





#For summaries
#Time series
ectimem<-read_csv('Data/Data.csv') %>%
  filter(Data=='750nm_f') %>%
  gather(Time_sec,OD,ends_with(".0")) %>%
  mutate(Time_sec=as.numeric(as.character(Time_sec)),
         Time_min=Time_sec/60,
         Time_h=Time_min/60,
         Metformin_mM=as.factor(Metformin_mM))

tsum<-ectimem %>%
  group_by(Index,Type,Metformin_mM,Time_sec,Time_h) %>%
  summarise(Mean=mean(OD),SD=sd(OD),SE=SD/sqrt(length(OD))) %>%
  left_join(info) %>%
  mutate(Organism='E. coli')



results.in<-read_csv('Celegans/Summary/Celegans_results_withNGM.csv') %>%
  select(-c(Group,Metabolite_Uniq,Metabolite,EcoCycID,KEGG_ID,CAS_ID,Description)) %>%
  left_join(info)
  
data.Isum<-readRDS('Celegans/Celegans_fluorescence_distributions.rds') %>%
  filter(Measure=='Absolute') %>%
  select(-c(Group,Metabolite_Uniq,Metabolite)) %>%
  mutate(Metformin_mM=ifelse(Type=='Control',0,50) %>% as.factor) %>%
  left_join(info)


mlevels<-c('Negative Control_C',
              'alpha-D-Glucose',
              'D-Ribose',
           'D-Arabinose',
           'Glycerol',
              'L-Serine',
           'Adenosine',
              'Acetoacetic Acid',
           'Itaconic Acid')

mlabels<-c('NGM',
           'D-Glucose',
           'D-Ribose',
           'D-Arabinose',
           'Glycerol',
           'L-Serine',
           'Adenosine',
           'Acetoacetate',
           'Itaconate')


indxsel<-info %>%
  filter(MetaboliteU %in% mlevels) %>%
  pull(Index)

indxsel


info %>%
  filter(MetaboliteU %in% mlevels) %>%
  mutate(Col=str_remove_all(Well,'[:alpha:]') %>% as.numeric,
         Row=str_remove_all(Well,'[:digit:]'),
         Row=factor(Row,levels=LETTERS[1:8]),
         N=Col+(as.numeric(Row)-1)*12,
         IIndex=paste(Plate,formatC(N, width=3, flag="0"), sep="_")) %>%
  pull(IIndex)
           



indxsel<-c(indxsel,'Controls-A1')


#Distributions

Metcols <- c("#FF0000","#32006F")
names(Metcols) <- c(0,50)
Metlab<-'Metformin, mM'


tsum %>%
  filter(Index %in% indxsel) %>%
  mutate( MetaboliteU=factor(MetaboliteU,levels=mlevels,labels=mlabels)) %>%
  ggplot(aes(x=Time_h,y=Mean,fill=Metformin_mM,color=Metformin_mM))+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),color=NA,alpha=0.2)+
  scale_colour_manual(name = Metlab,values =Metcols)+
  scale_fill_manual(name = Metlab,values =Metcols)+
  geom_line()+
  scale_x_continuous(breaks=seq(0,24,by=12))+
  ylab("OD")+
  xlab("Time, h")+
  labs(fill="Type")+
  facet_grid(Organism~MetaboliteU)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(face = 'bold.italic',size = 10))


#New growth curves

ggsave(file=paste0(odir,"/Summary_Ecoli_growth_curves_new.pdf"),
       width=110,height=26,units='mm',scale=2,device=cairo_pdf,family="Arial")



data.Isel<-data.Isum %>%
  filter(Index %in% indxsel ) %>%
  mutate(MetaboliteU=ifelse(Index=='Controls-A1','Negative Control_C',MetaboliteU),
         MetaboliteU=factor(MetaboliteU,levels=mlevels,labels=mlabels),
         Metformin_mM=as.factor(Metformin_mM)) 




res.sel<-results.in %>%
  filter(Index %in% indxsel) %>% 
  mutate(Type=as.factor(ifelse(Index=='Controls-A1','Control','Treatment')),
         MetaboliteU=as.character(MetaboliteU),
         MetaboliteU=ifelse(Index=='Controls-A1','Negative Control_C',MetaboliteU),
         MetaboliteU=factor(MetaboliteU,levels=mlevels,label=mlabels),
         Metformin_mM=factor(Type,levels=c("Control","Treatment"),labels=c("0","50")),
         Organism='C. elegans')



res.sel %>%
  ggplot(aes(color=Metformin_mM))+
  geom_vline(aes(xintercept=res.sel[res.sel$Index=='PM1-A1',]$Raw_Q90_Mean),color='blue4',linetype='longdash',alpha=0.5)+
  geom_vline(aes(xintercept=res.sel[res.sel$Index=='Controls-A1',]$Raw_Q90_Mean),color='red4',linetype='longdash',alpha=0.5)+
  geom_rect(aes(xmin=Raw_Q90_Mean-Raw_Q90_SE,xmax=Raw_Q90_Mean+Raw_Q90_SE,fill=Metformin_mM),ymin=-Inf,ymax=Inf,alpha=0.2,color=NA)+
  geom_vline(aes(xintercept=Raw_Q90_Mean,color=Metformin_mM))+
  geom_ribbon(data=data.Isel,aes(x=LogBrightness_num,ymin=(Mean-SD)*100,ymax=(Mean+SD)*100,fill=Metformin_mM),alpha=0.5,color=NA)+
  geom_line(data=data.Isel,aes(x=LogBrightness_num,y=Mean*100,color=Metformin_mM))+
  scale_colour_manual(name = Metlab,values =Metcols)+
  scale_fill_manual(name = Metlab,values =Metcols)+
  ylab('Density, %')+
  xlab('log2 Fluorescence brightness, A.U.')+
  scale_x_continuous(breaks=seq(-20,20,by=2))+
  coord_cartesian(ylim=c(0,15))+
  facet_grid(Organism~MetaboliteU)+
  theme(legend.position = "none",
        #strip.background = element_blank(),
        strip.text.x =  element_blank(),
        strip.text.y = element_text(face = 'bold.italic',size = 10)
        )

ggsave(file=paste(odir,"/Summary_Celegans_brightness_distribution_new.pdf",sep=''),
             width=110,height=22.7,units='mm',scale=2,device=cairo_pdf,family="Arial")




#Comparisons summary


res.Ecsel<-read_csv(paste0(odir,"/Ecoli_unfiltered_results.csv")) %>%
  filter(Contrast %in% c("C","T")  & Measure=="G")  %>%
  # filter(Metabolite=="Negative Control") %>% 
  group_by(Group,Plate) %>%
  mutate(logFC=ifelse(Metabolite=="Negative Control",logFC-logFC[Metabolite=="Negative Control" & Contrast=='C'],logFC)) %>%
  group_by(Contrast,Group,Plate) %>%
  mutate(logFC=ifelse(Metabolite!="Negative Control",logFC+logFC[Metabolite=="Negative Control"],logFC),
         NE=logFC-SE,
         PE=logFC+SE ) %>%
  ungroup %>%
  filter(Index %in% indxsel) %>%
  mutate(MetaboliteU=factor(MetaboliteU,levels=mlevels,labels=mlabels),
         Metformin_mM=factor(Contrast,levels=c("C","T"),labels=c("0","50")),
         Organism="E. coli") %>%
  select(Plate,Well,Index,Contrast,Group,Metabolite,MetaboliteU,Metformin_mM,logFC,PE,NE,SE,Organism) 
  


res.Ecsel %>%
  ggplot(aes(x=MetaboliteU,y=logFC,color=Metformin_mM))+
  geom_hline(aes(yintercept=res.Ecsel[res.Ecsel$Index=='PM1-A1' & res.Ecsel$Contrast=='C',]$logFC),color='red4',linetype='longdash',alpha=0.5)+
  geom_hline(aes(yintercept=res.Ecsel[res.Ecsel$Index=='PM1-A1' & res.Ecsel$Contrast=='T',]$logFC),color='blue4',linetype='longdash',alpha=0.5)+
  geom_errorbar(aes(ymin=NE,ymax=PE),position = position_dodge(width = dwidth),width=0.25)+
  geom_point(position = position_dodge(width = dwidth),size=2)+
  scale_colour_manual(name = Metlab,values =Metcols)+
  scale_y_continuous(breaks=seq(-20,20,by=1))+
  coord_cartesian(ylim=c(-2,3))+
  labs(color="Metformin, mM")+
  ylab('Growth AUC vs NGM control, logFC')+
  xlab('Metabolite')+
  theme(axis.text.x = element_text(angle = 90,hjust=1))

dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Summary_Ecoli_growth_comparison_new.pdf",sep=''),
             width=4,height=4)


res.sel.adj<-res.sel %>%
  mutate(Raw_Q90_Mean=Raw_Q90_Mean-Raw_Q90_Mean[Index=='PM1-A1'])

dwidth<-0.5
res.sel.adj %>%
  ggplot(aes(x=MetaboliteU,y=Raw_Q90_Mean,color=Metformin_mM))+
  geom_hline(aes(yintercept=res.sel.adj[res.sel.adj$Index=='PM1-A1',]$Raw_Q90_Mean),color='blue4',linetype='longdash',alpha=0.5)+
  geom_hline(aes(yintercept=res.sel.adj[res.sel.adj$Index=='Controls-A1',]$Raw_Q90_Mean),color='red4',linetype='longdash',alpha=0.5)+
 
  geom_errorbar(aes(ymin=Raw_Q90_Mean-Raw_Q90_SE,ymax=Raw_Q90_Mean+Raw_Q90_SE),position = position_dodge(width = dwidth),width=0.25)+
  geom_point(position = position_dodge(width = dwidth),size=2)+
  scale_colour_manual(name = Metlab,values =Metcols)+
  scale_y_continuous(breaks=seq(-20,20,by=1))+
  coord_cartesian(ylim=c(-3,1))+
  ylab('Fluorescence brightness vs NGM treatment, logFC')+
  labs(color="Metformin, mM")+
  theme(axis.text.x = element_text(angle = 90,hjust=1))


dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Summary_Celegans_brightness_comparison_new.pdf",sep=''),
             width=4,height=4)



reflines<-data.frame(Organism=c("E. coli","E. coli","C. elegans","C. elegans"),
                     Metformin_mM=c(0,50,0,50),
                     Ref=c(0,-1.125,-2.261,0)) %>%
  mutate(Metformin_mM=as.factor(Metformin_mM))




#Combine two datasets
res.comb<-res.sel.adj %>%
  select(Index,Metformin_mM,MetaboliteU,logFC=Raw_Q90_Mean,SE=Raw_Q90_SE) %>%
  mutate(NE=logFC-SE,
         PE=logFC+SE,
         Organism="C. elegans") %>%
  rbind(res.Ecsel %>% select(Index,Organism,Metformin_mM,MetaboliteU,logFC,SE,NE,PE) ) %>%
  mutate(Organism=factor(Organism,levels=c("E. coli","C. elegans")))


res.comb %>%
  ggplot(aes(x=MetaboliteU,y=logFC,color=Metformin_mM))+
  geom_hline(data=reflines,aes(yintercept=Ref,color=Metformin_mM),linetype='longdash',alpha=0.5)+
  scale_y_continuous(breaks=seq(-20,20,by=1))+
  geom_errorbar(aes(ymin=NE,ymax=PE),position = position_dodge(width = dwidth),width=0.25)+
  scale_colour_manual(name = Metlab,values =Metcols)+
  geom_point(position = position_dodge(width = dwidth),size=1)+
  labs(color="Metformin, mM")+
  xlab('Metabolite')+
  #ylab('log2 AUC, OD*h')+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  facet_grid(Organism~.,scale="free_y")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(face = 'bold.italic'))


ggsave(file=paste0(odir,"/Summary_CelegansEcoli_comparison_new.pdf"),
             width=60,height=72,units='mm',scale=2,device=cairo_pdf,family="Arial")






#Venn diagrams
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




venn.expl<-results.ecocel %>%
  filter(OMC %in% c("Ec_G_T-C","Ce_F_T-C")) %>%
  mutate(Rescue=ifelse(logFC>0 & FDR<0.05,'TRUE',NA),
         Aggravate=ifelse(logFC<0 & FDR<0.05,'TRUE',NA) ) %>%
  gather(Type,Pass,Aggravate,Rescue) %>%
  unite(OT,Organism_short,Type) %>%
  select(Plate:KEGG_ID,OT,Pass) %>%
  spread(OT,Pass)

venn.expl %>%
  write_csv(paste0(odir,"/Venn_explanation.csv"),na="")




rescue<-venn.pass %>%
  filter(OMC %in% c("Ec_G_T-C","Ce_F_T-C") & Type!="All") %>%
  unite(OT,Organism_short,Type) %>%
  select(OT,List) %>%
  spread(OT,List) %>%
  as.list %>%
  unlist(recursive=FALSE)


vcols<-c("red","red4","blue","blue4")
grid::grid.draw(VennDiagram::venn.diagram(rescue, col=vcols,cat.col=vcols ,NULL))

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Celegans_Ecoli_updated_check.pdf",sep=''),
             width=3,height=3, useDingbats=FALSE)
dev.off()


#EcoCyc class enrichment

ecli.co<-read_csv('../Biolog/Ecoli_All_class_instances_non_redundant.csv') %>%
  select(-X1)

ecli.co %>%
  group_by(Class) %>%
  summarise()


#EcoCyc to Biolog links
pmcs<-read_csv('../Biolog/EColi_net/All_Media_Mix_clean_explicit.csv')




#Enrichment procedure
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




results.ecocel %>%
  left_join(pmcs %>% select(Plate,Well,Instance)) %>%
  left_join(ecli.co %>% filter(NoInstances>2)) %>%
  filter(!is.na(FDR) & !is.na(Class))



EC_mappings<-info %>%
  select(Plate,Well,Index) %>%
  left_join(pmcs %>% select(Plate,Well,Instance)) %>%
  left_join(ecli.co %>% filter(NoInstances>2)) %>%
  filter(!is.na(Class)) %>%
  group_by(Plate,Well,Index) %>%
  summarise(EcoCyc_Instances=paste(unique(Instance),collapse = ';'),
            EcoCyc_Classes=paste(unique(Class),collapse = ';'))



EC_mappings %>%
  write_csv(paste0(odir,'/EcoCyc_Class_mappings.csv'))




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


classorder<-c("Carbohydrates","Aldehydes-Or-Ketones",
              "Carboxylates","Amino-Sugars","Peptides","Alcohols",
              "All-Nucleosides","Ribonucleosides","Purines",
              "Alpha-Amino-Acids","Amino-Acids","L-Amino-Acids","Dicarboxylate",
              "Esters","Lactones","Glycosides")



EcoCyc.enrichment<-read_csv(paste0(odir,'/EcoCyc_enrichment.csv')) %>%
  mutate(OMC=ifelse(Organism=='E. coli',"Ec_G_T-C","Ce_F_T-C"))



EcoCyc.enrichment %>%
  filter(Test!="All" & OMC %in% c("Ec_G_T-C","Ce_F_T-C") ) %>% 
  group_by(Class) %>%
  filter( any(FDR<0.05)) %>%
  ungroup %>%
  clustorder('Class',c("Organism","Test"),'logFDR',descending=FALSE) %>% 
  mutate(Class=factor(Class,levels=rev(classorder) )) %>%
  PlotEnrichment("Test","Class",ncols = length(enrbrks))+
  facet_grid(~Organism)+
  ylab("EcoCyc metabolite class")


#Biolog class enrichment
#Classes preparation
classes<-info %>%
  select('Index','Plate','Well','MetaboliteU') %>%
  left_join(readxl::read_xlsx('../Biolog/Biolog_metabolites_EcoCyc.xlsx',sheet = 'Classes') %>% select('Index','Class1','Class2')) %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A')) %>%
  gather(Redundancy,Class,contains("Class")) %>%
  filter(Class!="") 




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
#Rocover possible links via classes
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


EcoCycenr %>%
  group_by(Organism,Measure,Contrast,Description,Type) %>%
  summarise(Count=n())


EcoCycenr %>%
  filter(Measure %in% c("G","F") & Contrast=="T-C" & p.value<0.05) %>%
  View


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

print(paste('Total KEGG pathways to plot:',length(allpaths)))

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



KEGGenrich %>%
  filter(Measure %in% c("G","F") & Contrast=="T-C" & Type !="All" ) %>%
  mutate(Interaction=ifelse(Type=='Up','Syn','Ant'),
           OMCT=paste(ifelse(Organism=='E. coli','Ec','Ce'),Interaction,sep="_")) %>%
  select(PathwayID,Pathway,OMCT,FDR) %>%
  spread(OMCT,FDR) %>%
  select(PathwayID,Pathway,Ec_Ant,Ec_Syn,Ce_Ant,Ce_Syn) %>%
  write_csv(paste0(odir,"/KEGG_pathway_enrichment_side-by-side.csv"))




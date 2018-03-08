library(tidyverse)
library(ggrepel)


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

theme_set(theme_Publication())
#theme_set(theme_light())



cwd<-"~/Dropbox/Projects/2015-Metformin/Biolog_Met_NGM/"
setwd(cwd)

#load('Biolog_combined.RData')
#save.image('Biolog_combined.RData')

odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



info<-read_csv('../Biolog/Biolog_metabolites_EcoCyc_Unique_PM1-PM5.csv') %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A'))


PM1PM2<-subset(info,Plate %in% c('PM1','PM2A'))

results.in<-read_csv('Celegans/Summary/Celegans_results_withNGM.csv')

results.castcomb<-read_csv('Summary/Combined_results_sidebyside_full.csv')

data.Isum<-readRDS('Celegans/Celegans_fluorescence_distributions.rds')





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

selectcast<-results.castcomb %>%
  filter(MetaboliteU %in% selmets$MetaboliteU) %>%
  arrange(`T-C_logFC`)%>%
  mutate(MetaboliteU=factor(MetaboliteU,levels=MetaboliteU),
         Order=row_number()) %>%
  mutate_at(c('T-C_logFDR_bin','Cel_logFDR_bin'),funs(factor(.,levels=logFDRbins,labels=logFDRbinsl)))


metorder<-as.character(selectcast$MetaboliteU)

results.cel<-read_csv('Celegans/Summary/Celegans_results_withNGM.csv')%>%
  filter(! Index %in% c('Controls-A1','Controls-B1')) %>%
  select(Index,logFC:logFDR_bin) %>%
  left_join(info) %>%
  mutate(Organism='C. elegans',
         logFCrev=-logFC,
         PErev=-NE,
         NErev=-PE) %>%
  select(Organism,Plate,Well,Index,MetaboliteU,logFC=logFCrev,SE,PE=PErev,NE=NErev,FDR,logFDR)

results.eco<-read_csv('Summary/Ecoli_results.csv') %>%
  filter(Measure=='G',Contrast=='T-C' )

results.ecocel<-results.eco %>%
  select(Plate,Well,Index,MetaboliteU,logFC,SE,PE,NE,FDR,logFDR) %>%
  mutate(Organism='E. coli') %>%
  rbind(results.cel) %>%
  filter(MetaboliteU %in% selmets$MetaboliteU) %>%
  mutate(MetaboliteU=factor(MetaboliteU,levels=metorder,labels=metorder),
         logFDR_bin=cut(logFDR, breaks=logFDRbreaks,labels=logFDRbinsl),
         Organism=factor(Organism,levels=c('E. coli','C. elegans') ) ) %>%
  group_by(Organism) %>%
  arrange(logFC) %>%
  mutate(Order=row_number()) %>%
  ungroup


View(selectcast)

#Just testing
selectcast %>%
  arrange(`T-C_logFC`) %>%
  select(MetaboliteU,Order,`T-C_logFC`,`T-C_FDR`,`T-C_logFDR`,`T-C_logFDR_bin`,Cel_logFC,Cel_logFDR_bin)


results.ecocel %>%
  #filter(Organism=='E. coli') %>%
  #arrange(Order) %>%
  #mutate(MetaboliteU=factor(MetaboliteU,levels=metorder,labels=metorder)) %>%
  ggplot(aes(y=MetaboliteU,x=logFC,color=logFDR_bin))+
  geom_vline(xintercept = 0,color='red',alpha=0.5)+
  geom_errorbarh(aes(xmin=NE,xmax=PE),height=0.5)+
  geom_point(size=1)+
  ylab('Metabolite')+
  xlab('logFC')+
  labs(color='FDR')+
  scale_colour_manual(values = logFDRcols)+
  scale_x_continuous(breaks=seq(-10,10) )+
  theme(axis.ticks.y=element_blank(),
        axis.line = element_line(colour = NA),
        axis.line.x = element_line(colour = NA),
        axis.line.y = element_line(colour = NA),
        #panel.border=element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())+
  facet_grid(~Organism)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/EcoliCelegans_logFC_tall_Treatment.pdf",sep=''),
             width=9,height=50)




results.ecocel %>%
  #filter(Organism=='E. coli') %>%
  #arrange(Order) %>%
  #mutate(MetaboliteU=factor(MetaboliteU,levels=metorder,labels=metorder)) %>%
  ggplot(aes(y=MetaboliteU,x=logFC,color=logFDR_bin))+
  geom_vline(xintercept = 0,color='red',alpha=0.5)+
  geom_errorbarh(aes(xmin=NE,xmax=PE),height=0,size=0.5)+
  geom_point(size=0.5)+
  ylab('Metabolites')+
  xlab('logFC')+
  labs(color='FDR')+
  scale_colour_manual(values = logFDRcols)+
  scale_x_continuous(breaks=seq(-10,10,by=2) )+
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = NA),
        axis.line.x = element_line(colour = NA),
        axis.line.y = element_line(colour = NA),
        #panel.border=element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())+
  facet_grid(~Organism)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/EcoliCelegans_logFC_tall_Treatment_tiny.pdf",sep=''),
             width=5,height=5)


selectcast  %>%
  ggplot(aes(y=MetaboliteU,x=`T-C_logFC`,color=`T-C_logFDR_bin`))+
    geom_vline(xintercept = 0,color='red',alpha=0.5)+
    geom_errorbarh(aes(xmin=`T-C_NE`,xmax=`T-C_PE`))+
    geom_point()+
    ylab('Metabolite')+
    xlab('Metabolite - metformin interaction as growth logFC vs NGM')+
    labs(color='FDR')+
    scale_colour_manual(values =logFDRcols)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ecoli_logFC_tall_Treatment.pdf",sep=''),
             width=9,height=50)



selectcast %>%
  arrange(Cel_logFC)%>%
  mutate(MetaboliteU=factor(MetaboliteU,levels=MetaboliteU)) %>%
  ggplot(aes(y=MetaboliteU,x=Cel_logFC,color=Cel_logFDR_bin))+
  geom_vline(xintercept = 0,color='red',alpha=0.5)+
  geom_errorbarh(aes(xmin=Cel_NE,xmax=Cel_PE))+
  geom_point()+
  xlab('Metabolite')+
  ylab('Metabolite - metformin interaction as growth logFC vs NGM')+
  labs(color='FDR')+
  scale_colour_manual(values = logFDRcols)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Celegans_logFC_tall_Fluorescence.pdf",sep=''),
             width=9,height=50)





head(selectcast)
#Scatter plots

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



#ncontrol<-subset(results.castcomb,Metabolite=='Negative Control' )

subset(selectcast,C_logFC<0 & C_FDR<0.05)


grp<-'Complete'

fit<-lm(T_logFC~C_logFC,selectcast)
lmeq<-lm_eqn(fit)
a<-fit$coefficients[[2]]
b<-fit$coefficients[[1]]

xmin<- -3
xmax<- 4
ymin<- -3
ymax<- 4

thres<-0.05

gradcols<-c('blue4','blue','gray80','red','red4')

ggplot(selectcast,aes(y=T_logFC,x=C_logFC,color=-Cel_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(intercept=1,slope=1),color='blue',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=T_NE,ymax=T_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=C_NE,xmax=C_PE),alpha=erralpha,color=errcolor,height=0)+
  ggtitle('Scatterplot of metformin and metabolite supplementation effects')+
  geom_point(size=3)+
  coord_cartesian(xlim=c(xmin,xmax),ylim = c(ymin,ymax))+
  xlab('Growth logFC vs NGM - Control')+
  ylab('Growth logFC vs NGM - +50mM Metformin')+
  geom_text_repel(aes(label=ifelse(abs(Cel_logFC)>0 & Cel_FDR<0.05, as.character(MetaboliteU),"" ))  )+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Scatter_Control-Treatment_Complete_Growth_Celegans_clean_refline_Cel<0.pdf",sep=''),
             width=10,height=6, useDingbats=FALSE)


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
selectcast %>%
  ggplot(aes(x=`T-C_logFC`,y=-Cel_logFC,color=pole(Cel_logFC,`T-C_logFC`)))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_errorbar(aes(ymin=-Cel_PE,ymax=-Cel_NE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=`T-C_NE`,xmax=`T-C_PE`),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=sqrt( Cel_logFC^2 + `T-C_logFC`^2)))+
  coord_cartesian(xlim=c(xmin,xmax),ylim = c(ymin,ymax))+
  scale_x_continuous(breaks=seq(xmin,xmax,by=1))+
  scale_y_continuous(breaks=seq(ymin,ymax,by=1))+
  xlab('Growth rescue in E. coli (metabolite effect normalised), logFC')+
  ylab('Phenotype rescue in C. elegans (acs-2 GFP), logFC')+
  geom_text_repel(aes(label=ifelse(Cel_FDR<0.05 &  `T-C_FDR`<0.05, as.character(MetaboliteU),'')),
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












indxsel<-c('PM1-A1','PM1-A2','PM4A-E6','PM1-G7') #'Controls-A1',

indxsel<-c('Controls-A1','PM1-A1','PM1-A2','PM1-G7','PM4A-E6')

mlevels<-c('Negative Control','L-Arabinose','Acetoacetic Acid','Phosphono Acetic Acid')
mlabels<-c('NGM','L-Arabinose','Acetoacetic Acid','Phosphono Acetic Acid')





tsum %>%
  filter(Index %in% indxsel) %>%
  mutate( Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels)) %>%
  ggplot(aes(x=Time_h,y=Mean,fill=Type,color=Type))+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),color=NA,alpha=0.2)+
  geom_line()+
  scale_x_continuous(breaks=seq(0,24,by=6))+
  ylab("OD")+
  xlab("Time, h")+
  labs(fill="Type")+
  facet_grid(~Metabolite)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Summary_Ecoli_growth_curves.pdf",sep=''),
             width=9,height=3)



data.Isel<-data.Isum %>%
  filter(Index %in% indxsel & Measure=='Absolute') %>%
  mutate(Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels))


res.sel<-results.in %>%
  filter(Index %in% indxsel) %>%
  mutate(Type=as.factor(ifelse(Index=='Controls-A1','Control','Treatment')),
         Metabolite=as.character(Metabolite),
         Metabolite=ifelse(Index=='Controls-A1','Negative Control',Metabolite),
         Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels))


ggplot(res.sel,aes(color=Type))+
  geom_vline(aes(xintercept=res.sel[res.sel$Index=='PM1-A1','Raw_Q90_Mean']),color='grey80',linetype='longdash')+
  geom_rect(aes(xmin=Raw_Q90_Mean-Raw_Q90_SE,xmax=Raw_Q90_Mean+Raw_Q90_SE,fill=Type),ymin=-Inf,ymax=Inf,alpha=0.2,color=NA)+
  geom_vline(aes(xintercept=Raw_Q90_Mean,color=Type))+
  geom_ribbon(data=data.Isel,aes(x=LogBrightness_num,ymin=(Mean-SD)*100,ymax=(Mean+SD)*100,fill=Type),alpha=0.5,color=NA)+
  geom_line(data=data.Isel,aes(x=LogBrightness_num,y=Mean*100,color=Type))+
  #geom_point()+
  ylab('Frequency, %')+
  xlab('log2 Brightness')+
  scale_x_continuous(breaks=seq(-20,20,by=1))+
  coord_cartesian(ylim=c(0,15))+
  facet_grid(.~Metabolite)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Summary_Celegans_brightness_distribution.pdf",sep=''),
             width=9,height=3)



dwidth<-0.5

ggplot(res.sel,aes(x=Metabolite,y=Raw_Q90_Mean,color=Type))+
  geom_hline(aes(yintercept=res.sel[res.sel$Index=='PM1-A1','Raw_Q90_Mean']),color='grey80',linetype='longdash')+
  scale_y_continuous(breaks=seq(-20,20,by=1),limits=c(-4,0))+
  geom_errorbar(aes(ymin=Raw_Q90_Mean-Raw_Q90_SE,ymax=Raw_Q90_Mean+Raw_Q90_SE),position = position_dodge(width = dwidth),width=0.25)+
  geom_point(position = position_dodge(width = dwidth),size=3)+
  ylab('log2 Fluoresecence brightness')+
  theme(axis.text.x = element_text(angle = 90,hjust=1))

dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Summary_Celegans_brightness_comparison.pdf",sep=''),
             width=5,height=5)









#Venn diagram

EC.ant<-as.character(subset(selectcast,`T-C_FDR`<0.05 & `T-C_logFC`>0 )$MetaboliteU)
EC.syn<-as.character(subset(selectcast,`T-C_FDR`<0.05 & `T-C_logFC`<0 )$MetaboliteU)
EC.non<-as.character(subset(selectcast,`T-C_FDR`>0.05 )$MetaboliteU)
EC.nonstrict<-as.character(subset(selectcast,`T-C_FDR`>0.05 & `C_FDR`>0.05 &`T_FDR`>0.05)$MetaboliteU)
Cel.rescue<-as.character(subset(selectcast,`Cel_FDR`<0.05 & -`Cel_logFC`>0 )$MetaboliteU)
Cel.aggr<-as.character(subset(selectcast,`Cel_FDR`<0.05 & -`Cel_logFC`<0 )$MetaboliteU)



rescue.non<-list('C. elegans aggravate'=Cel.aggr,
                 'E. coli non interacting'=EC.non,
                 'C. elegans rescue'=Cel.rescue)


rescue.nons<-list('C. elegans aggravate'=Cel.aggr,
                  'E. coli no effect'=EC.nonstrict,
                  'C. elegans rescue'=Cel.rescue)



intersect(EC.syn,Cel.rescue)
setdiff(Cel.rescue,EC.ant)



rescue<-list('C. elegans rescue'=Cel.rescue,
             'C. elegans aggravate'=Cel.aggr,
             'E. coli rescue'=EC.ant,
             'E. coli aggravate'=EC.syn)

vcols<-c("red","red4","blue","blue4")
grid.draw(VennDiagram::venn.diagram(rescue, col=vcols,cat.col=vcols ,NULL))

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Celegans_Ecoli_test.pdf",sep=''),
             width=8,height=8, useDingbats=FALSE)
dev.off()

#plot(Venn(rescue.non),show = list(Faces = FALSE),doWeights = TRUE)
grid.draw(VennDiagram::venn.diagram(rescue.non,NULL))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Celegans_Ecoli_noninteracting.pdf",sep=''),
             width=8,height=8, useDingbats=FALSE)
dev.off()


#plot(Venn(rescue.nons),show = list(Faces = FALSE),doWeights = TRUE)
grid.draw(VennDiagram::venn.diagram(rescue.nons, NULL))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Celegans_Ecoli_noeffect.pdf",sep=''),
             width=8,height=8, useDingbats=FALSE)
dev.off()


intersect(Cel.rescue,EC.syn)
intersect(Cel.rescue,EC.ant)


intersect(Cel.aggr,EC.ant)
intersect(Cel.aggr,EC.syn)

intersect(Cel.rescue,EC.non)
intersect(Cel.aggr,EC.non)


intersect(Cel.rescue,EC.nonstrict)
intersect(Cel.aggr,EC.nonstrict)



#Biolog class enrichemnt
#Enrichment by class

# head(results)
# 
# results.ec<-results.eco %>%
#   filter(Contrast=='T-C') %>%
#   select('Index','Plate','Metabolite','MetaboliteU','logFC','FDR')
# 
# results.ec$Organism<-'Ec'
# head(results.ec)
# dim(results.ec)
# 
# head(results.cel.r)
# 
# results.celegans<-results.cel.r[,c('Index','Plate','Metabolite','MetaboliteU','logFC','FDR')]
# results.celegans$logFC<- - results.celegans$logFC
# results.celegans$Organism<-'Ce'
# head(results.celegans)
# dim(results.celegans)
# results.enr<-rbind(results.ec,results.celegans)
#results.enr.f<-subset(results.enr,!(Plate %in% c('PM3B','PM4A') & Metabolite %in% PM1PM2$Metabolite ) & !Metabolite %in% c('Negative Control',coxy ) )


results.ecocel

#337

head(results.enr.f)

#Classes preparation
classes<-info %>%
  select('Index','Plate','Well','MetaboliteU') %>%
  left_join(readxl::read_xlsx('../Biolog/Biolog_metabolites_EcoCyc.xlsx',sheet = 'Classes') %>% select('Index','Class1','Class2')) %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A')) %>%
  gather(Redundancy,Class,contains("Class")) %>%
  filter(Class!="") 
  


dim(classes)


class.all<-data.frame(table(as.character(classes$Class))) %>%
  rename('Class'='Var1','Class_size'='Freq') 

#sum(classes$Class_size)
# class.all
# classes.c<-classes.m[,c('Index','Class')]
# classes.data<-merge(classes.c,results.enr.f,by='Index',all.x=TRUE)




classes.data<-results.ecocel%>%
  left_join(classes) %>%
  filter(!is.na(FDR)) %>%
  #select thresholds
  mutate(Ant=logFC>0 & FDR <0.05,
         Syn=logFC<0 & FDR <0.05,
         All=FDR <0.05) %>%
  #gather different thresholds
  gather(Type,Pass,Ant,Syn,All) %>%
  group_by(Organism,Type) %>%
  mutate(Total_pass=sum(Pass,na.rm = TRUE)) %>%
  group_by(Organism,Type,Class,Total_pass) %>%
  summarise(Class_size=n(),
         Class_pass=sum(Pass,na.rm = TRUE))



classes.data %>%
  group_by(Organism,Type,Total_pass) %>%
  summarise

classes.data %>%
  View()

# 
# int.counts<-ddply(classes.data,.(Class,Organism),summarise,Ant=sum(logFC>0 & FDR <0.05),Syn=sum(logFC<0 & FDR <0.05),Chg=sum(FDR <0.05))
# int.m<-melt(int.counts,measure.vars=c('Ant','Syn','Chg'),variable.name = 'Type',value.name = 'Count')
# int.m

# grp.counts<-ddply(int.m,.(Type,Organism),summarise,FreqGroup=sum(Count,na.rm = TRUE))
# grp.counts$Group<-paste(met.counts$Organism,met.counts$Type,sep='_')
# grp.counts

# class.counts<-ddply(int.m,.(Organism,Class),summarise,FreqGroup=sum(Count,na.rm = TRUE))
# class.counts



met.cast<-dcast(int.m,Class~Organism+Type,value.var = 'Count')
met.freq.t<-merge(class.all,met.cast,by='Class',all=TRUE)

met.freq.t

met.freq.m<-melt(met.freq.t,id.vars = c('Class','FreqAll'),variable.name = 'Group', value.name = 'Freq')
met.freq.m

met.sum<-merge(met.freq.m,grp.counts[,c('Group','FreqGroup')],by='Group')
met.sum

sum(class.all$FreqAll)
freq.All<-sum(met.freq.t$FreqAll)
freq.All


#phyper(q, m, n, k)
# pop size : 5260
# sample size : 131
# Number of items in the pop that are classified as successes : 1998
# Number of items in the sample that are classified as successes : 62
#phyper(62-1, 1998, 5260-1998, 131, lower.tail=FALSE)

met.sum$p_val<-phyper(met.sum$Freq-1,met.sum$FreqGroup,freq.All-met.sum$FreqGroup,met.sum$FreqAll,lower.tail =FALSE)

#met.sum$FDR<-p.adjust(met.sum$p_val,method = 'fdr')
met.sum$FC<- (met.sum$Freq/met.sum$FreqGroup)/(met.sum$FreqAll/freq.All)

met.sum
met.sum.m<-melt(met.sum,id.vars = c('Class','Group','FreqAll'),measure.vars = c('Freq','p_val','FC'),variable.name = 'Stat',value.name = 'Value')
met.cast<-dcast(met.sum.m,Class+FreqAll~Group+Stat,value.var = 'Value')

met.cast




write.csv(met.cast,paste(odir,'/Metabolite_class_enrichment_nocarboxy.csv',sep=''),row.names = FALSE)


subset(met.cast,Ce_Syn_p_val<0.05)


unique(as.character(selectcast$Metabolite_class))






#EcoCyc enrichment

ecli.co<-read.csv('../Biolog/Ecoli_All_class_instances_non_redundant.csv',header=TRUE,row.names = 1)

#EcoCyc to Biolog links
pmcs<-read.csv('../Biolog/EColi_net/All_Media_Mix_clean_explicit.csv',stringsAsFactors = FALSE)


head(ecli.co)


ecli.co.u<-ecli.co[!duplicated(ecli.co$Class),]
ecli.co.u[,c('Instance','IsClass')]<-NULL




#Enrichment analysis start
encast.r<-selectcast
encast.r$Res_logFC<- -encast.r$Cel_logFC
encast.r$Res_FDR<-encast.r$Cel_FDR

encast<-encast.r[,c('Plate','Well','Index','MetaboliteU','EcoCycID',
                    'T-C_logFC','T-C_FDR',
                    'Res_logFC','Res_FDR')]

colnames(encast)<-gsub('Res','Celegans',colnames(encast))
colnames(encast)<-gsub('T-C','Ecoli',colnames(encast))

ecdata<-merge(encast,pmcs[,c('Plate','Well','Instance')], by=c('Plate','Well'))

dim(ecdata)
head(ecdata)


ecli.co1<-subset(ecli.co,NoInstances>1)
ecli.co2<-subset(ecli.co,NoInstances>2)


ecli.cos<-subset(ecli.co,NoInstances>2 )



cdata<-merge(ecli.cos[,c('Class','NoInstances','subclasses','supclasses','Hierarchy','Instance')],ecdata,by='Instance')

head(cdata)
dim(cdata)
cdata.c<-cdata[!duplicated(cdata[,c('Class','Index','EcoCycID')]),]

#remove redundant assignments to class
dim(cdata.c)

write.csv(cdata.c,paste0(odir,'/EcoCyc_all_classes.csv'))




head(cdata.c)

idvariables<-c('Class','NoInstances','subclasses','supclasses','Hierarchy','Instance','Index','EcoCycID','MetaboliteU')
selstats<-c('logFC','FDR')

cdata.m<-enrichment.melt(cdata,idvariables,selstats)

unique(cdata.m$Contrast)


#remove(enrichmentmelt)

cdata.en<-enrichment(cdata.m,terms = c('Class','NoInstances','subclasses','supclasses','Hierarchy'),IDs = 'Instance',comparisons = 'Contrast',change = 'logFC',sign = 'FDR')


head(cdata.en)
subset(cdata.en,FDR<0.05)


write.csv(cdata.en,paste0(odir,'/EcoCyc_enrichment.csv'))


#Get data ready for heatmap
cdata.cast<-dcast(subset(cdata.en, Term_total>3 & Test!='All'),Class+Hierarchy~Contrast+Test,value.var = 'FDR')

comps<-c('Celegans_Up','Celegans_Down','Ecoli_Up','Ecoli_Down')


rownames(cdata.cast)<-cdata.cast$Class



#131 term
#Leave only enrichment terms with at least one significant value in a row
dim(cdata.cast)
cdata.cast<-cdata.cast[apply(cdata.cast[,comps], 1, function(x) { any(as.numeric(x) < 0.05,na.rm=TRUE)}),]
dim(cdata.cast)


cdata.cast[order(cdata.cast$Hierarchy),]



cdata.cast$Code<-apply(cdata.cast[,comps],1,function(x) paste0(as.numeric(x< 0.05),collapse = '') )


#write.csv(cdata.cast,paste0(odir,'/EcoCyc_enrichment_heatmap.csv'))


#Some manual filtering
EChm<-read.xlsx(paste0(odir,'/EcoCyc_enrichment_heatmap.xlsx'),sheetName = 'Heatmap',header = TRUE)



cdata.cast.f<-subset(cdata.cast,Class %in% EChm[is.na(EChm$Remove),'Class'])


cdata.cast.f<-cdata.cast.f[order(cdata.cast.f$Class),]





reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
gyrs<-colorRampPalette(c("gray90","steelblue1","blue4"))(n = 5)
enbrks<-c(0,-log(0.05,10),2,3,4,5)



hdata<-cdata.cast.f[,comps]


rownames(hdata)<-gsub('\\|','',rownames(hdata))


hdata<- -log10(hdata)

hdatafill<-hdata
hdatafill[is.na(hdatafill)]<-0
hdatafill[hdatafill==Inf]<-30

#Use the dummy dataset to get ordering of rows (enrichment terms)
#Graphical settings here are irrelevent for the final plot
#Data needs to be supplied as a matrix - table with row/column names defined implicitly and only numeric values
hmap<-heatmap.2(data.matrix(hdatafill),
                key=TRUE,
                Colv=FALSE,
                trace='none',
                col=gyrs,
                xlab='Comparison',
                Rowv=TRUE,
                dendrogram="row",
                scale="none",
                na.color="white", 
                symkey=FALSE,
                breaks=enbrks,
                reorderfun=reorderfun_mean,
                cexRow=0.8,
                cexCol=0.5)

#Use ordering with the real data
heatmap.2(data.matrix(hdata),
          key=FALSE, #Colour key needs to be provided separately
          Colv=FALSE, #Columns should not be reordered
          Rowv=hmap$rowDendrogram, # Here comes in the ordering
          trace='none',
          col=gyrs, # Colour scale
          breaks=enbrks, # Colour breaks
          xlab='Comparison',
          dendrogram='none', #Row dendogram, but should be changed to none, as dendrogram represents our data with filled-in values
          scale="none", #Should values be normalised in rows or columns - No
          na.color="white", # What colour to use with not missing values
          symkey=FALSE, #Provided colour scale is not symetrical
          cexRow=0.7, #Some figure scaling parameters. Works only after a lot of experimentation
          cexCol=0.7,
          margin=c(10,20),
          lwid=c(0.2,0.8),
          lhei=c(0.05,0.95))

dev.copy2pdf(device=cairo_pdf,file=paste0(odir,'/Enrichment_heatmap_EcoCyc_UpDown_clean.pdf'),
             width=5,height=5,useDingbats=FALSE)



#Enrichment tables
fcodes<-unique(cdata.cast.f$Code)




all.en.all$p <- format(all.en.all$p, scientific = TRUE,digits = 2)
all.en.all$q <- format(all.en.all$q, scientific = TRUE,digits = 2)



grid.table(cdata.cast.f)






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








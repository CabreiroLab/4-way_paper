#Figure numbering might have been changed
library(tidyverse)
library(broom)
library(ggrepel)


#devtools::install_github("PNorvaisas/PFun")
library(PFun)


theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau


cwd<-"~/Dropbox/Projects/Metformin_project/Bacterial Growth Assays/"
setwd(cwd)


odir<-'Summary_Arg_mutants'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#save.image('BGA_Arg_mutants.RData')
#load('BGA_Arg_mutants.RData')



strainlist<-c('OP50-C','astA','adiA','speA','speB','adiAspeA','astAspeB','speAspeB','astAadiAspeA','astAadiAspeAspeB')
strainlabels<-c('OP50-C','astA','adiA','speA','speB','adiAspeA','astAspeB','speAspeB','TM','QM')



data<-read_csv('arginine degradation mutants/Data/Summary.csv') %>%
  select(Plate:TReplicate,logAUC=Int_600nm_log,Dph=a_log,-Data) %>%
  filter(!is.na(Strain) ) %>%
  mutate_at(vars(Plate:TReplicate),as.factor) %>%
  mutate(Strain=factor(Strain,levels=strainlist,labels=strainlabels),
         logDph=log2(Dph)) %>%
  select(-Dph) %>%
  gather(Measure,Value,logDph,logAUC) %>%
  group_by(Measure) %>%
  mutate(Ref_C=mean(Value[Strain=='OP50-C']),
         Norm_C=Value-Ref_C) %>%
  ungroup %>%
  select(-Ref_C) %>%
  gather(Normalisation,Value,Value,Norm_C) %>%
  #left_join(NormNames) %>%
  mutate(Value=ifelse(Value %in% c(Inf,-Inf),NA,Value ))

data.sum<-data %>%
  group_by(Measure,Normalisation,Strain) %>%
  summarise(Mean=mean(Value),
            SD=sd(Value),
            SE=SD/sqrt(n())) %>%
  mutate(NE=Mean-SE,
         PE=Mean+SE,
         Prc=ifelse(Normalisation=="Value",NA,2^Mean*100),
         PrcNE=ifelse(Normalisation=="Value",NA,2^NE*100),
         PrcPE=ifelse(Normalisation=="Value",NA,2^PE*100))



data %>%
  write_csv(paste0(odir,"/Raw_data_Summary.csv"))


data_ts<-read_csv('arginine degradation mutants/Data/Data.csv') %>%
  filter(!is.na(Strain) ) %>%
  filter(Data=='600nm_f') %>%
  gather(Time_s,OD,contains('.0')) %>%
  mutate_at(vars(Plate:TReplicate),as.factor) %>%
  select(-Data) %>%
  mutate(Strain=factor(Strain,levels=strainlist,labels=strainlabels),
         Time_s=as.numeric(Time_s),
         Time_h=Time_s/3600)



ggplot(data_ts,aes(x=Time_h,y=OD,color=Strain))+
  geom_line(aes(group=interaction(Replicate,TReplicate,Strain)))+
  xlab('Time, h')+
  scale_x_continuous(breaks=seq(0,24,by=6))+
  facet_wrap(~Strain,ncol=5)


ggsave(file=paste(odir,"/Growth_overview.pdf",sep=''),
             width=110,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")





datats.sum<-data_ts %>%
  group_by(Strain,Time_h) %>%
  summarise(OD_Mean=mean(OD),
            OD_SD=sd(OD))


datats.sum %>%
  ggplot(aes(x=Time_h,y=OD_Mean,color=Strain,fill=Strain))+
  geom_line()+
  geom_ribbon(aes(ymin=OD_Mean-OD_SD,
                  ymax=OD_Mean+OD_SD),alpha=0.5,color=NA)+
  xlab('Time, h')+
  ylab('OD')+
  scale_x_continuous(breaks=seq(0,18,by=6))
+
  facet_wrap(~Strain,ncol=5)

ggsave(file=paste0(odir,"/Growth_Summary_tiny.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



stats<-data %>%
  filter(Normalisation!="Value") %>%
  group_by(Measure,Normalisation) %>%
  do(tidy(lm(Value~Strain,data=.))) %>%
  filter(term!="(Intercept)") %>%
  rename(Strain=term,logFC=estimate,SE=std.error,t.value=statistic,p=p.value) %>%
  mutate(Comparison="Strain",
         Strain=str_replace(Strain,"Strain",""),
         NE=logFC-SE,
         PE=logFC+SE,
         Prc=ifelse(Normalisation=="Value",NA,2^logFC*100),
         PrcNE=ifelse(Normalisation=="Value",NA,2^NE*100),
         PrcPE=ifelse(Normalisation=="Value",NA,2^PE*100),
         pStars=pStars(p),
         Strain=factor(Strain,levels=strainlabels)) %>%
  select(Measure,Normalisation,Comparison,Strain,everything())



stats %>%
  write_csv(paste0(odir,"/All_results.csv"))


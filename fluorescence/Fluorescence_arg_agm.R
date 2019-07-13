library(tidyverse)
library(scales)
library(broom)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)



setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy/")

odir<-'Summary_ArgAgmGly'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")




#load("Fluorescence_ArgAgmGly.RData")
#save.image('Fluorescence_ArgAgmGly.RData')



#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau




condinfo<-read_csv("Conditions_all_ArgAgmGly.csv") %>%
  mutate_all(as.factor)


removals<-c("Agm_QM_50_2_12",
            "Agm_QM_0_1_2",
            "Agm_QM_0_1_9",
            "Agm_QM_0_1_5",
            "Agm_speB_50_1_4",
            "Agm_OP50-C_10_1_3",
            "Agm_OP50-C_0_2_4",
            "Agm_OP50-C_10_1_3",
            "Agm_OP50-C_25_2_1",
            "Agm_OP50-C_50_2_11")






#Translations for consistency

strans<-c("op50-c"="OP50-C")


data.all<-data.frame(Type=c("Agm","Arg","Gly"),File=c("Agmatine titration-OP50-C-speB-QM.txt","arginine mutants.txt","Glycine-serine mutants.txt")) %>%
  group_by(Type,File) %>%
  do(read_delim(paste('./Data/Arginine mutants/',.$File,sep='/'),delim=ifelse(str_detect(.$File,fixed('.csv')),',','\t')) %>% mutate_all(funs(str_replace(.,"\\*","")) ) %>% gather(Condition,Abs,everything()) ) %>%
  group_by(File,Type,Condition) %>%
  mutate(Abs=as.numeric(Abs),
         Rank=row_number(),
         Replicate=case_when(Rank>=31 ~ 2,
                             Rank>64 ~ 3,
                             TRUE ~ 1) %>% as.character(),
         Strain=str_trim(Condition) %>% str_split(' ') %>% unlist %>% head(1),
         Strain=ifelse(Strain %in% names(strans),strans[Strain],Strain),
         Agmatine_mM=str_extract(Condition,"[0-9]{1,2}.?mM") %>% str_replace('mM',"") %>% as.numeric ) %>%
  ungroup %>%
  filter(!is.na(Abs) ) %>%
  bind_rows( filter(.,Type=="Gly" | (Type=="Arg" & Strain %in% c("OP50-C","speB") ) ) %>%
              mutate(Replicate=ifelse(Type=="Gly",3,Replicate),
                     Type="ArgGly") %>%
               data.frame ) %>%
  mutate(Log=log2(Abs)) %>%
  mutate_at(c('Type','Strain','Agmatine_mM'),as.factor) %>%
  gather(Measure,Value,Abs,Log) %>%
  select(Type,Condition,Replicate:Value) %>%
  group_by(Measure,Type,Strain,Condition,Replicate) %>%
  mutate(Worm=row_number()) %>%
  unite(ID,Type,Strain,Agmatine_mM,remove = FALSE) %>%
  unite(WID,ID,Replicate,Worm,remove = FALSE) %>%
  filter(!WID %in% removals) %>%
  #Internal normalisation for each transgene
  group_by(Type,Measure) %>%
  mutate(Ref=ifelse(Type=="ArgGly",mean(Value[Strain=='OP50-C' & Replicate %in% c('1','2')]),mean(Value[Strain=='OP50-C']) ) ) %>%
  group_by(Type,Replicate,Measure) %>%
  mutate(RefR=mean(Value[Strain=='OP50-C']),
         Norm=Value-RefR+Ref,
         NormAbs=2^Norm) 


data.all %>%
  filter(Type =="ArgGly") %>%
  pull(Replicate)




unique(data.all$ID)

data.all %>%
  write_csv(paste0(odir,"/All_data_raw.csv"))

View(data.all)

data.all %>%
  group_by(Type,Condition,Strain,Agmatine_mM) %>%
  summarise(N=n()) %>%
  View




head(data.all)



data.all %>%
  filter(Measure=='Log') %>%
  group_by(Type,Replicate,Strain,Agmatine_mM) %>%
  summarise(Count=n()) %>%
  write_csv(paste0(odir,'/All_data_Worm_count.csv') )



data.all %>%
  group_by(Type,Condition,Strain,Agmatine_mM) %>%
  summarise %>%
  write_csv('Conditions_raw_all_ArgAgmGly.csv')



data.all %>%
  filter(is.na(ID))




#Variability
sum.c<-data.all %>%
  filter(Measure=='Log') %>%
  group_by(ID,Strain,Agmatine_mM) %>%
  summarise(SD=sd(Value,na.rm = TRUE),Mean=mean(Value,na.rm = TRUE)) %>%
  mutate(Index=paste(ID),
         VarPrc=ifelse(is.na(SD) ,Inf, (2^(SD)-1)*100 ) ) %>%
  arrange(VarPrc) %>%
  data.frame %>%
  mutate(Index=factor(Index, levels=Index,labels=Index))


sum.c %>%
  filter(VarPrc>300)


ggplot(sum.c,aes(x=Index,y=VarPrc))+
  geom_point()+
  scale_y_continuous(breaks = seq(0,1000,by=25),limits=c(0,100))+
  ylab('In-group variation, %')+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste0(odir,"/Ingroup_variation_Percentage.pdf"),
             width=7,height=10, useDingbats=FALSE)





#Agmatine
data.agm<-data.all %>%
  filter(Type=="Agm" & Measure=="Log") %>%
  mutate(Strain=factor(Strain,levels=c("OP50-C","speB","QM") ) ) %>%
  ungroup


data.agm %>%
  ggplot(aes(x=Agmatine_mM,y=NormAbs,color=Strain))+
  geom_jitter(width=0.25)+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x))) + 
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  ylab('Mean fluorescence per worm, a.u.')+
  xlab("Agmatine, mM")+
  #geom_text(aes(label=paste(Replicate,Worm,sep="_")))+
  facet_grid(.~Strain) 

ggsave(file=paste0(odir,"/Agmatine_Boxplot.pdf"),
       device=cairo_pdf,family="Arial",
       scale=2,
       width=110,height=41,units='mm')



stat.agm.S <- data.agm %>%
  group_by(Agmatine_mM) %>%
  do(tidy(lm(Norm~Strain,data=.))) %>%
  filter(term!="(Intercept)") %>%
  rename(Strain=term) %>%
  mutate(Strain=str_replace(Strain,"Strain",""),
         Comparison="Strain") %>%
  ungroup


stat.agm.A <-data.agm %>%
  group_by(Strain) %>%
  do(tidy(lm(Norm~Agmatine_mM,data=.))) %>%
  filter(term!="(Intercept)") %>%
  rename(Agmatine_mM=term) %>%
  mutate(Agmatine_mM=str_replace(Agmatine_mM,"Agmatine_mM",""),
         Comparison="Agmatine") %>%
  ungroup


stat.agm.I <- data.agm %>%
  do(tidy(lm(Norm~Strain*Agmatine_mM,data=.))) %>%
  filter(str_detect(term,":")) %>%
  mutate(term=str_replace_all(term,"Strain|Agmatine_mM","")) %>%
  separate(term,c("Strain","Agmatine_mM"),sep=":") %>%
  mutate(Comparison="Interaction") %>%
  ungroup


stat.agm <- stat.agm.S %>%
  rbind(stat.agm.A) %>%
  rbind(stat.agm.I) %>%
  rename(logFC=estimate,SE=std.error,t.value=statistic,p=p.value) %>%
  mutate(pStars=pStars(p),
         Strain=factor(Strain,levels=c("OP50-C","speB","QM") )) %>%
  select(Comparison,everything())
  


stat.agm %>%
  write_csv(paste0(odir,"/Results_Agmatine.csv"))


View(stat.agm)

sum.agm<-data.agm %>%
  group_by(Strain,Agmatine_mM) %>%
  summarise(Mean=mean(Norm),
            SD=sd(Norm),
            SE=SD/sqrt(n()),
            PD=Mean+SD,
            ND=Mean-SD,
            PE=Mean+SE,
            NE=Mean-SE,
            MeanAbs=2^Mean,
            SDAbs=2^SD,
            SEAbs=2^SE,
            PEAbs=2^PE,
            NEAbs=2^NE,
            PDAbs=2^PD,
            NDAbs=2^ND)



quartz()
#hj<-0.8
vj<-2
nx<--0.2


sum.agm %>%
  ggplot(aes(x=Agmatine_mM,y=MeanAbs,color=Strain,fill=Strain))+
  geom_ribbon(aes(ymin=NDAbs,ymax=PDAbs,group=Strain),color=NA,alpha=0.5)+
  
  geom_line(aes(group=Strain))+
  geom_point(size=5)+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x)),
                     limits=c(2^7,2^12.5)) + 
  geom_text(data=stat.agm %>% filter(Comparison=="Interaction" ),aes(label=pStars,y=2^12,color=Strain),show.legend = FALSE,size=5,angle=45)+
  #geom_text(data=stat.agm %>% filter(Comparison=="Strain" ),aes(label=pStars,y=2^(11.5+as.numeric(Strain)*0.5),color=Strain),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj,angle=45)+
  ylab('Mean Pacs-2::gfp fluorescence per worm, a.u.')+
  xlab("Agmatine, mM")


ggsave(file=paste0(odir,"/Agmatine_Summary.pdf"),
       device=cairo_pdf,family="Arial",
       scale=2,
       width=55,height=41,units='mm')





#Arginine

argstrains<-c("OP50-C","astA","adiA","speA","speB","adiAspeA","speAspeB","astAspeB","astAadiAspeA","astAadiAspeAspeB")

data.arg<-data.all %>%
  filter(Type=="Arg" & Measure=="Log") %>%
  mutate(Strain=factor(Strain,levels=argstrains)) %>%
  ungroup



stat.arg <-data.arg %>%
  do(tidy(lm(Norm~Strain,data=.))) %>%
  rename(Strain=term) %>%
  filter(Strain!="(Intercept)") %>%
  mutate(Strain=str_replace_all(Strain,"Strain","")) %>%
  rename(logFC=estimate,SE=std.error,t.value=statistic,p=p.value) %>%
  mutate(pStars=pStars(p),
         Strain=factor(Strain,levels=argstrains )) 

stat.arg %>%
  write_csv(paste0(odir,"/Results_Arginine_degradation.csv"))



quartz()
#hj<-0.8
vj<-2
nx<--0.2



data.arg %>%
  ggplot(aes(x=Strain,y=NormAbs))+
  geom_jitter(width=0.25)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x)),
                     limits=c(2^7,2^12.5)) + 
  ylab('Mean Pacs-2::gfp fluorescence per worm, a.u.')+
  geom_text(data=stat.arg,aes(label=pStars,y=2^12),show.legend = FALSE,size=5,angle=45)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggsave(file=paste0(odir,"/Arginine_degradation_Summary.pdf"),
       device=cairo_pdf,family="Arial",
       scale=2,
       width=55,height=41,units='mm')




#Glycine

glystrains<-c("OP50-C","speB","argG","speF","gcvH","gcvP","gcvT","tdcB","tdcC","hpt")
data.gly<-data.all %>%
  filter(Type=="ArgGly" & Strain %in% glystrains & Measure=="Log") %>%
  mutate(Strain=factor(Strain,levels=glystrains)) %>%
  ungroup


stat.gly<-data.gly %>%
  do(tidy(lm(Norm~Strain,data=.))) %>%
  rename(Strain=term,logFC=estimate,SE=std.error,t.value=statistic,p=p.value) %>%
  filter(Strain!="(Intercept)") %>%
  mutate(Strain=str_replace_all(Strain,"Strain","")) %>%
  mutate(pStars=pStars(p),
         Strain=factor(Strain,levels=glystrains )) 

stat.gly %>%
  write_csv(paste0(odir,"/Results_Other.csv"))


data.gly %>%
  View




data.gly %>%
  ggplot(aes(x=Norm,fill=Strain))+
  geom_histogram()


quartz()
#hj<-0.8
vj<-2
nx<--0.2


data.gly %>%
  ggplot(aes(x=Strain,y=NormAbs))+
  geom_jitter(width=0.25)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x)),
                     limits=c(2^7,2^12)) +
  ylab('Mean Pacs-2::gfp fluorescence per worm, a.u.')+
  geom_text(data=stat.gly,aes(label=pStars,y=2^12),show.legend = FALSE,nudge_y = -0.2,size=5,angle=45)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#

ggsave(file=paste0(odir,"/Other_Summary.pdf"),
       device=cairo_pdf,family="Arial",
       scale=2,
       width=55,height=41,units='mm')









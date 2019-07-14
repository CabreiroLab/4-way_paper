# Script used for several figures. Figure numbering might have changed.

library(tidyverse)
library(PFun)
library(broom)

#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau


Metcols <- c("#FF0000","#32006F")
names(Metcols) <- c("0","50")
Metlab<-'Metformin, mM'

setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy")


odir<-'Summary_Fluorescence_media'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

#load('Fluorescence_media.RData')
#save.image('Fluorescence_media.RData')

translations<-c('op50-c'='OP50-C', 'op50-mr'='OP50-MR')



#2j nhr-49 acs-2 fluor
wtrans<-c("nhr-49"="nhr-49 (nr2041)")

f2j<-read_delim('Data/redo of figures/2j_acs-2nhr-49.txt',delim="\t") %>%
  gather(Condition,Abs,everything()) %>%
  filter(!is.na(Abs)) %>%
  mutate(Condition=str_trim(Condition)) %>%
  separate(Condition,c("Worm","Metformin_mM"),sep = ':') %>%
  mutate(Log=log2(Abs),
         Worm=ifelse(Worm %in% names(wtrans),wtrans[Worm],Worm),
         Gene="acs-2::gfp",
         Strain="OP50-C") 


sf2j_tr<-f2j %>%
  group_by(Worm) %>%
  filter(!is.na(Log))%>%
  do(tidy(lm(Log~Metformin_mM,data=.))) %>%
  filter(term!="(Intercept)") %>%
  mutate(Metformin_mM=str_replace(term,'Metformin_mM',''),
         Contrast_type="Treatment",
         term=NULL)

sf2j_int<-f2j %>%
  filter(!is.na(Log)) %>%
  do(tidy(lm(Log~Metformin_mM*Worm,data=.))) %>%
  ungroup %>%
  filter(str_detect(term,':') & term!="(Intercept)" ) %>%
  mutate(term=str_replace_all(term,'Metformin_mM|Worm',''),
         Contrast_type="Interaction") %>%
  separate(term,c('Metformin_mM','Worm'),sep=':') 

sf2j<-sf2j_tr %>%
  ungroup %>%
  bind_rows(sf2j_int) %>%
  select(Contrast_type,Worm,Metformin_mM,everything()) %>%
  mutate_at(c("Metformin_mM",'Worm'),as.factor) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0","50")),
         pStars=pStars(p.value)) %>%
  rename(SE=std.error,
         logFC=estimate)



sf2j %>%
  write_csv(paste0(odir,'/2j_N2_nhr-49_Summary.csv'))


blank_data_f2j<-f2j %>%
  group_by(Worm) %>%
  summarise(Abs=2^(min(Log)+(max(Log)-min(Log))*1.3 ) )%>% 
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0",'50')))


hj<-0.8
vj<-2
nx<-0#-0.25


#,position = position_jitterdodge(dodge.width=1,jitter.width=0.4)
#,position = position_dodge(width=1)


quartz()
f2j %>%
  ggplot(aes(x=Worm,y=Abs,color=Metformin_mM))+
  geom_blank(data = blank_data_f2j,aes(x=Worm,y=Abs))+
  geom_jitter(aes(fill=Metformin_mM),size=1,alpha=0.5,width=0.25)+
  #geom_boxplot(data=bdata,aes(ymax=Max,ymin=Min,lower=Mean-SD,upper=Mean+SD,middle=Mean),position="identity",stat = "identity",alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot", alpha=0.5) +
  scale_y_continuous(trans = 'log2',
                     breaks = scales::trans_breaks('log2', function(x) 2^x),
                     labels = scales::trans_format('log2', scales::math_format(2^.x))) + 
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  scale_fill_manual(name = "Metformin, mM",values =Metcols)+
  geom_text(data=sf2j %>% filter(Contrast_type=="Treatment") ,aes(label=pStars,y=2^13),show.legend = FALSE,size=5)+
  geom_text(data=sf2j %>% filter(Contrast_type=="Interaction") ,aes(label=pStars,y=2^14),show.legend = FALSE,size=5,color="green4")+
  ylab('Mean fluorescence per worm, a.u.')+
  facet_grid(Gene~Strain,space="free_x",scales="free_x")+
  theme(legend.position="top")

ggsave(file=paste0(odir,'/2h_Fluorescence_N2_nhr-49_logScale2_tiny.pdf'),
       width=25,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")





#2k nhr-49 acs-2 fluor wha-6 mCherry fluorescence
strans<-c("nhr-49"="OP50-C","op50-mr"="OP50-MR","op50-c"="OP50-C")

f2k<-read_delim('Data/redo of figures/2k_vha-6nhr-49.txt',delim="\t") %>%
  gather(Condition,Abs,everything()) %>%
  filter(!is.na(Abs)) %>%
  mutate(Condition=str_trim(Condition)) %>%
  separate(Condition,c("Strain","Metformin_mM"),sep = ' ') %>%
  mutate(Log=log2(Abs),
         Worm=ifelse(Strain=="nhr-49","nhr-49 (nr2041)","N2"),
         Strain=ifelse(Strain %in% names(strans),strans[Strain],Strain),
         Gene="vha-6::mCherry") 

unique(f2k$Strain)
unique(f2k$Worm)


sf2k_tr<-f2k %>%
  group_by(Worm,Strain) %>%
  filter(!is.na(Log))%>%
  do(tidy(lm(Log~Metformin_mM,data=.))) %>%
  filter(term!="(Intercept)") %>%
  mutate(Metformin_mM=str_replace(term,'Metformin_mM',''),
         Contrast_type="Treatment",
         term=NULL) 

sf2k_intS<-f2k %>%
  filter(Worm=="N2") %>%
  group_by(Worm) %>%
  filter(!is.na(Log)) %>%
  do(tidy(lm(Log~Metformin_mM*Strain,data=.))) %>%
  filter(str_detect(term,':') & term!="(Intercept)" ) %>%
  mutate(term=str_replace_all(term,'Metformin_mM|Strain',''),
         Contrast_type="Interaction_Strain") %>%
  separate(term,c('Metformin_mM','Strain'),sep=':') 
  
sf2k_intW<-f2k %>%
  filter(Strain=="OP50-C") %>%
  group_by(Strain) %>%
  filter(!is.na(Log)) %>%
  do(tidy(lm(Log~Metformin_mM*Worm,data=.))) %>%
  filter(str_detect(term,':') & term!="(Intercept)" ) %>%
  mutate(term=str_replace_all(term,'Metformin_mM|Worm',''),
         Contrast_type="Interaction_Worm") %>%
  separate(term,c('Metformin_mM','Worm'),sep=':') 
  
  

sf2k<-sf2k_tr %>%
  bind_rows(sf2k_intS) %>%
  bind_rows(sf2k_intW) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0","50")),
         pStars=pStars(p.value))%>%
  rename(SE=std.error,
         logFC=estimate) %>%
  select(Contrast_type,Worm,Strain,Metformin_mM,everything()) 


sf2k %>%
  write_csv(paste0(odir,'/2k_vhs-6_N2_nhr-49_Summary.csv'))


blank_data_f2k<-f2k %>%
  group_by(Worm,Strain) %>%
  summarise(Abs=2^(min(Log)+(max(Log)-min(Log))*1.3 ) )%>% 
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0",'50')))


hj<-0.8
vj<-2
nx<-0#-0.25

#,position = position_jitterdodge(dodge.width=1,jitter.width=0.4)
#,position = position_dodge(width=1)

quartz()
f2k %>%
  ggplot(aes(x=Worm,y=Abs,color=Metformin_mM))+
  geom_blank(data = blank_data_f2k,aes(x=Worm,y=Abs))+
  geom_jitter(aes(fill=Metformin_mM),size=1,alpha=0.5,width=0.25)+
  #geom_boxplot(data=bdata,aes(ymax=Max,ymin=Min,lower=Mean-SD,upper=Mean+SD,middle=Mean),position="identity",stat = "identity",alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot", alpha=0.5) +
  scale_y_continuous(trans = 'log2',
                     breaks = scales::trans_breaks('log2', function(x) 2^x),
                     labels = scales::trans_format('log2', scales::math_format(2^.x))) + 
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  scale_fill_manual(name = "Metformin, mM",values =Metcols)+
  geom_text(data=sf2k %>% filter(Contrast_type=="Treatment") ,aes(label=pStars,y=2^10.25),show.legend = FALSE,size=5)+
  geom_text(data=sf2k %>% filter(Contrast_type!="Treatment") ,aes(label=pStars,y=2^10.5),show.legend = FALSE,size=5,color="green4")+
  #ylab('Mean vha-6::mCherry fluorescence per worm, a.u.')+
  facet_grid(Gene~Strain,space="free_x",scales="free_x")+
  theme(legend.position="top",
        axis.title.y = element_blank())


ggsave(file=paste0(odir,'/2k_Fluorescence_N2_nhr-49_logScale2_tiny.pdf'),
       width=30,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")




#3i nhr-49 acs-2 fluor wha-6 mCherry fluorescence
mtrans<-c("Bacto"="Bacto peptone","Soy"="Soy peptone")

medias<-c("Bacto peptone","Soy peptone","LB","MRS")

f3i<-read_delim('Data/redo of figures/3i_acs-2media.txt',delim="\t") %>%
  gather(Condition,Abs,everything()) %>%
  filter(!is.na(Abs)) %>%
  mutate(Condition=str_trim(Condition)) %>%
  separate(Condition,c("Media","Metformin_mM"),sep = ' ') %>%
  mutate(Log=log2(Abs),
         Media=ifelse(Media %in% names(mtrans),mtrans[Media],Media),
         Media=factor(Media,levels=medias,labels=medias),
         Metformin_mM=as.factor(Metformin_mM)) 


sf3i<-f3i %>%
  group_by(Media) %>%
  filter(!is.na(Log))%>%
  do(tidy(lm(Log~Metformin_mM,data=.))) %>%
  ungroup %>%
  filter(term!="(Intercept)") %>%
  mutate(Metformin_mM=str_replace(term,'Metformin_mM',''),
         Contrast_type="Treatment",
         term=NULL,
         Metformin_mM=factor(Metformin_mM),
         pStars=pStars(p.value)) %>%
  rename(SE=std.error,
         logFC=estimate) %>%
  select(Contrast_type,Media,Metformin_mM,everything()) 


sf3i %>%
  write_csv(paste0(odir,'/3i_Media_Summary.csv'))


blank_data_f3i<-f3i %>%
  group_by(Metformin_mM) %>%
  summarise(Abs=2^(min(Log)+(max(Log)-min(Log))*1.3 ) )%>% 
  mutate(Media="LB",
         Media=factor(Media,levels=medias))


hj<-0.8
vj<-2
nx<-0#-0.25

sumf3i<-f3i %>%
  group_by(Media,Metformin_mM) %>%
  summarise(Mean=mean(Abs),
            SD=sd(Abs),
            PD=Mean+SD,
            ND=Mean-SD)


Medcols <- c("red","blue4", viridis::viridis(4,option="magma"))
names(Medcols) <- medias
Medlab<-"Media"



sumf3i %>%
  ggplot(aes(x=Metformin_mM,y=Mean,fill=Media,color=Media))+
  geom_ribbon(aes(group=interaction(Media),ymin=ND,ymax=PD),alpha=0.5,show.legend=FALSE,color=NA)+
  geom_point()+
  geom_line(aes(group=interaction(Media)))+
  scale_colour_manual(name = Medlab,values = Medcols)+
  scale_fill_manual(name = Medlab,values = Medcols)+
  scale_y_continuous(trans = 'log2',
                     breaks = scales::trans_breaks('log2', function(x) 2^x),
                     labels = scales::trans_format('log2', scales::math_format(2^.x))) + 
  geom_text(data=sf3i %>% filter(Contrast_type=="Treatment") ,aes(label=pStars,y=as.numeric(Media)*4 ),nudge_y = 2^4.9,show.legend = FALSE,size=5)+
  ylab('Mean vha-6::mCherry fluorescence per worm, a.u.')+
  xlab("Metformin, mM")

ggsave(file=paste0(odir,'/3i_Fluorescence_Medias_logScale2.pdf'),
       width=65,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



#4e nhr-49 acs-2 fluor wha-6 mCherry fluorescence

strans<-c("Arabinose"="Ara","Glycerol"="Glyc","D-Ribose"="Rib")

sglvls<-c("OP50-C","OP50-C+Glycerol","glpK","glpK+Glycerol",
          "OP50-C+Arabinose","fucI","fucI+Arabinose",
          "OP50-C+D-Ribose","rbsK","rbsK+D-Ribose")

f4<-data.frame(File=list.files("Data/redo of figures/",pattern="4e_*|s4c_*|s4f_*")) %>%
  group_by(File) %>%
  do(read_delim(paste0('Data/redo of figures/',.$File),delim="\t") %>%
       gather(Condition,Abs,everything())) %>%
  filter(!is.na(Abs)) %>%
  mutate(Condition=str_trim(Condition),
         Condition=str_replace_all(Condition,"mM",""),
         Replicate=ifelse(str_detect(File,"rep1"),1,2),
         Log=log2(Abs),
         Figure=str_sub(File, 1, str_length(File)-8),
         Gene="acs-2") %>%
  separate(Condition,c("SM","Supplement"),sep="\\+") %>% 
  separate(SM,c("Strain","Metformin_mM"),sep=" ") %>%
  mutate_at(c("Strain","Supplement"),str_trim) %>%
  mutate(Supplement=ifelse(is.na(Supplement),"",Supplement),
         Suppl=ifelse(Supplement %in% names(strans),strans[Supplement],Supplement),
         SGroup=ifelse(Supplement=="",Strain,paste(Strain,Supplement,sep="+")),
         SGroup=factor(SGroup,levels=sglvls),
         Metformin_mM=factor(Metformin_mM,levels=c(0,50)),
         Strain=factor(Strain,levels=c("OP50-C","glpK","fucI","rbsK")),
         Supplement=factor(Supplement,levels=c("","Glycerol","Arabinose","D-Ribose")),
         Suppl=factor(Suppl,levels=c("","Glyc","Ara","Rib"))) %>%
  group_by(Figure) %>%
  mutate(Ref=mean(Log[Strain=='OP50-C' & Supplement==""]) ) %>%
  group_by(Figure,Replicate) %>%
  mutate(Norm=Log-mean(Log[Strain=='OP50-C' & Supplement==""])+Ref,
         NormAbs=2^Norm)


sf4_T<-f4 %>%
  filter(!is.na(Norm))%>%
  group_by(Figure,Gene,Strain,Supplement,SGroup) %>%
  do(tidy(lm(Norm~Metformin_mM,data=.))) %>%
  ungroup %>%
  filter(term!="(Intercept)") %>%
  mutate(Metformin_mM=str_replace(term,'Metformin_mM',''),
         Contrast_type="Treatment",
         term=NULL) 
  
sf4_Su<-f4 %>%
  filter(!is.na(Norm))%>%
  group_by(Figure,Gene,Metformin_mM,Strain) %>%
  do(tidy(lm(Norm~Supplement,data=.))) %>%
  ungroup %>%
  filter(term!="(Intercept)") %>%
  mutate(Supplement=str_replace(term,'Supplement',''),
         Contrast_type="Supplement",
         term=NULL,
         SGroup=NA)


sf4_St<-f4 %>%
  filter(!is.na(Norm))%>%
  group_by(Figure,Gene,Metformin_mM,Supplement) %>%
  do(tidy(lm(Norm~Strain,data=.))) %>%
  ungroup %>%
  filter(term!="(Intercept)") %>%
  mutate(Strain=str_replace(term,'Strain',''),
         Contrast_type="Strain",
         term=NULL,
         SGroup=NA) 


sf4_StIM<-f4 %>%
  filter(!is.na(Norm))%>%
  group_by(Figure,Gene,Supplement) %>%
  do(tidy(lm(Norm~Metformin_mM*Strain,data=.))) %>%
  ungroup %>%
  filter(str_detect(term,':') & term!="(Intercept)" ) %>%
  mutate(term=str_replace_all(term,'Metformin_mM|Strain',''),
         Contrast_type="Interaction_Strain",
         SGroup=NA) %>%
  separate(term,c('Metformin_mM','Strain'),sep=':')

sf4_SuIM<-f4 %>%
  filter(!is.na(Norm))%>%
  group_by(Figure,Gene,Strain) %>%
  do(tidy(lm(Norm~Metformin_mM*Supplement,data=.))) %>%
  ungroup %>%
  filter(str_detect(term,':') & term!="(Intercept)" ) %>%
  mutate(term=str_replace_all(term,'Metformin_mM|Supplement',''),
         Contrast_type="Interaction_Supplement",
         SGroup=NA) %>%
  separate(term,c('Metformin_mM','Supplement'),sep=':') 


sf4_StISu<-f4 %>%
  filter(!is.na(Norm))%>%
  group_by(Figure,Gene,Metformin_mM) %>%
  do(tidy(lm(Norm~Strain*Supplement,data=.))) %>%
  ungroup %>%
  filter(str_detect(term,':') & term!="(Intercept)" ) %>%
  mutate(term=str_replace_all(term,'Strain|Supplement',''),
         Contrast_type="Interaction_StrainSupplement",
         SGroup=NA) %>%
  separate(term,c('Strain','Suppl'),sep=':') 
  
  
sf4<-sf4_T %>%
  bind_rows(sf4_Su) %>%
  bind_rows(sf4_St) %>%
  bind_rows(sf4_SuIM) %>%
  bind_rows(sf4_StIM) %>%
  bind_rows(sf4_StISu) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0","50")),
         SGroup=ifelse(Supplement=="",Strain,paste(Strain,Supplement,sep="+")),
         SGroup=factor(SGroup,levels=sglvls),
         pStars=pStars(p.value))%>%
  rename(SE=std.error,
         logFC=estimate) %>%
  select(Figure,Contrast_type,Gene,Strain,SGroup,Metformin_mM,everything()) 
  
sf4 %>%
  write_csv(paste0(odir,'/4ecf_acs-2_Supplements_Summary.csv'))



plot_boxstats<-function(data,allstats) {
  
  hj<-0.8
  vj<-2
  nx<--0.1
  
  stats<-allstats %>%
    filter(Figure==unique(data$Figure))
  
  blank_dataf<-f4 %>%
    filter(Figure==unique(data$Figure)) %>%
    group_by(Figure,Gene,Strain,Supplement,Suppl,SGroup) %>%
    summarise(NormAbs=2^(min(Norm)+(max(Norm)-min(Norm))*1.2 ) ) %>%
    mutate(Metformin_mM=0,
           Metformin_mM=factor(Metformin_mM,levels=c("0","50")))
  
  data %>%
    ggplot+
    aes(x=SGroup,y=NormAbs,color=Metformin_mM)+
    geom_jitter(aes(fill=Metformin_mM),width=0.25,size=1,alpha=0.5)+
    stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
    geom_blank(data = blank_dataf,aes(x=SGroup,y=NormAbs),show.legend=FALSE)+
    scale_y_continuous(trans = 'log2',
                       breaks = scales::trans_breaks('log2', function(x) 2^x),
                       labels = scales::trans_format('log2', scales::math_format(2^.x))) + 
    scale_colour_manual(name = "Metformin, mM",values =Metcols)+
    scale_fill_manual(name = "Metformin, mM",values =Metcols)+
    geom_text(data=stats %>% filter(Contrast_type=="Treatment"),aes(label=pStars,y=Inf),show.legend = FALSE,size=5, vjust=vj+0.5)+ #,nudge_x=nx,angle=45
    geom_text(data=stats %>% filter(Contrast_type=="Interaction_Supplement"),aes(label=pStars,y=Inf),
              show.legend = FALSE,size=5,color="green4", vjust=vj)+#,nudge_x=nx,hjust=hj,angle=45
    ylab('Mean acs-2::gfp fluorescence per worm, a.u.')+
    xlab('Strain & Supplement')+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 9))+
    theme(legend.position="top")#+ #,axis.text.x = element_text(angle = 45, hjust = 1)
    #facet_wrap(~Gene,scale = "free",ncol=3)
}

#quartz()
f4plots<-f4 %>%
  group_by(Figure) %>%
  do(plot=plot_boxstats(.,sf4))
  

map2(paste0(odir,"/",as.character(f4plots$Figure),"_Fluorescence_log2.pdf"),
     f4plots$plot,
     width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial",ggsave)




#5 crp oe
f6strains<-c("OP50-C","crp, Pcrp (oe)")

f6oe<-data.frame(File=list.files("Data/crp (oe)/")) %>%
  group_by(File) %>%
  do(read_delim(paste0('Data/crp (oe)/',.$File),delim="\t") %>%
       gather(Strain,Abs,everything())) %>%
  ungroup %>%
  filter(!is.na(Abs)) %>%
  mutate(Strain=str_trim(Strain),
         Strain=ifelse(Strain=="crp (oe)","crp, Pcrp (oe)",Strain),
         Strain=factor(Strain,levels=f6strains,labels=f6strains),
         Replicate=ifelse(str_detect(File,"rep1"),1,2),
         Gene="acs-2") %>%
  mutate(Log=log2(Abs)) %>%
  mutate(Ref=mean(Log[Strain=='OP50-C']) ) %>%
  group_by(Replicate) %>%
  mutate(Norm=Log-mean(Log[Strain=='OP50-C'])+Ref,
         NormAbs=2^Norm) %>%
  ungroup
  

sf6<-f6oe %>%
  filter(!is.na(Norm))%>%
  do(tidy(lm(Norm~Strain,data=.))) %>%
  filter(term!="(Intercept)") %>%
  mutate(pStars=pStars(p.value),
         Strain=str_replace_all(term,"Strain",""),
         Contrast_type="Strain",
         term=NULL) %>%
  rename(SE=std.error,
         logFC=estimate) %>%
  select(Contrast_type,Strain,everything())


sf6 %>%
  write_csv(paste0(odir,'/6h_crp_oe_Summary.csv'))



f6oe %>%
  ggplot(aes(x=Strain,y=NormAbs,color=Strain))+
  #geom_blank(data = blank_data_f2k,aes(x=Worm,y=Abs))+
  geom_point(aes(fill=Strain),size=1,alpha=0.5,position = position_jitterdodge(dodge.width=1,jitter.width=0.4))+
  #geom_boxplot(data=bdata,aes(ymax=Max,ymin=Min,lower=Mean-SD,upper=Mean+SD,middle=Mean),position="identity",stat = "identity",alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = position_dodge(width=1), alpha=0.5) +
  scale_y_continuous(trans = 'log2',
                     breaks = scales::trans_breaks('log2', function(x) 2^x),
                     labels = scales::trans_format('log2', scales::math_format(2^.x))) + 
  geom_text(data=sf5 ,aes(label=pStars,y=2^12.5),show.legend = FALSE,size=5)+
  ylab('Mean acs-2::gfp fluorescence per worm, a.u.')+
  theme(legend.position="top")

ggsave(file=paste0(odir,'/6h_acs-2_crp_oe_logScale2.pdf'),
       width=50,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



f6oe %>%
  ggplot(aes(x="Strain",y=NormAbs,color=Strain))+
  #geom_blank(data = blank_data_f2k,aes(x=Worm,y=Abs))+
  geom_point(aes(fill=Strain),size=1,alpha=0.5,position = position_jitter(width=0.4))+
  #geom_boxplot(data=bdata,aes(ymax=Max,ymin=Min,lower=Mean-SD,upper=Mean+SD,middle=Mean),position="identity",stat = "identity",alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot", alpha=0.5) +
  scale_y_continuous(trans = 'log2',
                     breaks = scales::trans_breaks('log2', function(x) 2^x),
                     labels = scales::trans_format('log2', scales::math_format(2^.x))) + 
  geom_text(data=sf5 ,aes(label=pStars,y=2^12.5),show.legend = FALSE,size=5)+
  ylab('Mean acs-2::gfp fluorescence per worm, a.u.')+
  theme(legend.position="top",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(file=paste0(odir,'/6h_acs-2_crp_oe_logScale2_tiny.pdf'),
       width=25,height=44,units='mm',scale=2,device=cairo_pdf,family="Arial")



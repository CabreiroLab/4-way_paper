#Figure numbering might have been changed.
library(tidyverse)
library(scales)
library(broom)
library(forcats)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)


setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy/")

odir<-'Summary_vha6'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#load("Fluorescence_vha6.RData")
#save.image('Fluorescence_vha6.RData')

lmtests<-function(data,form) {
  
  val<-form %>% str_split("~") %>% unlist(.) %>% first()
  facs<-form %>% str_split("~") %>% unlist(.) %>% nth(2) %>% str_split("\\*") %>% unlist
  fac1<-facs[1]
  fac2<-facs[2]
  
  facregex<-paste(fac1,fac2,sep="|")
  
  groupings<-group_vars(data)
  
  #Effect of first factor
  stat1<-data %>%
    group_by_(.dots=c(groupings,fac2)) %>%
    do(tidy(lm(as.formula(paste(val,fac1,sep="~")),data=.)) ) %>%
    ungroup %>%
    filter(term!='(Intercept)') %>%
    mutate(term=str_replace_all(term,facregex,""),
           Contrast_type=fac1 )%>%
    rename_(.dots=setNames("term",fac1))
  
  
  #Effect of second factor
  stat2<-data %>%
    group_by_(.dots=c(groupings,fac1)) %>%
    do(tidy(lm(as.formula(paste(val,fac2,sep="~")),data=.)) ) %>%
    ungroup %>%
    filter(term!='(Intercept)') %>%
    mutate(term=str_replace_all(term,facregex,""),
           Contrast_type=fac2)%>%
    rename_(.dots=setNames("term",fac2))
  
  #Effect of interaction
  stat3<-data %>%
    do(tidy(lm(as.formula(form),data=.)) ) %>%
    ungroup %>%
    filter(str_detect(term,":")) %>%
    mutate(term=str_replace_all(term,facregex,""),
           Contrast_type="Interaction")%>%
    separate(term,c(fac1,fac2))
  
  stat<-stat1 %>%
    rbind(stat2) %>%
    rbind(stat3) %>%
    rename(SE=std.error,
           logFC=estimate) %>%
    mutate(PE=logFC+SE,
           NE=logFC-SE,
           Prc=2^logFC*100,
           PrcNE=2^NE*100,
           PrcPE=2^PE*100,
           pStars=pStars(p.value)) %>%
    mutate_at(c(fac1,fac2,"Contrast_type"),as.factor) %>%
    select(one_of(groupings),"Contrast_type",fac1,fac2,everything())
  
  return(stat)
}


#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau


data.crp<-read_delim('./Data/vha-6/pvha-6 c v crp only.txt',delim='\t') %>%
  gather(Group,Abs,c0:crp50) %>%
  filter(!is.na(Abs)) %>%
  mutate(Metformin_mM=ifelse(str_detect(Group,"50"),50,0) %>% factor(levels=c(0,50),labels=c("0","50")),
         Strain=ifelse(str_detect(Group,"crp"),"crp","OP50-C") %>% factor(levels=c("OP50-C","crp")),
         Log=log2(Abs)) %>%
  gather(Measure,Value,Abs,Log)


data.glu<-read_delim('./Data/vha-6/pvha-6 c v glucose only.txt',delim='\t') %>%
  gather(Group,Abs,c0:glu50) %>%
  filter(!is.na(Abs)) %>%
  mutate(Metformin_mM=ifelse(str_detect(Group,"50"),50,0) %>% factor(levels=c(0,50),labels=c("0","50")),
         Glucose=ifelse(str_detect(Group,"glu"),"0.2","0") %>% factor(levels=c("0","0.2")),
         Log=log2(Abs)) %>%
  gather(Measure,Value,Abs,Log)



form<-"Value~Metformin_mM*Strain"
stats.crp<-data.crp %>%
  group_by(Measure) %>%
  do(lmtests(data=.,form))


stats.crp %>%
  write_csv(paste0(odir,"/Stats_Summary_vha-6_crp.csv"))



form<-"Value~Metformin_mM*Glucose"
stats.glu<-data.glu %>%
  group_by(Measure) %>%
  do(lmtests(data=.,form))


stats.glu %>%
  write_csv(paste0(odir,"/Stats_Summary_vha-6_Glucose.csv"))




Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- c("0","50")
Metlab<-'Metformin, mM'


blank_data<-data.crp %>%
  filter(Measure=='Log' ) %>%
  group_by(Strain) %>%
  summarise(Value=2^(min(Value)+(max(Value)-min(Value))*1.3 ) ) %>% #ymin+2^( log2( max(NormAbs)/ymin )*1.5 )  #ymin+(max(NormAbs)-ymin)*1.5
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0","50")))


showstats<-stats.crp %>%
  filter(Measure=="Log" )

hj<-0.8
vj<-2
nx<--0.5

quartz()
data.crp %>%
  filter(Measure=='Abs' ) %>%
  ggplot+
  aes(x=Strain,y=Value,color=Metformin_mM)+
  geom_jitter(width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_blank(data = blank_data,aes(x=Strain,y=Value))+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x))) + 
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Strain')+
  scale_colour_manual(name = "Metformin,\nmM",values =Metcols)+
  geom_text(data=filter(showstats,Contrast_type=="Interaction"),aes(label=pStars,y=Inf,x=Strain),show.legend = FALSE,size=5,color="green4",nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=filter(showstats,Contrast_type=="Metformin_mM"),aes(label=pStars,y=Inf,x=Strain),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)

ggsave(file=paste0(odir,'/Fluorescence_crp_logScale2.pdf'),
       width=50,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



blank_data<-data.glu %>%
  filter(Measure=='Log' ) %>%
  group_by(Glucose) %>%
  summarise(Value=2^(min(Value)+(max(Value)-min(Value))*1.3 ) ) %>% #ymin+2^( log2( max(NormAbs)/ymin )*1.5 )  #ymin+(max(NormAbs)-ymin)*1.5
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0","50")))


showstats<-stats.glu %>%
  filter(Measure=="Log" )


hj<-0.8
vj<-2
nx<--0.5



quartz()
data.glu %>%
  filter(Measure=='Abs' ) %>%
  ggplot+
  aes(x=Glucose,y=Value,color=Metformin_mM)+
  geom_jitter(width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_blank(data = blank_data,aes(x=Glucose,y=Value))+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x))) + 
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Glucose, %')+
  scale_colour_manual(name = "Metformin,\nmM",values =Metcols)+
  geom_text(data=filter(showstats,Contrast_type=="Interaction"),aes(label=pStars,y=Inf,x=Glucose),show.legend = FALSE,size=5,color="green4",nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=filter(showstats,Contrast_type=="Metformin_mM"),aes(label=pStars,y=Inf,x=Glucose),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)


filter(showstats,Contrast_type=="Metformin_mM")

ggsave(file=paste0(odir,'/Fluorescence_Glucose_logScale2.pdf'),
       width=50,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")

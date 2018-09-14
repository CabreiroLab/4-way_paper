library(tidyverse)
library(scales)
library(broom)
library(forcats)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)



setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy/")

odir<-'Summary_Acetoacetate'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")




#load("Fluorescence_Acetoacetate.RData")
#save.image('Fluorescence_Acetoacetate.RData')



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




data.trans<-data.frame(Gene=c("acs-2","atgl-1","cpt-2","cpt-5")) %>%
  mutate(File=paste0(Gene,'.txt')) %>%
  group_by(Gene, File) %>%
  do(read_delim(paste('./Data/acetoacetate',.$File,sep='/'),delim='\t')) %>%
  gather(Group,Abs,`0`:`50_1`) %>%
  filter(!is.na(Abs)) %>%
  mutate(Metformin_mM=ifelse(str_detect(Group,"50"),50,0) %>% factor(levels=c(0,50),labels=c("0","50")),
         Acetoacetate_mM=ifelse(str_detect(Group,"_1"),10,0) %>% factor(levels=c(0,10),labels=c("0","10")),
         Log=log2(Abs)) %>%
  gather(Measure,Value,Abs,Log)




Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- c("0","50")
Metlab<-'Metformin, mM'




form<-"Value~Metformin_mM*Acetoacetate_mM"
stats<-data.trans %>%
  group_by(Measure,Gene) %>%
  do(lmtests(data=.,form))


stats %>%
  write_csv(paste0(odir,"/Stats_Summary_transgenes.csv"))




blank_data<-data.trans %>%
  filter(Measure=='Log' ) %>%
  group_by(Gene,Acetoacetate_mM) %>%
  summarise(Value=2^(min(Value)+(max(Value)-min(Value))*1.3 ) ) %>% #ymin+2^( log2( max(NormAbs)/ymin )*1.5 )  #ymin+(max(NormAbs)-ymin)*1.5
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0","50")))


showstats<-stats %>%
  filter(Measure=="Log" )


glimpse(data.all)
glimpse(stats)

showstats$Contrast_type


hj<-0.8
vj<-2
nx<--0.5

quartz()
data.trans %>%
  filter(Measure=='Abs' ) %>%
  ggplot+
  aes(x=Acetoacetate_mM,y=Value,color=Metformin_mM)+
  geom_jitter(width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_blank(data = blank_data,aes(x=Acetoacetate_mM,y=Value))+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x))) + 
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Acetoacetate, mM')+
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  geom_text(data=filter(showstats,Contrast_type=="Interaction"),aes(label=pStars,y=Inf,x=Acetoacetate_mM),show.legend = FALSE,size=5,color="green4",nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=filter(showstats,Contrast_type=="Metformin_mM"),aes(label=pStars,y=Inf,x=Acetoacetate_mM),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  theme(legend.position="top")+
  facet_wrap(~Gene,scale = "free_y",ncol=7)

ggsave(file=paste0(odir,'/Fluorescence_logScale2.pdf'),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")





data.tit<-read_delim('./Data/acetoacetate/acs-2 titration.txt',delim='\t')%>%
  gather(Group,Abs,`0`:`20 + Met`) %>%
  filter(!is.na(Abs)) %>%
  mutate(Metformin_mM=ifelse(str_detect(Group,"Met"),50,0) %>% factor(levels=c(0,50),labels=c("0","50")),
         Acetoacetate_mM=str_remove_all(Group,fixed(" + Met") )%>% str_trim %>% factor(levels=c("0","1","5","10","20"),labels=c("0","1","5","10","20")) ,
         Log=log2(Abs)) %>%
  gather(Measure,Value,Abs,Log)


data.tit %>%
  View

form<-"Value~Metformin_mM*Acetoacetate_mM"
stats.tit<-data.tit%>%
  group_by(Measure) %>%
  do(lmtests(data=.,form))


stats.tit %>%
  write_csv(paste0(odir,"/Stats_Summary_titration.csv"))



blank_data<-data.tit %>%
  filter(Measure=='Log' ) %>%
  group_by(Acetoacetate_mM) %>%
  summarise(Value=2^(min(Value)+(max(Value)-min(Value))*1.3 ) ) %>% #ymin+2^( log2( max(NormAbs)/ymin )*1.5 )  #ymin+(max(NormAbs)-ymin)*1.5
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0","50")))


showstats<-stats.tit %>%
  filter(Measure=="Log" )

quartz()
data.tit %>%
  filter(Measure=='Abs' ) %>%
  ggplot+
  aes(x=Acetoacetate_mM,y=Value,color=Metformin_mM)+
  
  geom_jitter(width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_blank(data = blank_data,aes(x=Acetoacetate_mM,y=Value))+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x))) + 
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Acetoacetate, mM')+
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  geom_text(data=filter(showstats,Contrast_type=="Interaction"),aes(label=pStars,y=Inf,x=Acetoacetate_mM),show.legend = FALSE,size=5,color="green4",nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=filter(showstats,Contrast_type=="Metformin_mM"),aes(label=pStars,y=Inf,x=Acetoacetate_mM),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  theme(legend.position="top")

ggsave(file=paste0(odir,'/Fluorescence_titration_logScale2.pdf'),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")




data.tit.sum <- data.tit %>%
  filter(Measure=='Log') %>%
  group_by(Measure,Metformin_mM,Acetoacetate_mM) %>%
  summarise(Mean=mean(Value),
            SD=sd(Value),
            SE=SD/n()) %>%
  mutate(Abs=2^Mean,
         APE=2^(Mean+SE),
         ANE=2^(Mean-SE),
         APD=2^(Mean+SD),
         AND=2^(Mean-SD))



quartz()
data.tit.sum %>%
  ggplot+
  aes(x=Acetoacetate_mM,y=Abs,color=Metformin_mM,fill=Metformin_mM)+
  geom_ribbon(aes(ymin=APD,ymax=AND,group=Metformin_mM),color=NA,alpha=0.5)+
  geom_line(aes(group=Metformin_mM))+
  geom_point()+
  geom_blank(data = blank_data,aes(x=Acetoacetate_mM,y=Value))+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x))) + 
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Acetoacetate, mM')+
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  scale_fill_manual(name = "Metformin, mM",values =Metcols)+
  # geom_text(data=filter(showstats,Contrast_type=="Interaction"),aes(label=pStars,y=Inf,x=Acetoacetate_mM),show.legend = FALSE,size=5,color="green4",nudge_x=nx, vjust=vj,angle=45)+
  # geom_text(data=filter(showstats,Contrast_type=="Metformin_mM"),aes(label=pStars,y=Inf,x=Acetoacetate_mM),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  geom_text(data=filter(showstats,Contrast_type=="Interaction"),aes(label=pStars,y=2^32.25,x=Acetoacetate_mM),show.legend = FALSE,size=5,color="green4")+
  geom_text(data=filter(showstats,Contrast_type=="Metformin_mM"),aes(label=pStars,y=2^31.75,x=Acetoacetate_mM),show.legend = FALSE,size=5)+
  theme(legend.position="top")

ggsave(file=paste0(odir,'/Fluorescence_titration_logScale2_line.pdf'),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")

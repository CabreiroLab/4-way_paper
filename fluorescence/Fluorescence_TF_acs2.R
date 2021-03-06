#Figure numbering might have been changed
library(tidyverse)
library(PFun)
library(broom)

#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau

Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- c("0","50")
Metlab<-'Metformin, mM'

setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy")

odir<-'Summary_TF_acs-2'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

#load('Metformin_TF_acs-2.RData')
#save.image('Metformin_TF_acs-2.RData')


translations<-c('gvcA'='gcvA','op50-c'='OP50-C','op50'='OP50-C')

allfiles<-read_csv('TF_acs-2_segmentation/Allfiles.csv') %>%
  rbind(read_csv('TF_acs-2_segmentation/Allfiles_additional.csv')) %>%
  filter( !str_detect(Folder,'grouped') ) %>%
  mutate(Dataset=ifelse(str_detect(Replicate_folder,'additional'),"Additional","Original"),
         Type=ifelse(str_detect(Folder,'wo metf'),'C','T' ),
         Replicate=as.integer(ifelse(str_detect(Replicate_folder,'Rep1'),1,2)),
         Gene=str_trim(str_replace_all(Folder,'metf|wo metf|later|met repeat','')),
         Gene=ifelse(Gene %in% names(translations),translations[Gene],Gene ),
         FileNo=str_extract(File,'[[:digit:]]{1,4}.tiff') %>% str_replace('.tiff','') %>% str_pad(width=4,pad='0') ) %>%
  group_by(Dataset,Replicate,Gene,Type,Folder) %>%
  arrange(FileNo) %>%
  mutate(FileInd=row_number()) %>%
  ungroup %>%
  #filter(!FileNo=='0001') %>%
  group_by(Dataset,Replicate,Gene,Type,Folder) %>%
  mutate(Count=n()) %>%
  ungroup %>%
  group_by(Gene) %>%
  mutate(Max=as.integer(max(Count))) %>%
  select(Dataset,Replicate,Type,Gene,FileNo,FileInd,Replicate_folder:File,Max) %>%
  ungroup


allfiles %>%
  group_by(Dataset,Replicate,Gene,Type) %>%
  summarise(Count=n()) %>%
  ungroup %>%
  unite(RT,Replicate,Type,remove = FALSE) %>%
  select(Gene,RT,Count,Dataset) %>%
  spread(RT,Count)
  
allfiles_sum<-allfiles %>%
  group_by(Dataset,Replicate,Gene,Type,Folder) %>%
  summarise(Count=n()) %>%
  ungroup %>%
  group_by(Gene) %>%
  summarise(Max=as.integer(max(Count)))

write_csv(allfiles,paste0(odir,"/Allfiles_annotated.csv"))
write_csv(allfiles_sum,paste0(odir,'/Allfiles_summary.csv'))


strfix<-c("op50-c"="OP50-C")
genes<-c('OP50-C','crp','cra','argR','ntrC','arcA','csiR','fur','gcvA',"marA",'mlc',"nac")
remindex<-c('1_crp_T_8')


cresults<-data.frame(CFile=list.files(path = paste0("TF_acs-2_segmentation/TF_acs-2_control_groups/"))) %>%
  group_by(CFile) %>%
  do(read_delim(paste('./TF_acs-2_segmentation/TF_acs-2_control_groups',.$CFile,sep='/'),delim='\t') %>% gather(Gene,W_Int,everything()) %>% mutate_all(as.character) ) %>%
  ungroup %>%
  mutate(Gene=str_trim(Gene),
         Gene=ifelse(Gene %in% names(strfix),strfix[Gene],Gene),
         Dataset=ifelse(CFile %in% c("TF acs-2 controls rep 1 23-11-17.txt",
                                     "TF acs-2 controls rep 2 23-11-17.txt"),
                        "Additional","Original"),
         Replicate=ifelse(str_detect(CFile,"rep 1"),1,2),
         W_Int=as.numeric(W_Int),
         Worm=1,
         Type="C",
         B_N=NA,B_Sum=NA,B_Mean=11.5533,
         W_N=NA,W_Sum=NA,W_Mean=W_Int+B_Mean) %>%
  filter(!is.na(W_Int))
  


#Collect results
results.o<-read_csv('TF_acs-2_segmentation/All_results_onlyprc.csv') %>%
  mutate(X1=NULL) %>%
  mutate(FileName=File,
         File=paste('.',Replicate_folder,Folder,FileName,sep='/')) %>%
  select(-c(FileName,Replicate_folder,Folder,FileName,Replicate,Gene,Type))

results.a<-read_csv('TF_acs-2_segmentation/All_results_onlyprc_additional.csv') %>%
  mutate(X1=NULL)


results.all<-rbind(results.o,results.a)


data<-allfiles %>%
  mutate(FileName=File,
         File=paste('.',Replicate_folder,Folder,FileName,sep='/')) %>%
  left_join(results.all) %>%
  unite(Index,Replicate,Gene,Type,FileInd,remove = FALSE) %>%
  filter(!Index %in% remindex) %>%
  filter(W_N<100000,
         W_N>40000,
         FileNo!='0001',
         !(Replicate==1 & FileInd==2),
         Folder!='OP50-C metf later' & Type=="T") %>%
  mutate(W_Int=W_Mean-B_Mean) %>%
  select(-c(FileNo,FileInd,Replicate_folder,Folder,File,Max,FileName,Index)) %>%
  rbind(cresults %>% select(-CFile)) %>%
  mutate(Metformin_mM=ifelse(Type=="C",0,50)) %>%
  mutate_at(c('Replicate','Type','Gene','Worm',"Metformin_mM"),as.factor) %>%
  mutate(Gene=factor(Gene,levels=genes)) %>%
  group_by(Dataset,Replicate) %>%
  mutate(LogInt=log2(W_Int),
         TpRep_Mean=mean(LogInt[Gene=='OP50-C'])) %>%
  ungroup %>%
  mutate(NormLog=mean(TpRep_Mean)+LogInt-TpRep_Mean,
         NormAbs=2^NormLog) %>%
  ungroup 

data %>%
  write_csv(paste0(odir,"/All_data_raw.csv"))


form<-"NormLog~Metformin_mM*Gene"

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


stats<-data %>%
  ungroup %>%
  do(lmtests(data=.,"NormLog~Metformin_mM*Gene"))


stats %>%
  write_csv(paste0(odir,"/Stats_Summary.csv"))


showstats<-stats %>%
  filter(Contrast_type %in% c("Gene","Interaction"))

CT<-data %>%
  filter(Gene=="OP50-C") %>%
  group_by(Metformin_mM) %>%
  summarise(Ref=mean(NormAbs))


blank_data <- data %>%
  group_by(Gene) %>%
  summarise(NormAbs=2^(min(NormLog)+(max(NormLog)-min(NormLog))*1.3 ) )%>% 
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0",'50')))



hj<-0.8
vj<-2
nx<--0.25

#quartz()
data %>%
  ggplot+
  aes(x=Gene,y=NormAbs,color=Metformin_mM)+
  geom_hline(data=CT,aes(yintercept = Ref,color=Metformin_mM))+
  geom_jitter(aes(fill=Metformin_mM),width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +#
  geom_blank(data = blank_data,aes(x=Gene,y=NormAbs))+
  scale_y_continuous(trans = 'log2',
                     breaks = scales::trans_breaks('log2', function(x) 2^x),
                     labels = scales::trans_format('log2', scales::math_format(2^.x))) + 
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  scale_fill_manual(name = "Metformin, mM",values =Metcols)+
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Gene')+
  geom_text(data=showstats %>% filter(Contrast_type=="Gene" & Metformin_mM=="50" ) ,aes(label=pStars,y=2^8.3),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=showstats %>% filter(Contrast_type=="Gene"& Metformin_mM=="0") ,aes(label=pStars,y=2^3),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=showstats %>% filter(Contrast_type=="Interaction"),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,color="green4",nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  # geom_text(data=showstats %>% filter(Contrast_type=="Metformin_mM") ,aes(label=pStars,y=Inf),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj,angle=45)+
  # geom_text(data=showstats %>% filter(Contrast_type=="Interaction"),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,color="red",nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  theme(legend.position="top",axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file=paste(odir,'/Fluorescence_logScale2.pdf',sep = ''),
       width=85,height=40,units='mm',scale=2,device=cairo_pdf,family="Arial")


#Treatment figure:
Metcols2 <- c("#32006F","#FF0000")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols2) <- c("50","0")
Metlab2<-'Metformin,\nmM'


quartz()
data %>%
  filter(Gene %in% c('OP50-C','crp','cra','argR','ntrC')) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=c(50,0))) %>%
  ggplot+
  aes(x=Gene,y=NormAbs,color=Metformin_mM)+
  #geom_hline(data=CT,aes(yintercept = Ref,color=Metformin_mM))+
  geom_jitter(aes(fill=Metformin_mM),width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +#
  #geom_blank(data = blank_data,aes(x=Gene,y=NormAbs))+
  scale_y_continuous(trans = 'log2',
                     breaks = scales::trans_breaks('log2', function(x) 2^x),
                     labels = scales::trans_format('log2', scales::math_format(2^.x))) + 
  scale_colour_manual(name = Metlab2,values =Metcols2)+
  scale_fill_manual(name = Metlab2,values =Metcols2)+
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Gene')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_blank())+
  facet_grid(Metformin_mM~.,scales = "free_y")


ggsave(file=paste(odir,'/Fluorescence_logScale2_axisbreak_Thesis.pdf',sep = ''),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



#quartz()
data %>%
  ggplot+
  aes(x=Gene,y=NormAbs,color=Metformin_mM)+
  geom_hline(data=CT,aes(yintercept = Ref,color=Metformin_mM))+
  geom_jitter(aes(fill=Metformin_mM),width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +#
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  scale_fill_manual(name = "Metformin, mM",values =Metcols)+
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Gene')+
  geom_text(data=showstats %>% filter(Contrast_type=="Gene" & Metformin_mM=="50" ) ,aes(label=pStars,y=130),show.legend = FALSE,size=5)+#,nudge_x=nx, vjust=vj,angle=45
  #geom_text(data=showstats %>% filter(Contrast_type=="Gene"& Metformin_mM=="0") ,aes(label=pStars,y=120),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=showstats %>% filter(Contrast_type=="Interaction"),aes(label=pStars,y=140),show.legend = FALSE,size=5,color="green4")+#,nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj
  # geom_text(data=showstats %>% filter(Contrast_type=="Metformin_mM") ,aes(label=pStars,y=Inf),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj,angle=45)+
  # geom_text(data=showstats %>% filter(Contrast_type=="Interaction"),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,color="red",nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  theme(legend.position="top",axis.text.x = element_text(angle = 45, hjust = 1))#

ggsave(file=paste(odir,'/Fluorescence_AbsScale_55.pdf',sep = ''),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



amp<-2

minv<- -amp
maxv<- amp

nstep<-maxv-minv
nstep<-10

clrbrks<-seq(-amp,amp,by=1)
#patclrscale <- colorRampPalette(c("purple", "gray50","green"))(n = nstep)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)

comps<-c("Gene 0"="0mM Metf","Gene 50"="50mM Metf","Interaction 50"="Difference")

stats %>%
  filter(Contrast_type %in% c("Gene","Interaction")) %>%
  mutate(Comparison=paste(Contrast_type,Metformin_mM),
         Label=comps[Comparison],
         Gene=factor(Gene,levels=rev(unique(as.character(Gene))) )  ,
         Label=factor(Label,levels=c("0mM Metf","50mM Metf","Difference")) ) %>%
  ggplot(aes(x=Label,y=Gene))+
  geom_tile(aes(fill=logFC))+
  ylab("Strain")+
  xlab("Comparison")+
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  geom_text(aes(label=as.character(pStars)))+
  theme_Heatmap()#+theme(axis.text.x = element_text(angle=45))
  

ggsave(file=paste(odir,'/Heatmap_Fluorescence_vertical.pdf',sep = ''),
       width=34,height=63,units='mm',scale=2,device=cairo_pdf,family="Arial")



stats %>%
  filter(Contrast_type %in% c("Gene","Interaction")) %>%
  mutate(Comparison=paste(Contrast_type,Metformin_mM),
         Label=comps[Comparison],
         Gene=factor(Gene,levels=unique(as.character(Gene)) )  ,
         Label=factor(Label,levels=c("Difference","50mM Metf","0mM Metf"),labels=c("Interaction","50mM Metf","0mM Metf")) ) %>%
  filter(Gene %in% c('crp','cra','argR','ntrC')) %>%
  ggplot(aes(x=Gene,y=Label))+
  geom_tile(aes(fill=logFC))+
  ylab("Comparison")+
  xlab("DIfference")+
  labs(fill='Pacs-2::gfp\nlogFC vs\nOP50-C')+
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  geom_text(aes(label=as.character(pStars)))+
  theme_Heatmap()+
  theme(axis.text.x = element_text(angle=45))


ggsave(file=paste(odir,'/Heatmap_Fluorescence_horizontal.pdf',sep = ''),
       width=55,height=31,units='mm',scale=2,device=cairo_pdf,family="Arial")


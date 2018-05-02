library(tidyverse)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)
library(scales)


setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy/")


odir<-'Summary_Transgenes'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#load("Fluorescence_Transgenes.RData")
#save.image('Fluorescence_Transgenes.RData')


#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau



strains<-c('OP50-C','OP50-MR','crp','cra')
Sstrains<-c('OP50-C','OP50-C-Glu','OP50-MR','crp','cra')


data <- read_csv('All_raw_data.csv') %>%
  filter(Type %in% c('CRP','RNAseq')) %>%
  mutate_at(c('Type','Gene','Strain','Metformin_mM','Supplement','SGroup','Condition','ID','Replicate'),as.factor) %>%
  mutate(Strain=factor(Strain,levels=strains),
         SGroup=factor(SGroup,levels=Sstrains))



data %>%
  filter(Measure=='Log' & Strain=="OP50-C" & Gene=='dhs-23' & Metformin_mM==50 & NormAbs<2^6)


unique(data$SGroup)

unique(data$Strain)


Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- levels(data$Metformin_mM)
Metlab<-'Metformin,\nmM'




# data %>%
#   filter(Measure=='Log') %>%
#   group_by(Type,Gene,Replicate,Strain,Supplement,Metformin_mM) %>%
#   summarise(Count=n()) %>%
#   View()
  

unique(data$Gene)
  

# data %>%
#   group_by(Condition) %>%
#   summarise %>%
#   write_csv('Conditions_raw_transgenes.csv')






plotBox<-function(data,yval,ylb){
  data %>%
  ggplot+
    aes_string(x="SGroup",y=yval,color="Metformin_mM")+
    stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
    geom_jitter(width=0.25)+
    ylab(ylb)+
    xlab('Bacterial strain + Supplement')+
    scale_colour_manual(name = Metlab,values =Metcols)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_grid(Gene~Type+Replicate,scale="free_y")
}


plotBoxC<-function(data,yval,ylb){
  data %>%
    ggplot+
    aes_string(x="SGroup",y=yval,color="Metformin_mM")+
    stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
    geom_jitter(width=0.25)+
    ylab(ylb)+
    xlab('Bacterial strain + Supplement')+
    scale_colour_manual(name = Metlab,values =Metcols)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_grid(Gene~.,scale="free_y")
}




Fluorplots<-data %>%
  group_by(Measure) %>%
  do(plot=plotBox(.,"Value",as.character(unique(.$Measure))) )


map2(paste0(odir,"/Fluorescence_raw_",as.character(Fluorplots$Measure),".pdf"),
     Fluorplots$plot,
     width=10,height=30, useDingbats=FALSE, ggsave)


Fluorplots<-data %>%
  group_by(Measure) %>%
  do(plot=plotBox(.,"Norm",as.character(unique(.$Measure))) )


map2(paste0(odir,"/Fluorescence_norm_",as.character(Fluorplots$Measure),".pdf"),
     Fluorplots$plot,
     width=10,height=30, useDingbats=FALSE, ggsave)
  

Fluorplots<-data %>%
  group_by(Measure) %>%
  do(plot=plotBoxC(.,"Norm",as.character(unique(.$Measure))) )


map2(paste0(odir,"/Fluorescence_norm_joined_",as.character(Fluorplots$Measure),".pdf"),
     Fluorplots$plot,
     width=6,height=30, useDingbats=FALSE, ggsave)






#
sum.c<-data %>%
  filter(Measure=='Log') %>%
  group_by(Gene,ID,Strain,Supplement,Metformin_mM,SGroup) %>%
  summarise(SD=sd(Value,na.rm = TRUE),Mean=mean(Value,na.rm = TRUE)) %>%
  mutate(Index=paste(Gene,ID),
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
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=7,height=20, useDingbats=FALSE)



#Linear modelling tests
#
#
# model<-lm(LogValue~0+ID, data=data)
# lmod_glht <- multcomp::glht(model, linfct = contr.matrix)
# result<-multcomp:::summary.glht(lmod_glht,test=multcomp::adjusted("none"))
# res<-data.frame(result$test[c('coefficients','sigma','tstat','pvalues')])


allgroups<-as.character(unique(data$ID))
allgroups

contrasts<-read.contrasts('!Contrasts_fluorescence.xlsx')



contrasts$Contrasts.table
contrasts.desc<-contrasts$Contrasts.table%>%
  select(Description:Metformin_mM)


contr.matrix<-contrasts$Contrasts.matrix
contr.matrix




results.all<-data %>%
  filter(Measure=='Log') %>%
  group_by(Gene) %>%
  do(hypothesise(.,"Norm~0+ID",contr.matrix)) %>%
  getresults(contrasts.desc)



results<-results.all$results
results.castfull<-results.all$castfull
results.cast<-results.all$cast
results.multi<-results.all$multi

head(results.castfull)
head(results.cast)

View(results)


View(results.new)



write.csv(results,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)







#Generate table for heatmap



heatsum<-results %>%
  filter(Contrast_type %in% c('Treatment','Interaction' ) ) %>%
  select(Description,Gene,logFC) %>%
  spread(Description,logFC) %>%
  data.frame(check.names = FALSE,check.rows = FALSE)


rownames(heatsum)<-heatsum$Gene
heatsum$Gene<-NULL



max(heatsum)
min(heatsum)

amp<-5

minv<- -amp
maxv<- amp

nstep<-maxv-minv
nstep<-8

clrbrks<-seq(-amp,amp,by=2)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)

d<-dist(as.matrix(heatsum),method = "euclidean")
h<-hclust(d)
ordmet<-rownames(heatsum[h$order,])

if (length(ordmet)!=length(unique(ordmet))){
  print("Non unique metabolites!")
}





results %>%
  filter(Contrast %in% c('OP50_T','OP50-MR_T','OP50-MR_I') & Gene %in% MRgene) %>%
  mutate(Gene=factor(Gene,levels=rev(MRgene) ))%>%
  ggplot(aes(x=Description,y=Gene) )+
  geom_tile(aes(fill=logFC))+
  geom_text(aes(label=as.character(pStars)),size=5)+
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  xlab("Comparison")+
  theme_Heatmap()

ggsave(file=paste(odir,'/Comparison_Heatmap_Fluorescence.pdf',sep = ''),
       width=40,height=70,units='mm',scale=2,device=cairo_pdf,family="Arial")




descriptions<-c('SM-S'='Treatment effect on OP50','RM-R'='Treatment effect on OP50-MR','SM-S-(RM-R)'='Interaction between treatment and OP50-MR')
translations<-c('F44G3.2'='argk-1','C05D11.7'='atgl-1')

genes<-unique(as.character(results$Gene))



# RNAseq.old<-read_csv('~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/Results/All_results.csv') %>%
#   rename(Contrast=Comparison,Gene=gene_name) %>%
#   filter(Gene %in% c(genes,'F44G3.2','C05D11.7') & Contrast %in% c('SM-S','RM-R','SM-S-(RM-R)')) %>%
#   select(Contrast,Gene,logFC,FDR) %>%
#   mutate(Description=descriptions[Contrast],
#          Type='RNAseq old',
#          Gene=ifelse(Gene %in% names(translations),translations[Gene],Gene),
#          Stars=pStars(FDR),
#          logFC=ifelse(Contrast=='SM-S-(RM-R)',-logFC,logFC))



RNAseq.new<-read_csv('~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/Results_1thrs_newannot/All_results.csv') %>%
  rename(Contrast=Comparison,Gene=external_gene_name) %>%
  filter(Gene %in% c(genes,'F44G3.2','C05D11.7') & Contrast %in% c('SM-S','RM-R','SM-S-(RM-R)')) %>%
  select(Contrast,Gene,logFC,FDR) %>%
  mutate(Description=descriptions[Contrast],
         Type='RNAseq',
         Gene=ifelse(Gene %in% names(translations),translations[Gene],Gene),
         Stars=pStars(FDR),
         logFC=ifelse(Contrast=='SM-S-(RM-R)',-logFC,logFC))


#& Gene !='atgl-1'

RNAresults<-results %>%
  filter(Contrast %in% c('OP50_T','OP50-MR_T','OP50-MR_I') ) %>%
  select(Gene, Contrast, Description,logFC,FDR,pStars) %>%
  rename(Stars=pStars) %>%
  mutate(Type='Fluorescence') %>%
  rbind(RNAseq.new)





amp<-10

minv<- -amp
maxv<- amp

nstep<-maxv-minv
nstep<-10

clrbrks<-seq(-amp,amp,by=4)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)

RNAgenes<-c('acs-2','F37H8.3', 'dhs-23', 'fat-7', 'cpt-2', 'cpt-5','atgl-1')

RNAresults %>%
  #filter(Gene %in% RNAgenes) %>%
  mutate(Description=recode(Description,
                            "Interaction between treatment and OP50-MR"="Difference",
                            "Treatment effect on OP50-MR"="OP50-MR + Metf",
                            "Treatment effect on OP50"="OP50-C + Metf"),
         Description=factor(Description,levels=rev(c("OP50-C + Metf","OP50-MR + Metf","Difference"))),
         Type=factor(Type,levels=c("RNAseq","Fluorescence")))%>%
  ggplot(aes(x=Gene,y=Description))+
  geom_tile(aes(fill=logFC))+
  geom_text(aes(label=as.character(Stars)),size=5,nudge_y = -0.15)+
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  xlab("Gene")+
  ylab("Comparison")+
  facet_grid(Type~.)+
  theme_Heatmap()+
  theme(axis.text.x = element_text(angle=45))

ggsave(file=paste0(odir,'/Comparison_Heatmap_RNAseq_Fluorescence_complete_horizontal.pdf'),
             width=83,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



#Additional fluorescence figure
nlpgene<-base::setdiff(unique(data$Gene),RNAgenes)



crpgene<-c('acs-2','atgl-1','cpt-2','cpt-5','dhs-23')
crpstrain<-c('OP50-C','crp','cra')


MRgene<-c('acs-2','atgl-1','cpt-2','cpt-5','dhs-23','F37H8.3','fat-7')
MRstrain<-c('OP50-C','OP50-MR')
MRcont<-c('OP50_T','OP50-MR_T','OP50-MR_I')



selgene<-MRgene
selstrain<-MRstrain
selcont<-MRcont


blank_data<-data %>%
  filter(Gene %in% selgene & Strain %in% selstrain & Measure=='Log' & Supplement=="None") %>%
  group_by(Gene,Strain,SGroup) %>%
  summarise(NormAbs=2^(min(Norm)+(max(Norm)-min(Norm))*1.3 ) )%>% #ymin+2^( log2( max(NormAbs)/ymin )*1.5 )  #ymin+(max(NormAbs)-ymin)*1.5
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0",'50')))

showstats<-results %>%
  filter(Contrast %in% selcont & Gene %in% selgene & Strain %in% selstrain)
  

hj<-0.8
vj<-2
nx<--0.5

quartz()
data %>%
  filter(Gene %in% selgene & Strain %in% selstrain & Measure=='Log' & Supplement=="None") %>%
  ggplot+
  aes(x=SGroup,y=NormAbs,color=Metformin_mM)+
  geom_jitter(width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_blank(data = blank_data,aes(x=Strain,y=NormAbs))+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x))) + 
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Strain')+
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  geom_text(data=showstats %>% filter(Contrast_type=="Interaction"),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,color="green4",nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=showstats %>% filter(Contrast_type=="Treatment"),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Gene,scale = "free_y",ncol=7)

ggsave(file=paste(odir,'/Fluorescence_logScale2_Fig2_horizontal_scaled_ncol7.pdf',sep = ''),
       width=100,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")

#facet_grid(.~Gene,scales = "free")
#scale="free_y",

ggsave(file=paste(odir,'/Fluorescence_logScale2_Fig2.pdf',sep = ''),
       width=55,height=123,units='mm',scale=2,device=cairo_pdf,family="Arial")

ggsave(file=paste(odir,'/Fluorescence_logScale2_Fig2_horizontal_scaled_ncol3.pdf',sep = ''),
       width=84,height=81,units='mm',scale=2,device=cairo_pdf,family="Arial")









selgene<-nlpgene

showstats<-results %>%
  filter(Contrast %in% selcont & Gene %in% selgene & Strain %in% selstrain)


blank_data<-data %>%
  filter(Gene %in% selgene & Strain %in% selstrain & Measure=='Log' & Supplement=="None") %>%
  group_by(Gene,Strain,SGroup) %>%
  summarise(NormAbs=2^(min(Norm)+(max(Norm)-min(Norm))*1.3 ) )%>% #ymin+2^( log2( max(NormAbs)/ymin )*1.5 )  #ymin+(max(NormAbs)-ymin)*1.5
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0",'50')))


quartz()

hj<-0.8
vj<-2
nx<--0.35

data %>%
  filter(Gene %in% selgene & Strain %in% selstrain & Measure=='Log' & Supplement=="None") %>%
  ggplot+
  aes(x=SGroup,y=NormAbs,color=Metformin_mM)+
  geom_jitter(width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_blank(data = blank_data,aes(x=Strain,y=NormAbs))+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x))) + 
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Strain')+
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  geom_text(data=showstats %>% filter(Contrast_type=="Interaction" ),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,color="green4",nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=showstats %>% filter(Contrast_type=="Treatment" ),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Gene,scale = "free_y",ncol=7)

ggsave(file=paste(odir,'/Fluorescence_logScale2_Fig2_nonlipid_horizontal_scaled_ncol5.pdf',sep = ''),
       width=80,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")







#Fig5 Transgenes

selgene<-crpgene
selstrain<-c("OP50-C","crp")
selcont<-c("OP50_T","CRP_T","CRP_I")

showstats<-results %>%
  filter(Contrast %in% selcont & Gene %in% selgene & Strain %in% selstrain)


blank_data<-data %>%
  filter(Gene %in% selgene & Strain %in% selstrain & Measure=='Log' & Supplement=="None") %>%
  group_by(Gene,Strain,SGroup) %>%
  summarise(NormAbs=2^(min(Norm)+(max(Norm)-min(Norm))*1.3 ) )%>% #ymin+2^( log2( max(NormAbs)/ymin )*1.5 )  #ymin+(max(NormAbs)-ymin)*1.5
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0",'50')))



quartz()

hj<-0.8
vj<-2
nx<--0.35

data %>%
  filter(Gene %in% crpgene & Strain %in% selstrain & Measure=='Log' & Supplement=="None") %>%
  ggplot+
  aes(x=SGroup,y=NormAbs,color=Metformin_mM)+
  geom_jitter(width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_blank(data = blank_data,aes(x=Strain,y=NormAbs))+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x))) + 
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Strain')+
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  geom_text(data=showstats %>% filter(Contrast_type=="Interaction" ),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,color="green4",nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=showstats %>% filter(Contrast_type=="Treatment" ),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Gene,scale = "free_y",ncol=5)

ggsave(file=paste(odir,'/Fluorescence_logScale2_Fig5_CRP_transgenes.pdf',sep = ''),
       width=62,height=40,units='mm',scale=2,device=cairo_pdf,family="Arial")





amp<-4

minv<- -amp
maxv<- amp

nstep<-maxv-minv
nstep<-10

clrbrks<-seq(-amp,amp,by=1)
#patclrscale <- colorRampPalette(c("purple", "gray50","green"))(n = nstep)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)


# showstats %>%
#   filter(Contrast_type %in% c("Treatment","Interaction")) %>%
#   pull(logFC) %>%
#   min



labeling<-c("OP50_T"="OP50-C + 50mM Metf","CRP_T"="crp + 50mM Metf","CRP_I"="Interaction")



showstats %>%
  filter(Contrast_type %in% c("Treatment","Interaction")) %>%
  mutate(Labels=labeling[as.character(Contrast)],
         Labels=factor(Labels,levels=rev(c("OP50-C + 50mM Metf","crp + 50mM Metf","Interaction")) ) ) %>%
  ggplot(aes(x=Gene,y=Labels))+
  geom_tile(aes(fill=logFC))+
  xlab("Strain")+
  ylab("Comparison")+
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  geom_text(aes(label=as.character(pStars)))+
  theme_Heatmap()+
  theme(axis.text.x = element_text(angle=45))

ggsave(file=paste(odir,'/Heatmap_Fluorescence_Fig5_CRP_transgenes_horizontal.pdf',sep = ''),
       width=55,height=25,units='mm',scale=2,device=cairo_pdf,family="Arial")




quartz()


selgene<-c("acs-2")
selstrain<-c("OP50-C")

showstats<-results %>%
  filter(Contrast %in% c("OP50_T","OP50_Glu","OP50Glu_I") & Gene %in% selgene & Strain %in% selstrain)



blank_data<-data %>%
  filter(Gene %in% selgene & Strain %in% selstrain & Measure=='Log') %>%
  group_by(Gene,Strain,SGroup) %>%
  summarise(NormAbs=2^(min(Norm)+(max(Norm)-min(Norm))*1.3 ) )%>% #ymin+2^( log2( max(NormAbs)/ymin )*1.5 )  #ymin+(max(NormAbs)-ymin)*1.5
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0",'50')))



quartz()
hj<-0.8
vj<-2
nx<--0.1

data %>%
  filter(Gene %in% c("acs-2") & Strain %in% c("OP50-C") & Measure=='Log') %>%
  ggplot+
  aes(x=SGroup,y=NormAbs,color=Metformin_mM)+
  geom_jitter(width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_blank(data = blank_data,aes(x=Strain,y=NormAbs))+
  scale_y_continuous(trans = 'log2',
                     breaks = trans_breaks('log2', function(x) 2^x),
                     labels = trans_format('log2', math_format(2^.x))) + 
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Strain & Supplement')+
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  geom_text(data=showstats %>% filter(Contrast_type=="Interaction" ),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,color="green4",nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=showstats %>% filter(Contrast_type=="Treatment" ),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  theme(legend.position = "none")

#+facet_wrap(~Gene,scale = "free_y",ncol=5)



ggsave(file=paste(odir,'/Fluorescence_logScale2_Fig6_OP50_Glucose_tiny.pdf',sep = ''),
       width=25,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")












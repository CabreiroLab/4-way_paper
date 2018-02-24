library(tidyverse)

devtools::install_github("PNorvaisas/PFun")
library(PFun)


setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy/")


odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#load("FluorescenceMicroscopy.RData")
save.image('FluorescenceMicroscopy.RData')


theme_set(theme_light())

library(ggthemes)

theme_Publication <- function(base_size=14) {
  
  (theme_foundation(base_size=base_size)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           #panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           #legend.position = "bottom",
           #legend.direction = "horizontal",
           #legend.key.size= unit(0.2, "cm"),
           #legend.margin = unit(0, "cm"),
           #legend.title = element_text(face="italic"),
           #plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}

theme_set(theme_Publication())


info<-read_csv('Conditions.csv')


#Translations for consistency
conditions<-c("glucose 0"="glu0","glucose 50"="glu50","r0"="mr0","r50"="mr50","C0"="c0","C50"="c50","R0"="mr0","R50"="mr50")


data<-data.frame(Type=c("RNAseq","CRP"),Folder=c("RNAseq","crp-cra-glucose")) %>%
  group_by(Folder,Type) %>%
  do(data.frame(File=list.files(path = paste0("./Data/",.$Folder)))) %>%
  mutate(Gene=str_replace_all(File,'_rep[[:digit:]]|.txt|.csv','')) %>%
  group_by(Folder,Type,File,Gene) %>%
  do(read_delim(paste('./Data',.$Folder,.$File,sep='/'),delim=ifelse(str_detect(.$File,fixed('.csv')),',','\t')) %>% gather(Condition,Abs,everything()) ) %>%
  group_by(Folder,Type,File,Gene,Condition) %>%
  mutate(Rank=row_number(),
         Replicate=as.integer(ifelse(Rank<31,1,2)),
         Replicate=ifelse(str_detect(File,fixed('rep2')),2,Replicate)) %>%
  ungroup %>%
    mutate(Condition=str_trim(Condition),
           Condition=ifelse(Condition %in% names(conditions),conditions[Condition],Condition)) %>%
  filter(!is.na(Abs) & !(Type=='RNAseq' & Gene %in% c('cpt-2','cpt-5','atgl-1'))) %>%
  left_join(info) %>%
  mutate(SGroup=ifelse(Supplement=='None',Strain,paste(Strain,str_sub(Supplement, 1, 3),sep='-'))) %>%
  mutate_at(c('Type','Gene','Strain','Metformin_mM','Supplement','SGroup','Condition','ID','Replicate'),as.factor) %>%
  mutate(SGroup=factor(SGroup,levels=c('OP50','OP50-MR','crp','cra','OP50-Glu')),
         Log=log2(Abs)) %>%
  gather(Measure,Value,Abs,Log) %>%
  select(Type,Gene,Condition,Replicate:Value) %>%
  group_by(Gene,Type,Replicate,Measure) %>%
  mutate(Norm_Rep=Value-mean(Value[Strain=='OP50' & Supplement=='None'])) %>%
  group_by(Gene,Type,Measure) %>%
  mutate(Norm=Norm_Rep-mean(Norm_Rep[Strain=='OP50' & Supplement=='None'])) 


data %>%
  filter(Measure=='Log') %>%
  group_by(Type,Gene,Replicate,Strain,Supplement,Metformin_mM) %>%
  summarise(Count=n()) %>%
  View()
  

unique(data$Gene)
  

data %>%
  group_by(Condition) %>%
  summarise %>%
  write_csv('Conditions_raw.csv')



plotBox<-function(data,yval,ylb){
  data %>%
  ggplot+
    aes_string(x="SGroup",y=yval,color="Metformin_mM")+
    stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
    geom_jitter(width=0.25)+
    ylab(ylb)+
    xlab('Bacterial strain + Supplement')+
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



#Linear modelling
#
#

# model<-lm(LogValue~0+ID, data=data)
# lmod_glht <- multcomp::glht(model, linfct = contr.matrix)
# result<-multcomp:::summary.glht(lmod_glht,test=multcomp::adjusted("none"))
# res<-data.frame(result$test[c('coefficients','sigma','tstat','pvalues')])


allgroups<-as.character(unique(data$ID))
allgroups

contrasts<-read.contrasts2('!Contrasts.xlsx')

contrasts$Contrasts.table
contrasts.desc<-contrasts$Contrasts.table%>%
  select(Description:Supplement)


contr.matrix<-contrasts$Contrasts.matrix
contr.matrix




lmdata<-data %>%
  filter(Measure=='Log') %>%
  group_by(Gene) %>%
  do(hypothesise2(.,"Norm~0+ID",contr.matrix)) %>%
  ungroup


results<-contrasts.desc %>%
  mutate(Description=factor(Description,levels=contrasts.desc$Description)) %>%
  left_join(lmdata) %>%
  #Adjustments within contrast
  adjustments %>%
  select(Gene,everything())



results %>%
  filter(Contrast=='OP50_T') %>%
  arrange(FDR) %>%
  View


results.m<-results %>%
  gather(Stat,Value,logFC:logFDR)

results.castfull<-results.m %>%
  arrange(Contrast,desc(Stat)) %>%
  unite(CS,Contrast,Stat) %>%
  select(Gene,CS,Value) %>%
  spread(CS,Value) %>%
  mutate_at(vars(-c(Gene),-contains('Stars')),as.numeric)

results.cast<-results.m %>%
  filter(Stat %in% c('logFC','FDR')) %>%
  #arrange(Contrast,desc(Stat)) %>%
  unite(CS,Contrast,Stat) %>%
  select(Gene,CS,Value) %>%
  spread(CS,Value) %>%
  mutate_at(vars(-c(Gene),-contains('Stars')),as.numeric)

head(results.castfull)
head(results.cast)

View(results.cast)



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
  #filter(Contrast %in% c('OP50_T','OP50-MR_T','OP50-MR_I')) %>%
  ggplot(aes(x=Description,y=Gene))+
  geom_tile(aes(fill=logFC))+
  #geom_point(aes(size=FDRc,colour=logFC),alpha=0.9)+
  theme_minimal()+
  geom_text(aes(label=as.character(pStars)))+
  #scale_size_discrete(range = c(2,4))+#,breaks=brks
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  #scale_fill_gradient2(low = "purple", mid = "gray", high = "red", midpoint = 0, breaks = clrbrks)+
  xlab("Comparison")+
  theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_All.pdf',sep = ''),
             width=8,height=6,useDingbats=FALSE)



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





amp<-5

minv<- -amp
maxv<- amp

nstep<-maxv-minv
nstep<-8

clrbrks<-seq(-amp,amp,by=1)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)



RNAresults %>%
  filter(Gene %in% c('acs-2','F37H8.3', 'dhs-23', 'fat-7', 'cpt-2', 'cpt-5')) %>%
  ggplot(aes(x=Description,y=Gene))+
  geom_tile(aes(fill=logFC))+
  geom_text(aes(label=as.character(Stars)))+
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  xlab("Comparison")+
  facet_grid(~Type)+
  theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = NA),
        axis.line.x = element_line(colour = NA),
        axis.line.y = element_line(colour = NA),
        strip.text = element_text(colour = 'black', face='bold',size=10),
        axis.text.x= element_text(face='bold', colour='black', size=10, angle = 90, hjust = 1),
        axis.text.y= element_text(face='bold', colour='black', size=10))


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_RNAseq_clean.pdf',sep = ''),
             width=4,height=5,useDingbats=FALSE)














library(tidyverse)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)


setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy/")


odir<-'Summary_Glycerol'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#load("Fluorescence_Glycerol.RData")
#save.image('Fluorescence_Glycerol.RData')


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


data <- read_csv('All_raw_data.csv') %>%
  filter(Type %in% c('Glycerol')) %>%
  mutate_at(c('Type','Gene','Strain','Metformin_mM','Supplement','SGroup','Condition','ID','Replicate'),as.factor) %>%
  mutate(SGroup=factor(SGroup,levels=c('OP50','OP50-MR','crp','cra','glpK','OP50-Glu','OP50-Gly','glpK-Gly')),
         Strain=factor(Strain,levels=c('OP50','glpK')),
         Supplement=factor(Supplement,levels=c('None','Glycerol'),labels=c('','+ Glycerol') ))



data %>%
  filter(Measure=='Log') %>%
  group_by(Type,Gene,Replicate,Strain,Supplement,Metformin_mM) %>%
  summarise(Count=n()) %>%
  View()
  

unique(data$Gene)
  




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
     width=10,height=12, useDingbats=FALSE, ggsave)


Fluorplots<-data %>%
  group_by(Measure) %>%
  do(plot=plotBox(.,"Norm",as.character(unique(.$Measure))) )


map2(paste0(odir,"/Fluorescence_norm_",as.character(Fluorplots$Measure),".pdf"),
     Fluorplots$plot,
     width=10,height=12, useDingbats=FALSE, ggsave)
  

Fluorplots<-data %>%
  group_by(Measure) %>%
  do(plot=plotBoxC(.,"Norm",as.character(unique(.$Measure))) )


map2(paste0(odir,"/Fluorescence_norm_joined_",as.character(Fluorplots$Measure),".pdf"),
     Fluorplots$plot,
     width=6,height=12, useDingbats=FALSE, ggsave)


data %>%
  filter(Measure=='Log') %>%
  ggplot(aes(x=interaction(Strain,Supplement,sep = ' '),y=Norm,color=Metformin_mM))+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_jitter(width=0.25)+
  ylab('C. elegans fluorescence (normalised), log2 A.U.')+
  xlab('Bacterial strain + Supplement')+
  labs(color='Metformin, mM')+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(~Gene,scale="free_y")


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Fluorescence_comparison.pdf",sep=''),
             width=10,height=6, useDingbats=FALSE)



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

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=7,height=10, useDingbats=FALSE)



sum.c %>%
  filter(VarPrc>300)


ggplot(sum.c,aes(x=Index,y=VarPrc))+
  geom_point()+
  scale_y_continuous(breaks = seq(0,1000,by=25),limits=c(0,100))+
  ylab('In-group variation, %')+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=7,height=10, useDingbats=FALSE)



#Linear modelling


contrasts<-read.contrasts2('!Contrasts_fluorescence.xlsx')

contrasts$Contrasts.table
contrasts.desc<-contrasts$Contrasts.table%>%
  select(Description:Supplement)


contr.matrix<-contrasts$Contrasts.matrix
contr.matrix

results.all<-data %>%
  filter(Measure=='Log') %>%
  group_by(Gene) %>%
  do(hypothesise2(.,"Norm~0+ID",contr.matrix)) %>%
  getresults(contrasts.desc)


results<-results.all$results
results.castfull<-results.all$castfull
results.cast<-results.all$cast
results.multi<-results.all$multi

head(results.castfull)
head(results.cast)

View(results)

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


d<-dist(as.matrix(heatsum),method = "euclidean")
h<-hclust(d)
ordmet<-rownames(heatsum[h$order,])

if (length(ordmet)!=length(unique(ordmet))){
  print("Non unique metabolites!")
}


max(heatsum)
min(heatsum)

amp<-6

minv<- -amp
maxv<- amp

nstep<-maxv-minv
nstep<-8

clrbrks<-seq(-amp,amp,by=2)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)




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
             width=6,height=4,useDingbats=FALSE)



descriptions<-c('SM-S'='Treatment effect on OP50','RM-R'='Treatment effect on OP50-MR','SM-S-(RM-R)'='Interaction between treatment and OP50-MR')
translations<-c('F44G3.2'='argk-1','C05D11.7'='atgl-1')

genes<-unique(as.character(results$Gene))














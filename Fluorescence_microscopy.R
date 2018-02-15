library(tidyverse)



devtools::install_github("PNorvaisas/PFun")
library(PFun)


setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy/")



odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



#load("FluorescenceMicroscopy.RData")
save.image('FluorescenceMicroscopy.RData')


theme_set(theme_light())

transgenes<-c('cpt-2')

info<-read_csv('Conditions.csv')

data<-data.frame(Gene=transgenes) %>%
  group_by(Gene) %>%
  do(read_delim(paste0(.$Gene,' example.txt'),delim='\t') ) %>%
  ungroup %>%
  gather(Condition,Value,-Gene) %>%
  mutate(Condition=trimws(Condition)) %>%
  filter(!is.na(Value)) %>%
  left_join(info) %>%
  mutate(SGroup=ifelse(Supplement=='None',Strain,paste(Strain,str_sub(Supplement, 1, 3),sep='-'))) %>%
  mutate_at(c('Gene','Strain','Metformin_mM','Supplement','SGroup','Condition','ID'),as.factor) %>%
  mutate(SGroup=factor(SGroup,levels=c('OP50','OP50-MR','crp','cra','OP50-Glu')),
         LogValue=log2(Value)) %>%
  select(Gene:Condition,ID:SGroup,Value,LogValue)



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
    facet_grid(~Gene)
}

data %>%
  plotBox("Value","Normalised fluorescence, A.U.")

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Fluorescence.pdf",sep=''),
             width=12,height=9)

data %>%
  plotBox("LogValue","Log Normalised fluorescence, A.U.")

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/LogFluorescence.pdf",sep=''),
             width=12,height=9)


#Linear modelling

allgroups<-as.character(unique(data$ID))
allgroups

contrasts<-read.contrasts('!Contrasts.xlsx','Contrasts_values',allgroups)

contrasts$Contrasts.table
contrasts.desc<-contrasts$Contrasts.table%>%
  select(Contrast,Description,Contrast_type,Strain,Supplement)


contr.matrix<-contrasts$Contrasts.matrix
contr.matrix




lmdata<-data %>%
  group_by(Gene) %>%
  do(hypothesise2(.,"LogValue~0+ID",contr.matrix)) %>%
  ungroup


model<-lm(LogValue~0+ID, data=data)
lmod_glht <- multcomp::glht(model, linfct = contr.matrix)
result<-multcomp:::summary.glht(lmod_glht,test=multcomp::adjusted("none"))
res<-data.frame(result$test[c('coefficients','sigma','tstat','pvalues')])


results<-contrasts.desc %>%
  left_join(lmdata) %>%
  #Adjustments within contrast
  group_by(Contrast) %>%
  mutate(FDR=p.adjust(p.value,method = 'fdr'),
         PE=logFC+SE,
         NE=logFC-SE,
         logFDR=-log10(FDR)) %>%
  ungroup %>%
  select(Gene,everything())



results %>%
  filter(Contrast %in% c('OP50_T','OP50Glu_I')) %>%
  arrange(FDR)



results.m<-results %>%
  gather(Stat,Value,logFC:logFDR)

results.castfull<-results.m %>%
  arrange(Contrast,desc(Stat)) %>%
  unite(CS,Contrast,Stat) %>%
  select(Gene,CS,Value) %>%
  spread(CS,Value)

results.cast<-results.m %>%
  filter(Stat %in% c('logFC','FDR')) %>%
  arrange(Contrast,desc(Stat)) %>%
  unite(CS,Contrast,Stat) %>%
  select(Gene,CS,Value) %>%
  spread(CS,Value)

head(results.castfull)
head(results.cast)



write.csv(results,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)



#Generate table for heatmap
heatsum<-results %>%
  filter(Description %in% comparisons & Metabolite %in% sel.mets$Metabolite) %>%
  select(Description,Metabolite,logFC) %>%
  spread(Description,logFC) %>%
  data.frame(check.names = FALSE,check.rows = FALSE)


rownames(heatsum)<-heatsum$Metabolite
heatsum$Metabolite<-NULL



max(heatsum)
min(heatsum)

amp<-8

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


results.sum<-results %>%
  filter(Metabolite %in% ordmet & Description %in% comparisons) %>%
  mutate(Metabolite=factor(Metabolite,levels=ordmet,labels=ordmet),
         Description=factor(Description,levels=comparisons,labels=comparisons),
         FDRstars=stars.pval(FDR))


ggplot(results.sum,aes(x=Description,y=Metabolite))+
  geom_tile(aes(fill=logFC))+
  #geom_point(aes(size=FDRc,colour=logFC),alpha=0.9)+
  theme_minimal()+
  geom_text(aes(label=as.character(FDRstars)))+
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


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Complete_tidy.pdf',sep = ''),
             width=6,height=16,useDingbats=FALSE)








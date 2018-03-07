library(tidyverse)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)


setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy/")


#odir<-'Summary_Transgenes'
#dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#load("Fluorescence_Collect.RData")
#save.image('Fluorescence_Collect.RData')




info<-read_csv('Conditions_all.csv')


#Translations for consistency
#CRP and RNAseq
conditions<-c("glucose 0"="c0_Glu","glucose 50"="c50_Glu",
              "glu0"="c0_Glu","glu50"="c50_Glu",
              "r0"="mr0","r50"="mr50","C0"="c0","C50"="c50","R0"="mr0","R50"="mr50",
              "glpK 0 mM"="glpK0","glpK 50 mM"="glpK50","OP50-C 0 mM"="c0","OP50-C 50 mM"="c50",
              "glpK 0 mM + Glycerol"="glpK0_Gly","glpK 50 mM + Glycerol"="glpK50_Gly",
              "OP50-C 0 mM + Glycerol"="c0_Gly","OP50-C 50 mM + Glycerol"="c50_Gly")


data<-data.frame(Type=c("RNAseq","CRP","Glycerol"),Folder=c("RNAseq","crp-cra-glucose","Glycerol")) %>%
  group_by(Folder,Type) %>%
  do(data.frame(File=list.files(path = paste0("./Data/",.$Folder)))) %>%
  mutate(Gene=str_replace_all(File,'_rep[[:digit:]]|.txt|.csv','')) %>%
  group_by(Folder,Type,File,Gene) %>%
  do(read_delim(paste('./Data',.$Folder,.$File,sep='/'),delim=ifelse(str_detect(.$File,fixed('.csv')),',','\t')) %>% gather(Condition,Abs,everything()) ) %>%
  group_by(Folder,File,Type,Gene,Condition) %>%
  mutate(Rank=row_number(),
         Replicate=as.integer(ifelse(ifelse( Type=='Glycerol',Rank>=22,Rank>=31 ),2,1)),
         Replicate=ifelse(Rank>64,3,Replicate),
         Replicate=ifelse(str_detect(File,fixed('rep2')),2,Replicate)) %>%
  ungroup %>%
    mutate(Condition=str_trim(Condition),
           Condition=ifelse(Condition %in% names(conditions),conditions[Condition],Condition)) %>%
  filter(!is.na(Abs) & !(Type=='RNAseq' & Gene %in% c('cpt-2','cpt-5','atgl-1'))) %>%
  left_join(info) %>%
  mutate(SGroup=ifelse(Supplement=='None',Strain,paste(Strain,str_sub(Supplement, 1, 3),sep='-'))) %>%
  mutate_at(c('Type','Gene','Strain','Metformin_mM','Supplement','SGroup','Condition','ID','Replicate'),as.factor) %>%
  mutate(SGroup=factor(SGroup,levels=c('OP50','OP50-MR','crp','cra','glpK','OP50-Glu','OP50-Gly','glpK-Gly')),
         Log=log2(Abs)) %>%
  gather(Measure,Value,Abs,Log) %>%
  select(Type,Gene,Condition,Replicate:Value) %>%
  group_by(Gene,Type,Replicate,Measure) %>%
  mutate(Norm=Value-mean(Value[Strain=='OP50' & Supplement=='None'])) 

data %>%
  write_csv('All_raw_data.csv')



data %>%
  filter(Measure=='Log') %>%
  group_by(Type,Gene,Replicate,Strain,Supplement,Metformin_mM) %>%
  summarise(Count=n()) %>%
  write_csv('All_data_Worm_count.csv')



# data %>%
#   group_by(Condition) %>%
#   summarise %>%
#   write_csv('Conditions_raw_all.csv')



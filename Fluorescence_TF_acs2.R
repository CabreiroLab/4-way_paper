library(tidyverse)
library(PFun)

theme_set(theme_light())

setwd("~/Dropbox/Projects/Metformin_TF_acs-2")

odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

#load('Metformin_TF_acs-2.RData')
#save.image('Metformin_TF_acs-2.RData')



translations<-c('gvcA'='gcvA','op50-c'='OP50-C','op50'='OP50-C')

allfiles<-read_csv('Allfiles.csv') %>%
  filter( !str_detect(Folder,'grouped') ) %>%
  mutate(Type=ifelse(str_detect(Folder,'wo metf'),'C','T' ),
         Replicate=as.integer(ifelse(str_detect(Replicate_folder,'Rep1'),1,2)),
         Gene=str_trim(str_replace_all(Folder,'metf|wo metf|later|met repeat','')),
         Gene=ifelse(Gene %in% names(translations),translations[Gene],Gene ),
         FileNo=str_extract(File,'[[:digit:]]{1,4}.tiff') %>% str_replace('.tiff','') %>% str_pad(width=4,pad='0') ) %>%
  group_by(Replicate,Gene,Type,Folder) %>%
  arrange(FileNo) %>%
  mutate(FileInd=row_number()) %>%
  ungroup %>%
  #filter(!FileNo=='0001') %>%
  group_by(Replicate,Gene,Type,Folder) %>%
  mutate(Count=n()) %>%
  ungroup %>%
  group_by(Gene) %>%
  mutate(Max=as.integer(max(Count))) %>%
  select(Replicate,Type,Gene,FileNo,FileInd,Replicate_folder:File,Max) %>%
  ungroup

head(allfiles)

#[[:digit:]]{1,4}

allfiles %>%
  filter(Gene=='OP50-C') %>%
  head


allfiles %>%
  group_by(Replicate,Gene,Type) %>%
  summarise(Count=n()) %>%
  ungroup %>%
  unite(RT,Replicate,Type,remove = FALSE) %>%
  select(Gene,RT,Count) %>%
  spread(RT,Count)
  
allfiles_sum<-allfiles %>%
  group_by(Replicate,Gene,Type,Folder) %>%
  summarise(Count=n()) %>%
  ungroup %>%
  group_by(Gene) %>%
  summarise(Max=as.integer(max(Count)))

allfiles_sum

write_csv(allfiles,'Allfiles_annotated.csv')

write_csv(allfiles_sum,'Allfiles_summary.csv')




genes<-c('OP50-C','arcA','argR','cra','crp','csiR','fur','gcvA','mlc','ntrC')

remindex<-c('1_crp_T_8')

data<-allfiles %>%
  mutate(FileName=File,
         File=paste('.',Replicate_folder,Folder,FileName,sep='/')) %>%
  left_join(read_csv('All_results_adaptive.csv') %>%
              mutate(X1=NULL) ) %>%
  mutate_at(c('Replicate','Type','Gene','Worm'),as.factor) %>%
  mutate(Gene=factor(Gene,levels=genes)) %>%
  filter(W_N<100000,
         W_N>40000,
         FileNo!='0001',
         !(Replicate==1 & FileInd==2),
         Folder!='OP50-C metf later') %>%
  group_by(Replicate) %>%
  mutate(W_Int=W_Mean-B_Mean,
         W_LogInt=log2(W_Int),
         TpRep_Mean=mean(W_LogInt[Gene=='OP50-C']),#
         W_NormLog=W_LogInt-TpRep_Mean) %>%
  unite(Index,Replicate,Gene,Type,FileInd,remove = FALSE) %>%
  filter(!Index %in% remindex)  %>%
  ungroup 




data %>%
  filter(Gene=='crp') %>%
  arrange(desc(W_NormLog))


data %>%
  filter(Gene=='OP50-C') %>%
  group_by(Gene,Folder) %>%
  summarise



data %>%
  #filter(Type=='T') %>%
  ggplot(aes(x=Gene,y=W_NormLog,color=Type))+ #,color=Replicate
  geom_jitter() +
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "dodge",alpha=0.2)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_T_Comparison_Log.pdf",sep=''),
             width=10,height=5)

data %>%
  filter(Type=='T') %>%
  pull(W_N) %>%
  mean

data %>%
  filter(Type=='T') %>%
  ggplot(aes(x=W_N))+ #,color=Replicate
  geom_histogram() 


fit<-lm(W_NormLog~Type*Gene,data)

summary(fit)


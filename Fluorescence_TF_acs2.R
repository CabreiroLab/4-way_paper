library(tidyverse)


translations<-c('gvcA'='gcvA','op50-c'='OP50-C','op50'='OP50-C')

allfiles<-read_csv('Allfiles.csv') %>%
  filter( !str_detect(Folder,'grouped') ) %>%
  mutate(Type=ifelse(str_detect(Folder,'wo metf'),'C','T' ),
         Replicate=as.integer(ifelse(str_detect(Replicate_folder,'Rep1'),1,2)),
         Gene=str_trim(str_replace_all(Folder,'metf|wo metf|later|met repeat','')),
         Gene=ifelse(Gene %in% names(translations),translations[Gene],Gene ),
         FileNo=str_pad(str_extract(File,'[[:digit:]]{1,2}'),width=2,pad='0')) %>%
  arrange(Replicate,Gene,Type,FileNo) %>%
  filter(!FileNo=='01') %>%
  mutate(FileInd=row_number()) %>%
  select(Replicate,Type,Gene,FileNo,FileInd,Replicate_folder:File) 

head(allfiles)



allfiles %>%
  filter(FileNo=='01') %>%
  View()


allfiles %>%
  group_by(Replicate,Gene,Type) %>%
  summarise(Count=n()) %>%
  ungroup %>%
  unite(RT,Replicate,Type,remove = FALSE) %>%
  select(Gene,RT,Count) %>%
  spread(RT,Count)
  
allfiles_sum<-allfiles %>%
  group_by(Replicate,Gene,Type) %>%
  summarise(Count=n()) %>%
  ungroup %>%
  group_by(Gene) %>%
  summarise(Max=as.integer(max(Count)))


write_csv(allfiles,'Allfiles_annotated.csv')

write_csv(allfiles_sum,'Allfiles_summary.csv')

library(tidyverse)
library(broom)
library(ggrepel)


# devtools::install_github("PNorvaisas/PFun")
library(PFun)


theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau


cwd <- "~/Dropbox/Projects/Metformin_project/Bacterial Growth Assays/"
setwd(cwd)


odir <- "Summary_IPTG_oeCRP_cyaA"
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


# save.image('BGA_IPTG_oeCRP.RData')
# load('BGA_IPTG_oeCRP.RData')

strainlist <- c("OP50-C", "crp", "cyaA")
overexp <- c("None", "crp")
iptgmm <- c("0", "10", "25", "50")








mutdata <- read_csv("Resistance_mutants_growth_assays/2017-Rosie_Mutants/Summary.csv") %>%
  filter(Strain %in% c("crp") & Metformin_mM == 0) %>%
  mutate(
    Overexpression = "None",
    IPTG_mM = 0,
    Media='LB'
  ) %>%
  select(File=Plate,Media,IPTG_mM,Overexpression,Well:TReplicate, logAUC = Int_600nm_log, ODph = a_log,-c(Metformin_mM,Data) )


mutdata %>%
  filter(IPTG_mM == 0 & Strain == "crp")

mutdata2 <- mutdata %>%
  mutate(Media='NGM')

colnames(data)

data <- read_csv("IPTG titration/Data/Summary.csv") %>%
  select(File:TReplicate, logAUC = `600nm_f_logAUC`, ODph = `600nm_dt_Max`, -Data) %>%
  bind_rows(mutdata) %>%
  bind_rows(mutdata2) %>% #Filler
  filter(!is.na(Strain)) %>%
  filter(!(Strain == "cyaA" & Overexpression == "None" & Replicate == 3 & Media=='LB')) %>%
  mutate_at(vars(File:TReplicate), as.factor) %>%
  mutate(
    Media = fct_relevel(Media,c('NGM','LB')),
    IPTG_mM = fct_relevel(IPTG_mM, iptgmm),
    Strain = fct_relevel(Strain, strainlist),
    Overexpression = fct_relevel(Overexpression, overexp),
    logDph = log2(ODph)
  ) %>%
  gather(Measure, Value, logDph, logAUC) %>%
  group_by(Media,Measure, IPTG_mM, Overexpression) %>%
  mutate(Ref_C = mean(Value[Strain == "OP50-C"])) %>%
  group_by(Media,Measure, Strain, Overexpression) %>%
  mutate(Ref_I = mean(Value[IPTG_mM == "0"])) %>%
  group_by(Media,Measure, Strain, IPTG_mM) %>%
  mutate(Ref_O = mean(Value[Overexpression == "None"])) %>%
  group_by(Media,Measure, Overexpression) %>%
  mutate(Ref_CI = mean(Value[ Strain == "OP50-C" & IPTG_mM == "0"])) %>%
  group_by(Media,Measure, IPTG_mM) %>%
  mutate(Ref_CO = mean(Value[Strain == "OP50-C" & Overexpression == "None"])) %>%
  group_by(Media,Measure, Strain) %>%
  mutate(Ref_IO = mean(Value[IPTG_mM == "0" & Overexpression == "None"])) %>%
  group_by(Media,Measure) %>%
  mutate(
    Ref_CIO = mean(Value[Strain == "OP50-C" & IPTG_mM == "0" & Overexpression == "None"]),
    Norm_C = Value - Ref_C,
    Norm_I = Value - Ref_I,
    Norm_O = Value - Ref_O,
    Norm_CI = Value - Ref_CI,
    Norm_CO = Value - Ref_CO,
    Norm_IO = Value - Ref_IO,
    Norm_CIO = Value - Ref_CIO
  ) %>%
  ungroup() %>%
  select(-contains("Ref_")) %>%
  gather(Normalisation, Value, Value, Norm_C, Norm_I, Norm_O, Norm_CI, Norm_CO, Norm_IO, Norm_CIO) %>%
  mutate(Value = ifelse(Value %in% c(Inf, -Inf), NA, Value))






data %>%
  write_csv(paste0(odir, "/Raw_data_Summary.csv"))






mutdatats <- read_csv("Resistance_mutants_growth_assays/2017-Rosie_Mutants/Data.csv") %>%
  filter(Strain %in% c("crp") & Metformin_mM == 0) %>% #Exclude OP50-C data
  mutate(
    Overexpression = "None",
    IPTG_mM = 0,
    Media='LB'
  ) %>%
  gather(Time_s, OD, `0.0`:`64800.0`) %>%
  select(Data,File=Plate,Media,IPTG_mM,Overexpression,Well:TReplicate, Time_s,OD )



mutdatats2<-mutdatats %>%
  mutate(Media='NGM')
  
  
colnames(data_ts)
unique(mutdatats$Data)  
unique(data_ts$Data)  
  
data_ts <- read_csv("IPTG titration/Data/Timeseries.csv") %>%
  filter(!is.na(Strain)) %>%

  filter(!(Strain == "cyaA" & Overexpression == "None" & Replicate == 3)) %>%
  gather(Time_s, OD, `0`:`64800`) %>%
  bind_rows(mutdatats) %>%
  bind_rows(mutdatats2) %>%
  filter(Data == "600nm_f") %>%
  mutate_at(vars(File:TReplicate), as.factor) %>%
  select(-Data) %>%
  mutate(
    Media=fct_relevel(Media,c('NGM','LB')),
    IPTG_mM = fct_relevel(IPTG_mM, iptgmm),
    Strain = fct_relevel(Strain, strainlist),
    Overexpression = fct_relevel(Overexpression, overexp),
    Time_s = as.numeric(Time_s),
    Time_h = Time_s / 3600
  )






unique(data$IPTG_mM)
unique(data$Overexpression)
unique(data$Strain)


data %>%
  filter(is.infinite(Value))


data.sum <- data %>%
  group_by(Media,Measure, Normalisation, Strain, Overexpression, IPTG_mM) %>%
  summarise(
    Mean = mean(Value),
    SD = sd(Value),
    SE = SD / sqrt(n())
  ) %>%
  mutate(
    PE = Mean + SE,
    NE = Mean - SE,
    Prc = ifelse(Normalisation == "Value", NA, 2^Mean * 100),
    PrcNE = ifelse(Normalisation == "Value", NA, 2^NE * 100),
    PrcPE = ifelse(Normalisation == "Value", NA, 2^PE * 100)
  )


datats.sum <- data_ts %>%
  group_by(Media,Strain, Overexpression, IPTG_mM, Time_h) %>%
  summarise(
    OD_Mean = mean(OD),
    OD_SD = sd(OD)
  )





data %>%
  group_by(Media,Measure, Normalisation, Strain, Overexpression, IPTG_mM) %>%
  summarise(N = n()) %>%
  View()





data_ts %>%
  group_by(Media,Strain, Overexpression, IPTG_mM) %>%
  summarise(Count = n()) %>%
  View()


IPTGcols <- colorRampPalette(c("red", "blue4"))(4)
names(IPTGcols) <- levels(data$IPTG_mM)
IPTGlab <- "IPTG,\nmM"

ggplot(data_ts, aes(x = Time_h, y = OD, color = IPTG_mM)) +
  geom_line(aes(group = interaction(File,Media, Replicate, TReplicate, IPTG_mM))) +
  xlab("Time, h") +
  scale_colour_manual(name = IPTGlab, values = IPTGcols) +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  facet_grid(Media+Overexpression ~ Strain, labeller = label_both)


ggsave(
  file = paste(odir, "/Growth_overview.pdf", sep = ""),
  width = 110, height = 160, units = "mm", scale = 2, device = cairo_pdf, family = "Arial"
)







datats.sum %>%
  ggplot(aes(x = Time_h, y = OD_Mean, color = IPTG_mM, fill = IPTG_mM)) +
  geom_line() +
  geom_ribbon(aes(
    ymin = OD_Mean - OD_SD,
    ymax = OD_Mean + OD_SD
  ), alpha = 0.5, color = NA) +
  xlab("Time, h") +
  ylab("OD") +
  scale_colour_manual(name = IPTGlab, values = IPTGcols) +
  scale_fill_manual(name = IPTGlab, values = IPTGcols) +
  scale_x_continuous(breaks = seq(0, 18, by = 6)) +
  facet_grid(Media+Overexpression ~ Strain, labeller = label_both)

ggsave(
  file = paste0(odir, "/Growth_Summary.pdf"),
  width = 110, height = 160, units = "mm", scale = 2, device = cairo_pdf, family = "Arial"
)



# stats <- data %>%
#   group_by(Measure, Normalisation) %>%
#   filter(!is.na(Value)) %>%
#   do(tidy(lm(Value ~ Strain, data = .))) %>%
#   filter(term != "(Intercept)") %>%
#   rename(
#     Strain = term,
#     SE = std.error,
#     logFC = estimate
#   ) %>%
#   mutate(
#     Strain = str_replace(Strain, "Strain", ""),
#     Strain = factor(Strain, levels = strainlist, labels = strainlist),
#     PE = logFC + SE,
#     NE = logFC - SE,
#     Prc = 2^logFC * 100,
#     PrcNE = 2^NE * 100,
#     PrcPE = 2^PE * 100,
#     pStars = pStars(p.value)
#   )



#Test recipes
# library(recipes)
# library(caret)
# 
# rec<-data %>%
#   filter(Measure=="logAUC" & Normalisation=="Value") %>%
#   recipe(Value ~ Strain+Overexpression+IPTG_mM,data=.) %>%
#   step_dummy(Strain,Overexpression,IPTG_mM) %>%
#   step_interact(~contains("Strain"):contains("Overexpression") ) %>%
#   step_interact(~contains("Strain"):contains("Overexpression"):contains("IPTG_mM") )%>%
#   prep(training = data,retain=TRUE)
# 
# summary(rec)
# 
# bake(rec,newdata = data)
# 
# fit<-train(rec,
#            data=data %>% filter(Measure=="logAUC" & Normalisation=="Value"),
#            method='lm')
# 
# summary(fit)



# Needs further updates
#summary(fit)


stats.O <- data %>%
  filter(Normalisation == "Value" & !is.na(Value) & Strain %in% c("OP50-C","cyaA")) %>%
  group_by(Media,Measure,Strain,IPTG_mM) %>%
  do(tidy(lm(Value ~ Overexpression, data = .,na.action = na.omit))) %>%
  mutate(
    Comparison = case_when(
      str_detect(term, "\\:") ~ "Interaction",
      str_detect(term, "Strain") ~ "Strain",
      str_detect(term, "Overexpression") ~ "Overexpression"
    ),
    Overexpression=str_replace_all(term,"Overexpression","")
  )


stats.S <- data %>%
  filter(Normalisation == "Value" & !is.na(Value) ) %>%
  group_by(Media,Measure,Overexpression,IPTG_mM) %>%
  do(tidy(lm(Value ~ Strain, data = .))) %>%
  mutate(
    Comparison = case_when(
      str_detect(term, "\\:") ~ "Interaction",
      str_detect(term, "Strain") ~ "Strain",
      str_detect(term, "Overexpression") ~ "Overexpression"
    ),
    Strain=str_replace_all(term,"Strain","")
  )

stats.I <- data %>%
  filter(Normalisation == "Value" & !is.na(Value) & !(Strain=='crp' & Overexpression=='None')) %>%
  group_by(Media,Measure,Overexpression,Strain) %>%
  do(tidy(lm(Value ~ IPTG_mM, data = .))) %>%
  mutate(
    Comparison = case_when(
      str_detect(term, "\\:") ~ "Interaction",
      str_detect(term, "IPTG_mM") ~ "IPTG_mM",
      str_detect(term, "IPTG_mM") ~ "IPTG_mM"
    ),
    IPTG_mM=str_replace_all(term,"IPTG_mM","")
  )


stats.SI <- data %>%
  filter(Normalisation == "Value" & !is.na(Value) ) %>%
  group_by(Media,Measure,Overexpression) %>%
  do(tidy(lm(Value ~ Strain*IPTG_mM, data = .))) %>%
  filter(str_detect(term, "\\:")) %>%
  mutate(
    Comparison = case_when(
      str_detect(term, "\\:") ~ "Interaction_Strain-IPTG_mM"
    ),
    term=str_replace_all(term,"Strain|IPTG_mM","")
    
  ) %>%
  separate(term,c("Strain","IPTG_mM"),remove = FALSE)


#Needs crp reference
data %>%
  filter(Normalisation == "Value" & !is.na(Value) & Measure=='logAUC' ) %>%
  group_by(Media,Measure,Strain,Overexpression,IPTG_mM) %>%
  summarise(n=n()) %>%
  View



stats.OI <- data %>%
  filter(Normalisation == "Value" & !is.na(Value) ) %>%
  group_by(Media,Measure,Strain) %>%
  do(tidy(lm(Value ~ Overexpression*IPTG_mM, data = .))) %>%
  filter(str_detect(term, "\\:")) %>%
  mutate(
    Comparison = case_when(
      str_detect(term, "\\:") ~ "Interaction_Overexpression-IPTG_mM"
    ),
    term=str_replace_all(term,"Overexpression|IPTG_mM","")
  ) %>%
  separate(term,c("Overexpression","IPTG_mM"),remove = FALSE)




stats.cyaA.I <- data %>%
  filter(Normalisation == "Value" & !is.na(Value) & Strain %in% c("OP50-C", "cyaA")) %>%
  group_by(Media,Measure, IPTG_mM) %>%
  do(tidy(lm(Value ~ Strain * Overexpression, data = .)))


stats.cyaA.O <- data %>%
  filter(Normalisation == "Value" & !is.na(Value) & Strain %in% c("OP50-C", "cyaA")) %>%
  group_by(Media,Measure, Strain, IPTG_mM) %>%
  do(tidy(lm(Value ~ Overexpression, data = .)))



colnames(stats.O)
colnames(stats.S)
colnames(stats.SI)
colnames(stats.OI)

stats <- stats.O %>%
  bind_rows(stats.S) %>%
  bind_rows(stats.I) %>%
  bind_rows(stats.SI) %>%
  bind_rows(stats.OI) %>%
  filter(term != "(Intercept)") %>%
  rename(
    SE = std.error,
    logFC = estimate,
    p = p.value
  ) %>%
  mutate(
    PE = logFC + SE,
    NE = logFC - SE,
    Prc = 2^logFC * 100,
    PrcNE = 2^NE * 100,
    PrcPE = 2^PE * 100,
    pStars = pStars(p)
  ) %>%
  select(Media,Measure, Comparison, Strain, Overexpression, IPTG_mM, everything(), -term)






View(stats)



stats %>%
  write_csv(paste0(odir, "/Stats_Summary.csv"))


View(stat)


stat %>%
  filter(Measure == "AUC" & Normalisation == "Norm_CM" & Metformin_mM == "150")




data.sum %>%
  filter(Measure == "logAUC" & Normalisation == "Norm_CIO") %>%
  ggplot(aes(x = IPTG_mM, y = Prc, color = Overexpression)) +
  geom_hline(aes(yintercept = 100), color = "red4", alpha = 0.5, linetype = "longdash") +
  scale_y_continuous(breaks = seq(0, 400, by = 25)) +
  coord_cartesian(ylim = c(0, 130)) +
  geom_line(aes(group = interaction(Strain, Overexpression))) +
  geom_errorbar(aes(ymin = PrcNE, ymax = PrcPE), width = 0.25) +
  geom_point() +
  # ggtitle("Comparison vs OP50-C in control")+
  # scale_colour_manual(name = Metlab,values = IPTGcols)+
  ylab("Growth logAUC vs OP50-C Control, %") +
  xlab("IPTG, mM") +
  # geom_text(data=stat %>%
  #             filter(Measure=='AUC' & Normalisation=="Norm_C"),
  #           aes(label=pStars),color='black',y = 125)+
  facet_grid(Media~ Strain)


ggsave(
  file = paste0(odir, "/Growth_Comparison_vs_OP50-C_Control.pdf"),
  width = 110, height = 41, units = "mm", scale = 2, device = cairo_pdf, family = "Arial"
)



ovcol <- ggthemes::tableau_color_pal(palette = "tableau20")(12)[1]



unique(stats$Comparison)

data.sum %>%
  filter(Media=="NGM",
         Measure == "logAUC" &
    Normalisation == "Norm_CIO" &
    Strain %in% c("OP50-C", "crp") &
    Overexpression %in% c("None", "crp")) %>%
  ggplot(aes(x = IPTG_mM, y = Prc, linetype = Overexpression, color = Strain)) +
  geom_hline(aes(yintercept = 100), color = "red4", alpha = 0.5, linetype = "longdash") +
  scale_y_continuous(breaks = seq(0, 400, by = 25)) +
  coord_cartesian(ylim = c(50, 130)) +
  geom_line(aes(group = interaction(Strain, Overexpression))) +
  geom_errorbar(aes(ymin = PrcNE, ymax = PrcPE), width = 0.25) +
  geom_point() +
  # geom_text(
  #   data = stats %>%
  #     filter(Measure == "logAUC" & Strain == "crp" & Comparison == "Strain"),
  #   aes(label = pStars), y = 130, show.legend = FALSE
  # ) +
  geom_text(
    data =filter(stats,Media=='NGM',
             Measure == "logAUC" & Strain !="cyaA" & Comparison == "Interaction_Overexpression-IPTG_mM"),
    aes(label = pStars), y = 125, show.legend = FALSE
  ) +
  scale_alpha_discrete(range = c(0.75, 1)) +
  ylab("Growth AUC vs OP50-C Control, %") +
  xlab("IPTG, mM")



ggsave(
  file = paste0(odir, "/Growth_Comparison_vs_OP50-C_Control_crp_NGM.pdf"),
  width = 55, height = 41, units = "mm", scale = 2, device = cairo_pdf, family = "Arial"
)



data.sum %>%
  filter(Media=='NGM',
         Measure == "logAUC" &
    Normalisation == "Norm_CIO" &
    Strain %in% c("OP50-C", "cyaA") &
    Overexpression %in% c("None", "crp")) %>%
  ggplot(aes(x = IPTG_mM, y = Prc, linetype = Overexpression, color = Strain)) +
  geom_hline(aes(yintercept = 100), color = "red4", alpha = 0.5, linetype = "longdash") +
  scale_y_continuous(breaks = seq(0, 400, by = 25)) +
  coord_cartesian(ylim = c(0, 130)) +
  geom_line(aes(group = interaction(Strain, Overexpression))) +
  geom_errorbar(aes(ymin = PrcNE, ymax = PrcPE), width = 0.25) +
  geom_point() +
  # geom_text(
  #   data = stats %>%
  #     filter(Measure == "logAUC" & Strain == "cyaA" & Comparison == "Strain"),
  #   aes(label = pStars), y = 130, show.legend = FALSE
  # ) +
  geom_text(
    data = stats %>%
      filter(Measure == "logAUC" & Strain != "crp" & Comparison == "Interaction_Overexpression-IPTG_mM"),
    aes(label = pStars), y = 125,show.legend = FALSE
  ) +
  ylab("Growth AUC vs OP50-C Control, %") +
  xlab("IPTG, mM")


ggsave(
  file = paste0(odir, "/Growth_Comparison_vs_OP50-C_Control_cyaA_NGM.pdf"),
  width = 55, height = 41, units = "mm", scale = 2, device = cairo_pdf, family = "Arial"
)

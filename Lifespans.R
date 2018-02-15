library(tidyverse)
library(survival)
library(survminer)

data("lung")
head(lung)

setwd("~/Dropbox/Projects/2015-Metformin/Lifespans")

#load("")
#save.image('Lifespans.RData')



theme_set(theme_light())

info<-read_csv('Groups.csv')

data<-read_delim('LS-11-9-19-genotypes.dat',delim=',') %>%
  filter(!is.na(Condition)) %>%
  mutate(Status=as.integer(recode(Censor,"0"="2","1"="1"))) %>%
  left_join(info) %>%
  mutate_at(c("Worm", "Bacteria", "Metformin_mM"),as.factor) %>%
  mutate(Worm=relevel(Worm,ref='N2'),
         Bacteria=relevel(Bacteria,ref='OP50'))


# data %>%
#   group_by(Condition) %>%
#   summarise %>%
#   write_csv('Groups_raw.csv')


data.sel<-data%>%
  filter(Worm=='N2') %>%
  mutate(Bacteria=relevel(Bacteria,ref='OP50'))


data.sel<-data%>%
  filter(Worm=='N2' ) %>%
  mutate(Bacteria=relevel(Bacteria,ref='OP50'),
         Worm=relevel(Worm,ref='N2'))

data %>%
  group_by(Condition,Worm,Bacteria,Metformin_mM) %>%
  summarise(Total=sum(Frequency))


res.cox <- coxph(Surv(Day, Status) ~Bacteria+Metformin_mM, weights = Frequency,data = data.sel)
res.cox


ph2frame<-function(fit){
  summary(fit)$coefficients %>%
    data.frame %>%
    rownames_to_column(var = "Test") %>%
    rename(exp_coef=exp.coef.,SE=se.coef.,p=Pr...z..)
}





fit <- survfit(Surv(Day, Status) ~ Condition, weights = Frequency,
               data = data.sel)

reg <- survreg(Surv(Day, Status) ~ 0+Condition, weights = Frequency,data = data.sel,dist="exponential")
reg

cox <- coxph(Surv(Day, Status) ~ Bacteria+Metformin_mM, weights = Frequency,data = data.sel)
cox


model.matrix(res.cox)


data("leuk", package = "MASS")

leuk.cox <- coxph(Surv(time) ~ ag + log(wbc), data = leuk)

library(multcomp)

lht <- glht(cox, linfct = diag(length(coef(cox))))

lht <- glht(cox)

summary(lht)




2.7015-3.1803

-0.4788








data.stat<-data %>%
  filter(Bacteria=='OP50') %>%
  group_by(Worm) %>%
  do(ph2frame(coxph(Surv(Day, Status) ~Metformin_mM, weights = Frequency,data = .) ))

data.stat<-data %>%
  filter(Metformin_mM=='0') %>%
  group_by(Worm) %>%
  do(ph2frame(coxph(Surv(Day, Status) ~Bacteria, weights = Frequency,data = .) ))

data.stat

install.packages('coxme')

library(coxme)


#library(multcomp)
#glht(res.cox,linfct=mcp(Condition="N2"))


fit2frame<-function(fit,tmin=0,tmax){
  data.frame('Time'=fit$time, 
                  'Survival'=fit$surv,
                  'Upper'=fit$upper,
                  'Lower'=fit$lower,
                  'At_risk'=fit$n.risk,
                  'Dying'=fit$n.event,
                  'SE'=fit$std.err) %>%
    add_row(Time=tmin,Survival=1) %>%
    add_row(Time=tmax,Survival=0) %>%
    mutate_at(group_vars(data), first)
}






data.freq<-data %>%
  group_by(Condition,Worm,Bacteria,Metformin_mM) %>%
  do(fit2frame(survfit(Surv(Day, Status) ~ Condition, weights = Frequency,
             data = .),0,40))



ggplot(data.freq,aes(x=Time,y=Survival,color=Bacteria,fill=Bacteria,linetype=Metformin_mM))+
  #geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.3,color=NA)+
  geom_step(direction='hv')+
  scale_x_continuous(breaks=seq(0,40,by=10),limits=c(0,40))+
  facet_wrap(~Worm)

dev.copy2pdf(device=cairo_pdf,
             file="Lifespans.pdf",
             width=12,height=6, useDingbats=FALSE)





data.freq %>%
  filter(Worm=='nhr-1' & Metformin_mM=='0')



data %>%
  filter(Worm=='nhr-1' & Metformin_mM=='0')





res<-summary(fit)






fit <- survfit(Surv(Day, Status) ~ Bacteria+Metformin_mM, weights = Frequency,
               data = data.sel)


ggsurvplot(fit,data = data.sel,risk.table = TRUE,
           ggtheme = theme_minimal())





############################################################################
# apply methods in Val's JCO using the supplmentary example
############################################################################

library(tidyverse)
library(readxl)
library(janitor)
library(gtsummary)
library(nnet)
library(survival)
library(broom)

############################################################################
#read in dataset

df<-read_xlsx("dat.xlsx") %>%
  clean_names() %>%
  mutate(domain =factor(domain))

glimpse(df)

############################################################################
#A: fit polytomous logistic regression with nnet:multinom() in JUST domains 1-3

df.mcidonly<-df %>%
  filter(domain!=0) %>% ##remove those with no MCID
  mutate(domain.D3ref = case_when(domain==3~0, ##recode D3 to the reference (0)
                                  domain==1~1, 
                                  domain==2~2) %>%as.factor())
#checks
glimpse(df)
glimpse(df.mcidonly)

#fit, using D3 as reference
mod1<-multinom(domain.D3ref~treat+log(time_to_mcid), 
               data=df.mcidonly)

#look at model 
summary(mod1)
exp(coef(mod1))
exp(confint(mod1))
##not identical to paper!

tbl_regression(mod1, exponentiate=TRUE)

############################################################################
#B: get probabilities


df.mcidonly %>%
  mutate(prob = fitted(mod1)) ->prob

#split by trt
prob %>%
  filter(treat==0) %>%
  group_by(time_to_mcid) %>%
  slice_head() ->p.trt0

prob %>%
  filter(treat==1) %>%
  group_by(time_to_mcid) %>%
  slice_head()->p.trt1

p.trt0
p.trt1

##also slightly different...

############################################################################
#C: cox model to fit the model on everyone:

df %>%
  mutate(event = case_when(domain==0~0, 
                           domain==1~1, 
                           domain==2~1, 
                           domain==3~1)) ->df.cox #recode event for cox

coxph(Surv(time_to_mcid, event)~ treat, 
        ties='breslow', #why not efron? esp given ties?
        data=df.cox) ->cox.mod

trt_HR<-exp(cox.mod$coefficients)

############################################################################
#D: get risk table to calc hazard
### why isnt this split by treatment?? -> assuming null is true

fit<-survfit(Surv(time_to_mcid, event)~ 1, 
        data=df.cox)

out_tidy <- tidy(fit) %>%
  mutate(hazard = n.event/n.risk)
out_tidy

############################################################################
#E: stick all together, perform calcs, by trt

out_tidy %>%
  select(time, n.risk,n.event, hazard) %>%
  add_column( D1_prob=p.trt0$prob[,2], 
              D2_prob=p.trt0$prob[,3],
              D3_prob=p.trt0$prob[,1]) %>%
  mutate(haz_D1 = hazard*D1_prob, 
         haz_D2 = hazard*D2_prob, 
         haz_D3 = hazard*D3_prob, 
         cum_D1 = cumsum(haz_D1), 
         cum_D2 = cumsum(haz_D2),
         cum_D3 = cumsum(haz_D3),
         cumprob_D1 = 1-exp(-cum_D1), 
         cumprob_D2 = 1-exp(-cum_D2), 
         cumprob_D3 = 1-exp(-cum_D3), 
         trt=0)->final_trt0

out_tidy %>%
  select(time, n.risk,n.event, hazard) %>%
  add_column( D1_prob=p.trt1$prob[,2], 
              D2_prob=p.trt1$prob[,3],
              D3_prob=p.trt1$prob[,1]) %>%
  mutate(haz_D1 = hazard*D1_prob*trt_HR, 
         haz_D2 = hazard*D2_prob*trt_HR, 
         haz_D3 = hazard*D3_prob*trt_HR, 
         cum_D1 = cumsum(haz_D1), 
         cum_D2 = cumsum(haz_D2),
         cum_D3 = cumsum(haz_D3),
         cumprob_D1 = 1-exp(-cum_D1), 
         cumprob_D2 = 1-exp(-cum_D2), 
         cumprob_D3 = 1-exp(-cum_D3), 
         trt=1)->final_trt1

final_trt0 %>%
  rbind(final_trt1) %>%
  select(time, trt, contains("cumprob"))->all

flextable::flextable(all)
ggplot(all, aes(x = time, y = cumprob_D1, group = as.factor(trt), colour=as.factor(trt))) + 
  geom_step() +
  theme_minimal()# geom_step creates a stairs plot

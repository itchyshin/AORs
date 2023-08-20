library(here)
library(asreml) 

dat<-readRDS(here("../katana_data/dat_MA.RDS"))

# head(dat)
# dim(dat) # 16652995
# sum(complete.cases(dat)) # 16652995
# length(unique(dat$SAMPLING_EVENT_IDENTIFIER)) # 16652995
# length(unique(dat$COUNTRY)) 
# length(unique(dat$STATE_CODE)) 

# correlation

mod <- asreml(Zr1 ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
              weights = weights, family = asr_gaussian(dispersion = 1),
              data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

mod1 <- asreml(Zr1 ~ I(sqrt(VZr)) + I(log(effort_time)- 5.480639), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE)  + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat, workspace = WS, na.action = na.method(x = 'omit'))

mod2 <- asreml(Zr1 ~ VZr + I(log(effort_time)- 5.480639), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))


mod3 <- asreml(Zr1 ~ I(sqrt(VZr)) + scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

mod4 <- asreml(Zr1 ~ VZr +  scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))


mod5 <- asreml(Zr1 ~ I(sqrt(VZr)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

mod6 <- asreml(Zr1 ~ scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

mod7 <- asreml(Zr1 ~ VZr, 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

# skewness

mod8 <- asreml(skewness ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
               weights = Wskew, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))


mod9 <- asreml(skewness ~ I(sqrt(Vskew)) + scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = Wskew, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

mod10 <- asreml(skewness ~ Vskew +  scale(log(effort_time)), 
                random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
                weights = Wskew, family = asr_gaussian(dispersion = 1),
                data = dat,  workspace = WS, na.action = na.method(x = 'omit'))


mod11 <- asreml(skewness ~ I(sqrt(Vskew)), 
                random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
                weights = Wskew, family = asr_gaussian(dispersion = 1),
                data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

mod12 <- asreml(skewness ~ scale(log(effort_time)), 
                random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
                weights = Wskew, family = asr_gaussian(dispersion = 1),
                data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

mod13 <- asreml(skewness ~ Vskew, 
                random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
                weights = Wskew, family = asr_gaussian(dispersion = 1),
                data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

# please save all

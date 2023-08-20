# doing some job - measuring how long a model take to run 

# meta-analysis

# packages
library(here)
#library(metafor) # this takes way too long
#library(mixmeta) # this takes too long
library(tidyverse)
#library(performance)
library(asreml) # this is what we use

# functions

# typical measurement error variance
tmev <- function(VZr, k){
  mev <- sum(1/VZr) * (k - 1)/(sum(1/VZr)^2 - sum((1/VZr)^2))
  return(mev)
}

# I2 
# for mod
I2 <- function(mod){
  sigma2_v <- sum(mod$mf$weights) * (mod$noeff[["units"]] - 1)/(sum(mod$mf$weights)^2 - 
                                                                  sum((mod$mf$weights)^2))
  I2_total <- (sum(mod$vparameters)-1)/((sum(mod$vparameters)-1) + sigma2_v)
  #I2_each <- mod$vparameters/((sum(mod$vparameters)-1) + sigma2_v)
  #names(I2_each) <- paste0("I2_", model$s.names)
  #I2s <- c(I2_total = I2_total, I2_each)
  return(I2_total)
}

# mod8
I2_2 <- function(mod){
  sigma2_v <- sum(mod$mf$Wskew) * (mod$noeff[["units"]] - 1)/(sum(mod$mf$Wskew)^2 - 
                                                                  sum((mod$mf$Wskew)^2))
  I2_total <- (sum(mod$vparameters)-1)/((sum(mod$vparameters)-1) + sigma2_v)
  #I2_each <- mod$vparameters/((sum(mod$vparameters)-1) + sigma2_v)
  #names(I2_each) <- paste0("I2_", model$s.names)
  #I2s <- c(I2_total = I2_total, I2_each)
  return(I2_total)
}

# R2 to be written

# needs to be written for each model
# mod1 + mod3
# 
R2_1 <- function(mod){
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(mod$coefficients$fixed[c(2,1),]) %*% t(as.matrix(mod$mf[, c(4, 5)]))))
  
  # marginal
  R2m <- fix / (fix + (sum(mod$vparameters)-1))
  
  # conditional
  R2c <- (fix + (sum(mod$vparameters)-1) - mod$vparameters[(length(mod$vparameters)-1)]) /(fix + (sum(mod$vparameters)-1))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  return(R2s)
}

# mod2 + mod4
R2_2 <- function(mod){
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(mod$coefficients$fixed[c(2,1),]) %*% t(as.matrix(mod$mf[, c(2, 5)]))))
  
  # marginal
  R2m <- fix / (fix + (sum(mod$vparameters)-1))
  
  # conditional
  R2c <- (fix + (sum(mod$vparameters)-1) - mod$vparameters[(length(mod$vparameters)-1)]) /(fix + (sum(mod$vparameters)-1))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  return(R2s)
}

# mod5 + mod6
R2_3 <- function(mod){
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(mod$coefficients$fixed[c(1),]) %*% t(as.matrix(mod$mf[, c(4)]))))
  
  # marginal
  R2m <- fix / (fix + (sum(mod$vparameters)-1))
  
  # conditional
  R2c <- (fix + (sum(mod$vparameters)-1) - mod$vparameters[(length(mod$vparameters)-1)]) /(fix + (sum(mod$vparameters)-1))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  return(R2s)
}

# mod7
R2_4 <- function(mod){
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(mod$coefficients$fixed[c(1),]) %*% t(as.matrix(mod$mf[, c(2)]))))
  
  # marginal
  R2m <- fix / (fix + (sum(mod$vparameters)-1))
  
  # conditional
  R2c <- (fix + (sum(mod$vparameters)-1) - mod$vparameters[(length(mod$vparameters)-1)]) /(fix + (sum(mod$vparameters)-1))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  return(R2s)
}



# mod
system.time(mod <- readRDS(here("Rdata/mod.RDS")))
mod$noeff
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed
I2(mod)

rm(mod)

# mod1
#system.time(mod1 <- readRDS(here("Rdata/mod1.RDS")))
mod1$noeff
summary(mod1)$varcomp
summary(mod1, coef = TRUE)$coef.fixed
R2_1(mod1)

rm(mod1)

# mod2
system.time(mod2 <- readRDS(here("Rdata/mod2.RDS"))) 
mod2$noeff
summary(mod2)$varcomp
summary(mod2, coef = TRUE)$coef.fixed
R2_2(mod2)

rm(mod2)

# mod3
#system.time(mod3 <- readRDS(here("Rdata/mod3.RDS")))
mod3$noeff
summary(mod3)$varcomp
summary(mod3, coef = TRUE)$coef.fixed
R2_1(mod3)

rm(mod3)

# mod4
system.time(mod4 <- readRDS(here("Rdata/mod4.RDS"))) 
mod4$noeff
summary(mod4)$varcomp
summary(mod4, coef = TRUE)$coef.fixed
R2_2(mod4)

rm(mod4)

# mod5
#system.time(mod5 <- readRDS(here("Rdata/mod5.RDS"))) 
mod5$noeff
summary(mod5)$varcomp
summary(mod5, coef = TRUE)$coef.fixed
R2_3(mod5)

rm(mod5)

# mod6
system.time(mod6 <- readRDS(here("Rdata/mod6.RDS")))
mod6$noeff
summary(mod6)$varcomp
summary(mod6, coef = TRUE)$coef.fixed
R2_3(mod6)

rm(mod6)

# mod7
system.time(mod7 <- readRDS(here("Rdata/mod7.RDS"))) 
mod7$noeff
summary(mod7)$varcomp
summary(mod7, coef = TRUE)$coef.fixed
R2_4(mod7)

rm(mod7)
# correlation

# mod8
#system.time(mod8 <- readRDS(here("Rdata/mod8.RDS"))) 
mod8$noeff
summary(mod8)$varcomp
summary(mod8, coef = TRUE)$coef.fixed
I2_2(mod8)

rm(mod8)
# mod9
#system.time(mod9 <- readRDS(here("Rdata/mod9.RDS"))) 
mod9$noeff
summary(mod9)$varcomp
summary(mod9, coef = TRUE)$coef.fixed
R2_1(mod9)

rm(mod9)

# mod10
system.time(mod10 <- readRDS(here("Rdata/mod10.RDS"))) 
mod10$noeff
summary(mod10)$varcomp
summary(mod10, coef = TRUE)$coef.fixed
R2_2(mod10)

rm(mod10)

# mod11
system.time(mod11 <- readRDS(here("Rdata/mod11.RDS"))) 
mod11$noeff
summary(mod11)$varcomp
summary(mod11, coef = TRUE)$coef.fixed
R2_3(mod11)

rm(mod11)

# mod12
system.time(mod12 <- readRDS(here("Rdata/mod12.RDS"))) 
mod12$noeff
summary(mod12)$varcomp
summary(mod12, coef = TRUE)$coef.fixed
R2_3(mod12)

rm(mod12)

# mod13
system.time(mod13 <- readRDS(here("Rdata/mod13.RDS"))) 
mod13$noeff
summary(mod13)$varcomp
summary(mod13, coef = TRUE)$coef.fixed
R2_4(mod13)

rm(mod13)


#############
# Do not run
#############
# the model prepared for Szymek

# not run
meta_data<-readRDS(here("../katana_data/meta_data.RDS"))
dat_full <- distinct(meta_data)
mod3 <- asreml(Zr1 ~ I(sqrt(VZr)) + scale(log(effort_time)),
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat_full, workspace = WS, na.action = na.method(x = 'omit'))

saveRDS(mod3, "mod3.RDS")

###

library(here)
library(asreml) 

dat<-readRDS(here("../katana_data/dat_MA.RDS"))




# # takes ~1 min
meta_data<-readRDS(here("Rdata/meta_data.RDS"))
# 
dat <- distinct(meta_data)

#cases <- complete.cases(dat)

max(log(dat$effort_time))

#
tmev <- tmev(dat$VZr, dim(dat)[1])

# mod
system.time(mod <- readRDS(here("Rdata/mod.RDS")))
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed
I2(mod)


# mod1
system.time(mod1 <- readRDS(here("Rdata/mod1.RDS")))
summary(mod1)$varcomp
summary(mod1, coef = TRUE)$coef.fixed

# mod2
system.time(mod2 <- readRDS(here("Rdata/mod2.RDS"))) 
summary(mod2)$varcomp
summary(mod2, coef = TRUE)$coef.fixed


# mod3
system.time(mod3 <- readRDS(here("Rdata/mod3.RDS")))
summary(mod3)$varcomp
summary(mod3, coef = TRUE)$coef.fixed

# mod4
system.time(mod4 <- readRDS(here("Rdata/mod4.RDS"))) 
summary(mod4)$varcomp
summary(mod4, coef = TRUE)$coef.fixed

# mod5
system.time(mod5 <- readRDS(here("Rdata/mod5.RDS"))) 
summary(mod5)$varcomp
summary(mod5, coef = TRUE)$coef.fixed


# mod6
system.time(mod6 <- readRDS(here("Rdata/mod6.RDS")))
summary(mod6)$varcomp
summary(mod6, coef = TRUE)$coef.fixed

# mod7
system.time(mod7 <- readRDS(here("Rdata/mod7.RDS"))) 
summary(mod7)$varcomp
summary(mod7, coef = TRUE)$coef.fixed

################
# skew
#################

mod <- asreml(Zr1 ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
              weights = weights, family = asr_gaussian(dispersion = 1),
              data = dat)
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed
I2(mod)



# visualization

# funnel
funnel <- ggplot(meta_data, aes(y = sqrt(1/VZr), x = Zr1)) +
  geom_point(color="black", fill="gray70") +
  geom_vline(xintercept= 0, linetype="dashed", col = "red", size = 1.5) +
  #geom_smooth(method="lm") +
  theme_bw() +
  xlab("Zr (Effect size)")+
  ylab("Precision (1/SE)")

funnel

# bubble

#TODO - ....
bubble1 <- ggplot(dat, aes(y = Zr1, x = log(effort_time))) +
  geom_point(aes(size = sqrt(1/VZr))) +
  geom_hline(yintercept= 0, linetype="dashed", col = "red", size = 0.8) +
  geom_smooth(method="lm") + 
  theme_bw() +
  xlab("ln(Effort time)") +
  ylab("Zr (Effect size)") +
  labs(size = "Precision (1/SE)") +
  theme(legend.position= "bottom")
  
bubble1

bubble2 <- ggplot(dat, aes(y = Zr1, x = sqrt(VZr))) +
  geom_point(aes(size = sqrt(1/VZr))) +
  geom_hline(yintercept= 0, linetype="dashed", col = "red", size = 0.8) +
  geom_smooth(method="lm") + 
  theme_bw() +
  xlab("Standard error (SE)") +
  ylab("Zr (Effect size)") +
  labs(size = "Precision (1/SE)") +
  theme(legend.position= "bottom")

bubble2

# models
# head(dat)
# dim(dat)
# sum(complete.cases(dat))
# length(unique(dat$SAMPLING_EVENT_IDENTIFIER))

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







# not to run

##################
# this is the size of RAM - if not 10gb - it won't work I think we need at least 5gb for the matrix
asreml.options(workspace = "20gb")
# save design matrix
asreml.options(design=TRUE)


###############
#asreml-R ---- 
###############

mod <- asreml(Zr1 ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
              weights = weights, family = asr_gaussian(dispersion = 1),
              data = dat)
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed
I2(mod)

mod1 <- asreml(Zr1 ~ I(sqrt(VZr)) + I(log(effort_time)- 5.480639), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE)  + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)
summary(mod1)$varcomp
summary(mod1, coef = TRUE)$coef.fixed

mod2 <- asreml(Zr1 ~ VZr + I(log(effort_time)- 5.480639), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)
summary(mod2)$varcomp
summary(mod2, coef = TRUE)$coef.fixed


# two more models

mod3 <- asreml(Zr1 ~ I(sqrt(VZr)) + scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)

mod4 <- asreml(Zr1 ~ VZr +  scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)


mod5 <- asreml(Zr1 ~ I(sqrt(VZr)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)

mod6 <- asreml(Zr1 ~ scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)

mod7 <- asreml(Zr1 ~ VZr, 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)

# skewness

mod8 <- asreml(skewness ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
              weights = Wskew, family = asr_gaussian(dispersion = 1),
              data = dat)


mod9 <- asreml(skewness ~ I(sqrt(Vskew)) + scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = Wskew, family = asr_gaussian(dispersion = 1),
               data = dat)

mod10 <- asreml(skewness ~ Vskew +  scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)


mod11 <- asreml(skewness ~ I(sqrt(Vskew)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)

mod12 <- asreml(skewness ~ scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)

mod13 <- asreml(skewness ~ Vskew, 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat)


saveRDS(mod3, here("Rdata/mod3.RDS"))
saveRDS(mod4, here("Rdata/mod4.RDS"))
saveRDS(mod5, here("Rdata/mod5.RDS"))
saveRDS(mod6, here("Rdata/mod6.RDS"))
saveRDS(mod7, here("Rdata/mod7.RDS"))

# we should save 3 models
saveRDS(mod, here("Rdata/mod.RDS"))
saveRDS(mod1, here("Rdata/mod1.RDS"))
saveRDS(mod2, here("Rdata/mod2.RDS"))

# not to run

########
# brms
########

brm1 <- brm(Zr1 | se(sqrt(VZr)) ~ 1 +(1 |STATE_CODE)  + (1|COUNTRY) 
                        data = dat,
                        chains = 4, iter = 4000,
                        cores = 4,
                        backend = "cmdstanr")

brm2 <- brm(Zr1 | se(sqrt(VZr)) ~ 
              I(sqrt(VZr)) + I(log(effort_time)- 5.480639) +(1 |STATE_CODE)  + (1|COUNTRY) 
            data = dat,
            chains = 4, iter = 4000,
            cores = 4,
            backend = "cmdstanr")

brm3 <- brm(Zr1 | se(sqrt(VZr)) ~ 
              VZr + I(log(effort_time)- 5.480639) +(1 |STATE_CODE)  + (1|COUNTRY) 
            data = dat,
            chains = 4, iter = 4000,
            cores = 4,
            backend = "cmdstanr")


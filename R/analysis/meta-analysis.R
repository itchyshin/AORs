# meta-analysis script

#################
# packages & data
#################

library(ggplot2)
library(here)
library(tidyverse)
library(patchwork)
#library(asreml) # one needs to buy this package

dat <- readRDS(here("Rdata/dat_MA.RDS"))

###################
# custom functions
##################

# typical measurement error variance (sampling error variance)
tmev <- function(VZr){
  k <- length(VZr)
  mev <- sum(1/VZr) * (k - 1)/(sum(1/VZr)^2 - sum((1/VZr)^2))
  return(mev)
}

# I2 
I2 <- function(mod){
  sigma2_v <- sum(mod$mf$weights) * (mod$noeff[["units"]] - 1)/(sum(mod$mf$weights)^2 - 
                                                                  sum((mod$mf$weights)^2))
  I2_total <- (sum(mod$vparameters)-1)/((sum(mod$vparameters)-1) + sigma2_v)
  return(I2_total)
}

# R2 function

R2 <- function(mod){
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(mod$vcoeff$fixed) %*% t(as.matrix(mod$design[, 1:length(mod$vcoeff$fixed)]))))
  
  # marginal
  R2m <- fix / (fix + (sum(mod$vparameters)-1))
  
  # conditional
  R2c <- (fix + (sum(mod$vparameters)-1) - mod$vparameters[(length(mod$vparameters)-1)]) /(fix + (sum(mod$vparameters)-1))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  return(R2s)
}


###############################################
# running models (will not run without asreml)
###############################################
# Do not run 

# meta-analysis
ma <- asreml(Zr1 ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
              weights = weights, family = asr_gaussian(dispersion = 1),
              data = dat,  workspace = WS, na.action = na.method(x = 'omit'))



mr1 <- asreml(Zr1 ~ scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

mr2 <- asreml(Zr1 ~ VZr, 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units,
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

mr3 <- asreml(Zr1 ~ VZr +  scale(log(effort_time)), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
               weights = weights, family = asr_gaussian(dispersion = 1),
               data = dat,  workspace = WS, na.action = na.method(x = 'omit'))

##############################
# loading results from asreml
##############################

# typical sampling variance
tmev <- tmev(dat$VZr, dim(dat)[1])

# meta-analysis
system.time(ma <- readRDS(here("Rdata/ma.RDS")))
summary(ma)$varcomp
summary(ma, coef = TRUE)$coef.fixed
I2(ma)


# meta-regression: effort time
system.time(mr1 <- readRDS(here("Rdata/mr1.RDS")))
summary(mr1)$varcomp
summary(mr1, coef = TRUE)$coef.fixed

# meta-regression: sampling variance
system.time(mr2 <- readRDS(here("Rdata/mr2.RDS"))) 
summary(mr2)$varcomp
summary(mr2, coef = TRUE)$coef.fixed

# meta-regression: sampling variance + effort time
system.time(mr3 <- readRDS(here("Rdata/mr3.RDS"))) 
summary(mr3)$varcomp
summary(mr3, coef = TRUE)$coef.fixed


############
# Figure 2
###########

# funnel
funnel1 <- ggplot(dat, 
                  aes(y = sqrt(1/VZr), x = Zr1)) +
  geom_point(color="gray30", alpha = 0.5) +
  geom_vline(xintercept= 0, linetype="dashed", 
             color = "red", alpha = 0.8, linewidth = 1) +
  #xlim(-5, 5) + 
  theme_bw() +
  xlab("Zr (effect size)")+
  ylab("Precision (1/SE)")

funnel2 <- ggplot(dat, aes(y = N-3, x = ln_estimate)) +
  geom_point(color="gray30", alpha = 0.5) +
  geom_vline(xintercept= 0, linetype="dashed", 
             color = "red", alpha = 0.8, linewidth = 1) +
  xlim(-1, 1) + 
  theme_bw() +
  xlab("Correlation")+
  ylab("No. of species - 3 (1/sampling varaince)")

patchwork <- funnel1 / funnel2  + 
  #plot_layout(widths = c(2, 1)) + 
  plot_annotation(tag_levels = 'A')

#patchwork

ggsave("fig2.png", width = 12, height = 20, units = "cm")
unlink("fig2.png")



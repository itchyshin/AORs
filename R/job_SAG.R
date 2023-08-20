# doing some job - measuring how long a model take to run 

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/GRAL/Bank_LMM/Queries/Nakagawa")


# meta-analysis

# packages
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# library(gridExtra)
# library(lme4)
# library(purrr)
# library(readr)
library(here)
#library(metafor) # this takes way too long
#library(mixmeta) # this takes too long
library(tidyverse)
#library(performance)
library(asreml) # this is what we use

# try with a lot more data

# # takes ~1 min
meta_data<-readRDS(here("meta_data.RDS"))
# 
dat <- distinct(meta_data)

dim(meta_data)
head(meta_data)
str(meta_data)

summary(meta_data$weights)
sum(is.na(meta_data$weights)) # no missing values
sum(is.na(meta_data$Zr1))    # no missing values


##################
# this is the size of RAM - if not 10gb - it won't work I think we need at least 5gb for the matrix
#asreml.options(workspace = "20gb")

# save design matrix
# Omited for now
# asreml.options(design=TRUE)


#################
# Making sure the definitions are correct
# And checking some of the data for integrity

meta_data$COUNTRY <- as.factor(meta_data$COUNTRY)
meta_data$STATE_CODE <- as.factor(meta_data$STATE_CODE)
length(levels(meta_data$COUNTRY))
length(levels(meta_data$STATE_CODE))
head(levels(meta_data$STATE_CODE))

freq.c <- as.matrix(table(meta_data$COUNTRY))
head(freq.c)
summary(freq.c)
hist(freq.c)
sum(freq.c[,1]<=1)
#View(freq.c)

freq.s <- as.matrix(table(meta_data$STATE_CODE))
head(freq.s)
summary(freq.s)
hist(freq.s)
sum(freq.s[,1]<=1)
#View(freq.c)

hist(meta_data$weights)
summary(meta_data$weights)
sum(meta_data$weights<=1) # 59,738/16,848,817 = 0.003545

# Niue might be a problem
View(meta_data[meta_data$COUNTRY == 'Niue',])
# Remove this country
S60481420


###############
#asreml-R ----
# SAG modified code
###############

asreml.options(ai.sing=FALSE)  # To have the varcomp data frame
asreml.options(fail='soft')
asreml.options(fail='hard')
asreml.options(uspd=TRUE)

asreml.options(ai.sing=TRUE)  # To have the varcomp data frame

# Testing with a subset of the data
t_meta_data <- meta_data[1:10000,]
t_meta_data <- t_meta_data[t_meta_data$SAMPLING_EVENT_IDENTIFIER != 'S60481420',]

mod <- asreml(100*Zr1 ~ 1, random = ~ COUNTRY + STATE_CODE + units, 
              weights = weights, 
              family = asr_gaussian(dispersion = 1),
              na.action=na.method(y="include",x="include"),
              workspace = 1e08,
              data = t_meta_data)
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed
#I2(mod)

# Model fitted using the sigma parameterization.
# ASReml 4.1.0 Tue Feb 22 13:11:36 2022
# LogLik        Sigma2     DF     wall    cpu
# 1     -118638.9           1.0   9999 13:11:36    0.0
# 2 singularities in the Average Information matrix.
# Error in asreml(100 * Zr1 ~ 1, random = ~COUNTRY + STATE_CODE + units,  : 
#                   Singularities in the average information matrix.

# Forced with asreml.options(ai.sing=TRUE)
# 
#            component std.error  z.ratio bound %ch
# COUNTRY     43.63436        NA       NA     S   0
# STATE_CODE  43.63436        NA       NA     S   0
# units      872.55243  12.34167 70.69972     P   0
# units!R      1.00000        NA       NA     F   0

# Note COUNTRY is Singular, that means it is confounded with STATE_COUNTRY.
# So there is some competition between COUNTRY and STATE_CODE as these are
# not crossed term

# Testing the following much simpler model....

mod <- asreml(fixed=100*Zr1 ~ 1, random = ~ idv(COUNTRY,init=600), # + STATE_CODE + units, 
              weights = weights, 
              family = asr_gaussian(dispersion = 1),
              workspace = 1e08,
              data = t_meta_data)
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed
#I2(mod)

##############
# I stopped here, not clear to me what is the issue with the data that
# gives me the singularity.
# It might be good to try when weights are all > 1.
#
# The singularity needs to be resolved first. Not sure the motive,
# but something is affecting it.


###############
#asreml-R ---- 
###############


mod <- asreml(Zr1 ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, 
              weights = weights, family = asr_gaussian(dispersion = 1),
              data = dat)
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed
#I2(mod)

mod1 <- asreml(Zr1 ~ I(sqrt(VZr)) + I(log(effort_time)- 5.480639), 
               random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY)  + units,
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

# we should save 3 models
saveRDS(mod, here("Rdata/mod.RDS"))
saveRDS(mod1, here("Rdata/mod1.RDS"))
saveRDS(mod2, here("Rdata/mod2.RDS"))


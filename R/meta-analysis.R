# meta-analysis

# packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(lme4)
library(purrr)
library(readr)
library(here)
#library(metafor) # this takes way too long
library(mixmeta) # this takes too long
library(tidyverse)
#library(performance)
#library(asreml) # this is what we use

library(patchwork)

#install.packages("data.table") 
#install.packages("ggplot2") 
#install.packages("jsonlite")

# custum functions

# typical measurement error variance
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
  #I2_each <- mod$vparameters/((sum(mod$vparameters)-1) + sigma2_v)
  #names(I2_each) <- paste0("I2_", model$s.names)
  #I2s <- c(I2_total = I2_total, I2_each)
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

########################
# 27 Apr 2023

dat <- readRDS(here("Rdata/dat_MA.RDS"))
#dat_full <- readRDS(here("Rdata/data_full.RDS"))

dim(dat)
names(dat)

##names(dat_full)

#tmev(dat$VZr)
# [1] 0.07265144

# TODO - try with 10,000 data points first and see how it looks like
# good for testing
#dat_10 <- dat_full[1:10000, ]


##################
# Figure 2
#################

# funnel
funnel1 <- ggplot(dat_full, 
                  aes(y = sqrt(1/VZr), x = Zr1)) +
  geom_point(color="gray30", alpha = 0.5) +
  geom_vline(xintercept= 0, linetype="dashed", 
             color = "red", alpha = 0.8, linewidth = 1) +
  #xlim(-5, 5) + 
  theme_bw() +
  xlab("Zr (effect size)")+
  ylab("Precision (1/SE)")

funnel2 <- ggplot(dat_full, aes(y = N-3, x = ln_estimate)) +
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
#######################################
######################################
# try with a lot more data
setwd(here("Data/correlation_results_chunks/"))
# 
df <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)

setwd(here("abundance_occupancy"))

length(unique(df$SAMPLING_EVENT_IDENTIFIER))

# - check what is r[pearson] and r[spearman]
df

# person correlations
# p1 <- ggplot(df, aes(x=ln_estimate))+
#   geom_histogram(color="black", fill="gray70")+
#   theme_bw()+
#   xlab("Correlation")+
#   ylab("Number of samples")
# 
# # spearman correlations
# p2<- ggplot(df, aes(x=ra_estimate))+
#   geom_histogram(color="black", fill="gray70")+
#   theme_bw()+
#   xlab("Correlation")+
#   ylab("Number of samples")
# 
# sum(is.na(df$ln_estimate))
# mean(df$ln_estimate, na.rm=TRUE)
# weighted.mean(df$ln_estimate, df$parameter-2, na.rm = TRUE)

# df %>%
#   group_by(SR) %>%
#   summarize(mean_correlation=mean(estimate, na.rm=TRUE),
#             sd_correlation=sd(estimate, na.rm=TRUE)) %>%
#   ggplot(., aes(x=SR, y=mean_correlation))+
#   geom_jitter()+
#   theme_bw()+
#   xlab("Species Richness")+
#   ylab("Mean correlation")

# # meta-analytic data
meta_data <- df   %>%  mutate(N = ln_parameter + 2,
                              Zr1 = atanh(ln_estimate), # person's correlation of logged
                              Zr2 = atanh(ra_estimate), # spearman's correlation
                              VZr = 1/(N - 3), # sampling variance
                              Vskew = (6*N*(N - 1))/((N - 2)*(N + 1)*(N + 3)),  # variance for skew
                              )
# 
# # country has some missing values
# # | !is.na(COUNTY) should not be used - many countries do not have counties
# each of list has N = > 10 but they could end up still 2 or 3
meta_data <- meta_data %>% filter(N > 3, !is.na(Zr1)) %>%
  mutate(effort_distance_km = if_else(is.na(EFFORT_DISTANCE_KM), 0, EFFORT_DISTANCE_KM),
    time_per_km = DURATION_MINUTES/(effort_distance_km + 1), # note this is 
    effort_time = DURATION_MINUTES, # this will be better (effort_time) - 
    weights = 1/VZr,
    Wskew = 1/Vskew)

data_full <- distinct(meta_data)


data_full %>% select(COUNTRY, STATE_CODE, N, Zr1, VZr, skewness, Vskew, time_per_km, effort_time, weights, weights, Wskew) -> dat

dim(dat)
sum(complete.cases(dat))


saveRDS(dat, here("Rdata/dat_MA.RDS"))


dim(data_full)

summary(data_full)

unique(levels(factor(data_full$COUNTRY)))
unique(levels(factor(data_full$STATE_CODE)))

saveRDS(data_full, here("Rdata/data_full.RDS"))

# skew plotting

# funnel
funnel <- ggplot(data_full, aes(y = sqrt(1/Vskew), x = skewness)) +
  geom_point(color="black", fill="gray70") +
  geom_vline(xintercept= 0, linetype="dashed", col = "red", size = 1.5) +
  #geom_smooth(method="lm") +
  theme_bw() +
  xlab("Skewness (Effect size)")+
  ylab("Precision (1/SE)")

#funnel

# bubble


bubble1 <- ggplot(data_full, aes(y = skewness, x = log(effort_time))) +
  geom_point(aes(size = sqrt(1/Vskew))) +
  geom_hline(yintercept= 0, linetype="dashed", col = "red", size = 0.8) +
  geom_smooth(method="lm") + 
  theme_bw() +
  xlab("ln(Effort time)") +
  ylab("Skewness (Effect size)") +
  labs(size = "Precision (1/SE)") +
  theme(legend.position= "bottom")

#bubble1

bubble2 <- ggplot(data_full, aes(y = skewness, x = sqrt(Vskew))) +
  geom_point(aes(size = sqrt(1/Vskew))) +
  geom_hline(yintercept= 0, linetype="dashed", col = "red", size = 0.8) +
  geom_smooth(method="lm") + 
  theme_bw() +
  xlab("Standard error (SE)") +
  ylab("Skewness (Effect size)") +
  labs(size = "Precision (1/SE)") +
  theme(legend.position= "bottom")

#bubble2





###########################################

#saveRDS(meta_data, here("Rdata/meta_data.RDS"))


# old code - do not run

# takes ~1 min
system.time(meta_data <- readRDS(here("Rdata/meta_data.RDS")))

# NO. country 
length(unique(meta_data$COUNTRY)) # 12 counties - correct?? - fixed

# NO. state code
length(unique(meta_data$STATE_CODE)) # 192 - fixed?

dat <- distinct(meta_data)
# No of county
#length(unique(meta_data$COUNTY)) # 4134

#meta_data2 <- meta_data %>%  filter(Zr2 > -5)
#dim(meta_data2)

length(unique(dat$SAMPLING_EVENT_IDENTIFIER)) # 14638121
nrow(dat) # ????

# wait for it
p3 <- ggplot(meta_data, aes(y = sqrt(1/VZr), x = Zr1)) +
  geom_point(color="black", fill="gray70") +
  geom_vline(xintercept= 0, linetype="dotted", col = "red", size = 3) +
  #geom_smooth(method="lm") +
  theme_bw()+
  xlab("Zr (effect size)")+
  ylab("Precision (sqrt(1/VZr)")

p3

# wait for it - just fo UK
# meta_data %>% filter(COUNTRY == "United Kingdom") %>% 
#   ggplot(aes(y = sqrt(1/VZr), x = Zr1)) +
#   geom_point(color="black", fill="gray70") +
#   geom_vline(xintercept= 0, linetype="dotted", col = "red", size = 3) +
#   #geom_smooth(method="lm") +
#   theme_bw()+
#   xlab("Zr (effect size)")+
#   ylab("Precision (sqrt(1/VZr)") -> p4

##################
# this is the size of RAM - if not 10gb - it won't work I think we need at least 5gb for the matrix
asreml.options(workspace = "20gb")
# save design matrix
asreml.options(design = TRUE) # maybe we can suppress this??


# typical measurement error variance
tmev <- function(VZr){
  k <- length(VZr)
  mev <- sum(1/VZr) * (k - 1)/(sum(1/VZr)^2 - sum((1/VZr)^2))
  return(mev)
}

tmev(meta_data$VZr, dim(meta_data)[1])

# I2 
I2 <- function(mod){
  sigma2_v <- sum(mod$mf$weights) * (mod$noeff[["units"]] - 1)/(sum(mod$mf$weights)^2 - 
                                                                  sum((mod$mf$weights)^2))
  I2_total <- (sum(mod$vparameters)-1)/((sum(mod$vparameters)-1) + sigma2_v)
  #I2_each <- mod$vparameters/((sum(mod$vparameters)-1) + sigma2_v)
  #names(I2_each) <- paste0("I2_", model$s.names)
  #I2s <- c(I2_total = I2_total, I2_each)
  return(I2_total)
}

# R2 to be written

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


###############
#asreml-R ---- 
###############


##########
# test 1
###########

# TODO - figure out what was the point of this - getting countries out>>>
# how much time 
# 


mod <- asreml(Zr1 ~ 1, random = ~ factor(STATE_CODE) + units, weights = weights, family = asr_gaussian(dispersion = 1),
              data = meta_data)
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed

mod1 <- asreml(Zr1 ~ scale(VZr) + I(log(effort_time)- 5.5), random = ~ factor(STATE_CODE) + units, weights = weights, family = asr_gaussian(dispersion = 1),
               data = meta_data)
summary(mod1)$varcomp
summary(mod1, coef = TRUE)$coef.fixed

#mod1b <- asreml(Zr1 ~ scale(VZr) + I(log(effort_time)- 6), random = ~ factor(STATE_CODE) + units, weights = weights, family = asr_gaussian(dispersion = 1),
#                data = meta_data)
#summary(mod1b)$varcomp
#summary(mod1b, coef = TRUE)$coef.fixed


mod2 <- asreml(Zr1 ~ scale(VZr) + I(log(time_per_km)- 5.5), random = ~ factor(STATE_CODE) + units, weights = weights, family = asr_gaussian(dispersion = 1),
               data = meta_data)
summary(mod2)$varcomp
summary(mod2, coef = TRUE)$coef.fixed



# TODO - do not use county

# Pearson
mod <- asreml(Zr1 ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + units, weights = weights, family = asr_gaussian(dispersion = 1),
              data = meta_data)
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed
#I2(mod)

mod1 <- asreml(Zr1 ~ scale(I(sqrt(VZr))) + I(log(effort_time)- 5.480639), random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY)  + units, weights = weights, family = asr_gaussian(dispersion = 1),
               data = meta_data)
summary(mod1)$varcomp
summary(mod1, coef = TRUE)$coef.fixed

mod2 <- asreml(Zr1 ~ scale(VZr) + I(log(effort_time)- 5.480639), random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY)  + units, weights = weights, family = asr_gaussian(dispersion = 1),
               data = meta_data)
summary(mod2)$varcomp
summary(mod2, coef = TRUE)$coef.fixed

# we should save 3 models
saveRDS(mod, here("Rdata/mod.RDS"))
saveRDS(mod1, here("Rdata/mod1.RDS"))
saveRDS(mod2, here("Rdata/mod2.RDS"))



# mod
system.time(mod <- readRDS(here("Rdata/old/mod.RDS"))) # how long does it take? - 5.5 min
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed

# mod1
system.time(mod1 <- readRDS(here("Rdata/mod1.RDS"))) # how long does it take? - 5.5 min
summary(mod1)$varcomp
summary(mod1, coef = TRUE)$coef.fixed

# mod2
system.time(mod2 <- readRDS(here("Rdata/mod2.RDS"))) # how long does it take? - 5.5 min
summary(mod2)$varcomp
summary(mod2, coef = TRUE)$coef.fixed



###################
# rank correlation
##################

# model <- asreml(Zr2 ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY) + units, weights = weights, family = asr_gaussian(dispersion = 1),
#                 data = meta_data)
# summary(model)$varcomp
# summary(model, coef = TRUE)$coef.fixed
# I2(model)
# 
# model1 <- asreml(Zr2 ~ scale(VZr) + I(log(effort_time)- 5.480639), random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY)  + units, weights = weights, family = asr_gaussian(dispersion = 1),
#                  data = meta_data)
# summary(model1)$varcomp
# summary(model1, coef = TRUE)$coef.fixed
# 
# model2 <- asreml(Zr2 ~ scale(VZr) + I(log(time_per_km)- 5.480639), random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY)  + units, weights = weights, family = asr_gaussian(dispersion = 1),
#                  data = meta_data)
# summary(model2)$varcomp
# summary(model2, coef = TRUE)$coef.fixed

# model, model1 and model2 are models we may not need to worry about??

# we should save 3 models
# saveRDS(model, here("Data/model.RDS"))
# saveRDS(model1, here("Data/model1.RDS"))
# saveRDS(model2, here("Data/model2.RDS"))

##################
# Extra analysis
##################

meta_data_UK <- meta_data %>% filter(COUNTRY == "United Kingdom")

mod_UK <- asreml(Zr1 ~ 1, random = ~ factor(STATE_CODE) + factor(COUNTY) + units, weights = weights, family = asr_gaussian(dispersion = 1),
                 data = meta_data_UK)
summary(mod_UK)$varcomp
summary(mod_UK, coef = TRUE)$coef.fixed -> fixed
I2(mod_UK)
summary(mod_UK, coef = TRUE)$coef.random -> ranef
#fixed
#ranef
# 5.480639 = 4 hours
mod_UK1 <- asreml(Zr1 ~ scale(VZr) + I(log(time_per_km)- 5.480639), random = ~ factor(STATE_CODE) + factor(COUNTY) + units, weights = weights, family = asr_gaussian(dispersion = 1),
                  data = meta_data_UK)
summary(mod_UK1, coef = TRUE)$coef.fix
#mod_UK1$design
R2(mod_UK1)

summary(mod_UK1, coef = TRUE)$coef.fixed
summary(mod_UK1, coef = TRUE)$coef.random -> ranef
summary(mod_UK1)



##################################################################
##################################################################
##########
# MCMCglmm---- 
###########
Bmod <- MCMCglmm(Zr1 ~ 1, random = ~ factor(STATE_CODE), weights = weights, mev = VZr, prior = ,
                 data = meta_data)
summary(mod)
summary(mod$SOL)
summary(mod$VCV)
#############
# lmer ---- 
##############
#### lmer analyisis
model <- lmer(Zr1 ~1 + (1|COUNTY), weights = 1/VZr, data = meta_data)
summary(model)

hist(log(meta_data$effort_time))
max(log(meta_data$effort_time))
hist(log(meta_data$time_per_km))
max(log(meta_data$time_per_km))


model2 <- lmer(Zr ~1 + scale(VZr) + I(log(effort_time)- 6)+ (1|COUNTY), weights = 1/VZr, data = meta_data)
summary(model2)


model2b <- lmer(Zr ~1 + scale(VZr) + I(log(effort_time)- 5.5)+ (1|COUNTY), weights = 1/VZr, data = meta_data)
summary(model2b)

model3 <- lmer(Zr ~1 + scale(VZr) +  I(log(time_per_km)- 5.5)+ (1|COUNTY), weights = 1/VZr, data = meta_data)
summary(model3)

model4 <- lmer(Zr ~1 + scale(sqrt(VZr)) + I(log(effort_time)- 5.5)+ (1|COUNTY), weights = 1/VZr, data = meta_data)
summary(model4)


#model3 <- lmer(Zr ~1 + VZr + scale(log(effort_time))+ (1|COUNTY), weights = 1/VZr, data = meta_data)
#summary(model3)

ggplot(meta_data, aes(y = Zr, x = log(effort_time))) +
  geom_point(aes(size = 1/VZr)) +
  geom_smooth(method="lm")

r2_nakagawa(model1)

system.time(model <- mixmeta(Zr ~ 1, random = ~1| COUNTY,  S = VZr, data = meta_data, bscov="diag", method="reml"))
system.time(summary(model))


##########################
# Test attemps for 2 places 
#  doing a meta-analysis 
# mixmeta
##########################
meta_data <- readRDS("Data/correlation_results_chunks/CA-AB.RDS") 

meta_data %>% mutate(Zr = atanh(estimate),
                     VZr = 1/(parameter -2))-> meta_data

meta_data <- meta_data %>% filter(parameter > 2) %>% 
  mutate(EFFORT_DISTANCE_KM = if_else(is.na(EFFORT_DISTANCE_KM), 0, EFFORT_DISTANCE_KM),
         time_per_km = DURATION_MINUTES/(EFFORT_DISTANCE_KM+1),
         effort_time = DURATION_MINUTES)

hist(log(meta_data$effort_time))
hist(meta_data$VZr)
#min(meta_data$parameter)

system.time(model <- mixmeta(Zr ~ 1, S = VZr, data = meta_data, bscov="diag", method="reml"))
system.time(summary(model))

system.time(model1 <- mixmeta(Zr ~ 1 + I(log(effort_time)-6), S = VZr, data = meta_data, method="reml"))
summary(model1)

# different dataset 
# mixmeta ---- 
meta_data <- readRDS("Data/correlation_results_chunks/US-KS.RDS") 

meta_data %>% mutate(Zr = atanh(estimate),
                     VZr = 1/(parameter -2))-> meta_data

meta_data <- meta_data %>% filter(parameter > 2) %>% 
  mutate(#time_per_km = DURATION_MINUTES/EFFORT_DISTANCE_KM,
    effort_time = DURATION_MINUTES)

# hist(log(meta_data$time_per_km))
# hist(meta_data$VZr)
# min(meta_data$parameter)

system.time(model <- mixmeta(Zr ~ 1, S = VZr, data = meta_data, bscov="diag", method="reml"))
system.time(summary(model))

str(model)


system.time(model1 <- mixmeta(Zr ~ 1 + VZr + I(log(time_per_km)- 6) , S = VZr, data = meta_data, method="reml"))
summary(model1)

system.time(model2 <- mixmeta(Zr ~ 1 + VZr + I(log(effort_time)- 6) , S = VZr, data = meta_data, method="reml"))
summary(model2)

###############################
#
###############################

# this may not work. 

# the first lot of correlations (probably wrong)
meta_data <- readRDS("Data/test_meta_test_data.RDS") 

meta_data %>% mutate(Zr = atanh(correlation),
                     VZr = 1/(df -2))-> meta_data

meta_data <- meta_data[complete.cases(meta_data), ] %>% filter(df > 2) %>% mutate(time_per_km = DURATION_MINUTES/EFFORT_DISTANCE_KM)%>% filter(time_per_km < Inf)

hist(log(meta_data$time_per_km))
hist(meta_data$VZr)
min(meta_data$df)

system.time(model <- mixmeta(Zr ~ 1, S = VZr, data = meta_data, bscov="diag", method="reml"))
summary(model)

system.time(model1 <- mixmeta(Zr ~ 1 + VZr, S = VZr, data = meta_data, method="reml"))
summary(model1)

system.time(model2 <- mixmeta(Zr ~ 1 + VZr + I(DURATION_MINUTES- 300) , S = VZr, data = meta_data, method="reml"))
summary(model2)


system.time(model3 <- mixmeta(Zr ~ 1 + VZr + I(log(time_per_km)- 8) , S = VZr, data = meta_data, method="reml"))
summary(model3)

system.time(model4 <- mixmeta(Zr ~ 1 + I(log(time_per_km)- 8) , S = VZr, data = meta_data, method="reml"))
summary(model4)

plot(residuals(model4), 1/meta_data$VZr)
plot(Zr ~log(time_per_km), data = meta_data)

# COEFFICIENTS AND (CO)VARIANCE MATRIX
coef(model)
vcov(model)

# RESIDUALS AND FITTED VALUES
plot(residuals(model), 1/meta_data$VZr)
fitted(model)

#ma <- rma(yi = Zr, vi = VZr, data = meta_data)

summary(ma)

##############################
##############################
# Some data wrangling....
##############################
##############################

# How many species and bird observations
# 7,635 species &  3,020,551,089
#library(purrr)
#library(dplyr)

setwd("Data/data_summary/")

df <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)

setwd("../..")

table_s1 <- df %>%
  group_by(COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarize(number_of_times_in_a_model=sum(N),
            number_of_total_individuals=sum(total_individuals))

length(unique(table_s1$COMMON_NAME))
sum(table_s1$number_of_total_individuals)


##############
# putting data together 
#################

# try with a lot more data
setwd(here("Data/correlation_results_chunks/"))
# 
df <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)

setwd(here("abundance_occumancy"))

length(unique(df$SAMPLING_EVENT_IDENTIFIER))

# - check what is r[pearson] and r[spearman]
df

# person correlations
# p1 <- ggplot(df, aes(x=ln_estimate))+
#   geom_histogram(color="black", fill="gray70")+
#   theme_bw()+
#   xlab("Correlation")+
#   ylab("Number of samples")
# 
# # spearman correlations
# p2<- ggplot(df, aes(x=ra_estimate))+
#   geom_histogram(color="black", fill="gray70")+
#   theme_bw()+
#   xlab("Correlation")+
#   ylab("Number of samples")
# 
# sum(is.na(df$ln_estimate))
# mean(df$ln_estimate, na.rm=TRUE)
# weighted.mean(df$ln_estimate, df$parameter-2, na.rm = TRUE)

# df %>%
#   group_by(SR) %>%
#   summarize(mean_correlation=mean(estimate, na.rm=TRUE),
#             sd_correlation=sd(estimate, na.rm=TRUE)) %>%
#   ggplot(., aes(x=SR, y=mean_correlation))+
#   geom_jitter()+
#   theme_bw()+
#   xlab("Species Richness")+
#   ylab("Mean correlation")

# # meta-analytic data
df   %>%  mutate(Zr1 = atanh(ln_estimate),
                 Zr2 = atanh(ra_estimate),
                 VZr = 1/(ln_parameter -2))-> meta_data
# 
# # country has some missing values
# # | !is.na(COUNTY) should not be used - many countries do not have counties
meta_data <- meta_data %>% filter(ln_parameter > 2, !is.na(Zr1)) %>%
  mutate(#EFFORT_DISTANCE_KM = if_else(is.na(EFFORT_DISTANCE_KM), 0, EFFORT_DISTANCE_KM),
    #time_per_km = DURATION_MINUTES/(EFFORT_DISTANCE_KM+1),
    effort_time = DURATION_MINUTES, # this will be better (effort_time) - 
    weights = 1/VZr)

dim(meta_data)


##### brms model
library(brms)

system.time(brm1 <- brm(Zr1 | se(sqrt(VZr)) ~ 1 +(1 |STATE_CODE), # + (1|COUNTRY) 
               data = dat[1000000:2000000,],
               chains = 1, iter = 2000, warmup = 1000,
               cores = 1,
               backend = "cmdstanr"))
                   #                   weights = weights, family = asr_gaussian(dispersion = 1),
                   #                   data = dat[1000000:2000000,])


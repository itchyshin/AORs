# This test script is to perform to get the abundance-occupancy relationship

library(purrr)
library(dplyr)
library(here)

setwd(here("Data/data_summary/"))

df <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)

setwd("../..")

table_s1 <- df %>%
  group_by(COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarize(number_of_times_in_a_model=sum(N),
            number_of_total_individuals=sum(total_individuals))

length(unique(table_s1$COMMON_NAME)) # the number of species - 7635
length(unique(table_s1$SCIENTIFIC_NAME)) # the number of species - 7635
sum(table_s1$number_of_total_individuals) # the number of birds involved - 3005668285


# packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(lme4)
library(purrr)
library(readr)
library(here)
#library(metafor)
library(mixmeta)
library(tidyverse)
library(performance)
library(asreml)
library(patchwork)

# try with a lot more data

setwd("Data/correlation_results_chunks/")

df <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)

setwd("..")

length(unique(df$SAMPLING_EVENT_IDENTIFIER))

# log and raw data but we should use correlation of ln-ed data. 

# - check what is r[pearson] and r[spearman]
df

p1 <- ggplot(df, aes(x=ln_estimate))+
  geom_histogram(color="black", fill="gray70")+
  theme_bw()+
  xlab("Correlation (Pearson)")+
  ylab("Number of samples")

p2 <- ggplot(df, aes(x=ra_estimate))+
  geom_histogram(color="black", fill="gray70")+
  theme_bw()+
  xlab("Correlation")+
  ylab("Number of samples (Spearman on ln)")

p1 + p2


sum(is.na(df$ln_estimate))

mean(df$ln_estimate, na.rm=TRUE)

weighted.mean(df$ln_estimate, df$ln_parameter-2, na.rm = TRUE)

# df %>%
#   group_by(SR) %>%
#   summarize(mean_correlation=mean(estimate, na.rm=TRUE),
#             sd_correlation=sd(estimate, na.rm=TRUE)) %>%
#   ggplot(., aes(x=SR, y=mean_correlation))+
#   geom_jitter()+
#   theme_bw()+
#   xlab("Species Richness")+
#   ylab("Mean correlation")

# meta-analytic data
df   %>%  mutate(Zr1 = atanh(ln_estimate),
              Zr2 = atanh(ra_estimate),
              VZr = 1/(ln_parameter -2))-> meta_data

# country has some missing values
# | !is.na(COUNTY)
meta_data <- meta_data %>% filter(ln_parameter > 2, !is.na(Zr1), !is.na(COUNTY)) %>%
  mutate(EFFORT_DISTANCE_KM = if_else(is.na(EFFORT_DISTANCE_KM), 0, EFFORT_DISTANCE_KM),
         time_per_km = DURATION_MINUTES/(EFFORT_DISTANCE_KM+1),
         effort_time = DURATION_MINUTES,
         weights = 1/VZr)

dim(meta_data)

meta_data2 <- meta_data %>%  filter(Zr2 > -5)

dim(meta_data2)

# wait for it
p3 <- ggplot(meta_data, aes(y = sqrt(1/VZr), x = Zr2)) +
  geom_point(color="black", fill="gray70") +
  geom_vline(xintercept= 0, linetype="dotted", col = "red", size = 3) +
  #geom_smooth(method="lm") +
  theme_bw()+
  xlab("Zr (effect size)")+
  ylab("Precision (sqrt(1/VZr)")

p4 <- ggplot(meta_data, aes(y = sqrt(1/VZr), x = Zr1)) +
  geom_point(color="black", fill="gray70") +
  geom_vline(xintercept= 0, linetype="dotted", col = "red", size = 3) +
  #geom_smooth(method="lm") +
  theme_bw()+
  xlab("Zr (effect size)")+
  ylab("Precision (sqrt(1/VZr)")

p3 + p4


##################

asreml.options(workspace = "10gb")

# typical ,
tmev <- function(VZr, k){
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

###############
#asreml-R
###############

# rank correlation
model0 <- asreml(Zr2 ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY) + units, weights = weights, family = asr_gaussian(dispersion = 1),
                data = meta_data2)
summary(model0)$varcomp
summary(model0, coef = TRUE)$coef.fixed
I2(model0)



# rank correlation
model <- asreml(Zr2 ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY) + units, weights = weights, family = asr_gaussian(dispersion = 1),
                data = meta_data)
summary(model)$varcomp
summary(model, coef = TRUE)$coef.fixed
I2(model)

model1 <- asreml(Zr2 ~ scale(VZr) + I(log(effort_time)- 5.480639), random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY)  + units, weights = weights, family = asr_gaussian(dispersion = 1),
               data = meta_data)
summary(model1)$varcomp
summary(model1, coef = TRUE)$coef.fixed

model2 <- asreml(Zr2 ~ scale(VZr) + I(log(time_per_km)- 5.480639), random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY)  + units, weights = weights, family = asr_gaussian(dispersion = 1),
               data = meta_data)
summary(model2)$varcomp
summary(model2, coef = TRUE)$coef.fixed

# we should save 3 models
saveRDS(model, here("Data/model.RDS"))
saveRDS(model1, here("Data/model1.RDS"))
saveRDS(model2, here("Data/model2.RDS"))

mod2 <- readRDS(here("Data/mod1.RDS"))
summary(mod2, coef = TRUE)$coef.fixed



# Pearson
mod <- asreml(Zr1 ~ 1, random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY) + units, weights = weights, family = asr_gaussian(dispersion = 1),
              data = meta_data)
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed
I2(mod)

mod1 <- asreml(Zr1 ~ scale(VZr) + I(log(effort_time)- 5.480639), random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY)  + units, weights = weights, family = asr_gaussian(dispersion = 1),
               data = meta_data)
summary(mod1)$varcomp
summary(mod1, coef = TRUE)$coef.fixed

mod2 <- asreml(Zr1 ~ scale(VZr) + I(log(time_per_km)- 5.480639), random = ~ factor(COUNTRY) + factor(STATE_CODE) + factor(COUNTY)  + units, weights = weights, family = asr_gaussian(dispersion = 1),
               data = meta_data)
summary(mod2)$varcomp
summary(mod2, coef = TRUE)$coef.fixed

# we should save 3 models
saveRDS(mod, here("Data/mod.RDS"))
saveRDS(mod1, here("Data/mod1.RDS"))
saveRDS(mod2, here("Data/mod2.RDS"))

mod2 <- readRDS(here("Data/mod1.RDS"))
summary(mod2, coef = TRUE)$coef.fixed

##########
# test 1
mod <- asreml(Zr1 ~ 1, random = ~ factor(STATE_CODE) + units, weights = weights, family = asr_gaussian(dispersion = 1),
             data = meta_data)
summary(mod)$varcomp
summary(mod, coef = TRUE)$coef.fixed

mod1 <- asreml(Zr1 ~ scale(VZr) + I(log(effort_time)- 5.5), random = ~ factor(STATE_CODE) + units, weights = weights, family = asr_gaussian(dispersion = 1),
              data = meta_data)
summary(mod1)$varcomp
summary(mod1, coef = TRUE)$coef.fixed

mod1b <- asreml(Zr1 ~ scale(VZr) + I(log(effort_time)- 6), random = ~ factor(STATE_CODE) + units, weights = weights, family = asr_gaussian(dispersion = 1),
               data = meta_data)
summary(mod1b)$varcomp
summary(mod1b, coef = TRUE)$coef.fixed


mod2 <- asreml(Zr1 ~ scale(VZr) + I(log(time_per_km)- 5.5), random = ~ factor(STATE_CODE) + units, weights = weights, family = asr_gaussian(dispersion = 1),
               data = meta_data)
summary(mod2)$varcomp
summary(mod2, coef = TRUE)$coef.fixed

##########
# MCMCglmm
###########
Bmod <- MCMCglmm(Zr1 ~ 1, random = ~ factor(STATE_CODE), weights = weights, mev = VZr, prior = ,
              data = meta_data)
summary(mod)
summary(mod$SOL)
summary(mod$VCV)


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

#################
#################

grids <- st_read("Data/grid_5_degree.geojson") %>%
  rename(grid_id=ID)

# species with range size....
# 7773 species (575 grids)
density_dat <- readRDS("Data/grid_specific_density_and_abundance.RDS")

sum_dat <- density_dat %>% summarise(max_den = median[which.max(median)],
                     min_den = median[which.min(median)],
                     ave_den = exp(mean(log(median))),
                     prop_grids = length(unique(grid_id))/575,
                     logit_prop = qlogis(prop_grids),
                     area_sum = sum(area_km2))

spp_AOR_data_main <- density_dat %>% summarise(range_size = sum(area_km2*prop_covered),
                                          max_den = median[which.max(median)],
                                           min_den = median[which.min(median)],
                                           ave_den = exp(mean(log(median))),  # ingoring weigthing
                                           wave_den = exp(weighted.mean(log(median), area_km2*prop_covered)),
                                           prop_grids = length(unique(grid_id))/575,
                                           logit_prop = qlogis(prop_grids),
                                           max_abund = median[which.max(abundance)],
                                           min_abund = median[which.min(abundance)],
                                           ave_abund = exp(mean(log(abundance))),
                                           wave_abund = exp(weighted.mean(log(abundance), area_km2*prop_covered)))

#saveRDS(spp_AOR_data2, "Data/spp_AOR_data2.RDS")
saveRDS(spp_AOR_data_main, "spp_AOR_data_main.RDS")


ggplot(spp_AOR_data_main, aes(x=log(range_size), y = log10(wave_den)))+
  geom_point()+
  ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Mean density per grid") + xlab("Extent of occurance")

model0 <- lm(log10(wave_den)~ log(range_size) , data = spp_AOR_data_main)
summary(model0)
cor.test(log10(spp_AOR_data_main$wave_den), log(spp_AOR_data_main$range_size))


# Using all the species
# 7773 species (575 grids)
density_dat1 <- readRDS("Data/grid_specific_densities.RDS")

sum_dat0 <- density_dat1 %>% summarise(max_den = median[which.max(median)],
                                     min_den = median[which.min(median)],
                                     ave_den = exp(mean(log(median))), # igenoring weigting
                                     prop_grids = length(unique(grid_id))/579,
                                     logit_prop = qlogis(prop_grids))

spp_AOR_data1 <- sum_dat0

saveRDS(spp_AOR_data2, "Data/spp_AOR_data1.RDS")
########################
#########################
density_dat2 <- sum_dat %>% left_join(., density_dat)

# get data for lookup table
setwd("Data/species_grid_range_lookup_temp")
species_grid_range_lookup <- list.files(pattern = ".RDS") %>%
  map(readRDS) %>% 
  bind_rows()

#setwd("../..")

species_lookup_species_list <- species_grid_range_lookup %>%
  dplyr::select(Species) %>%
  distinct()

clements <- read_csv("Data/clements_clean.csv")

species_joined_1 <- species_lookup_species_list %>%
  left_join(., clements %>%
              group_by(TipLabel) %>%
              slice(1) %>%
              rename(Species=TipLabel) %>%
              dplyr::select(ebird_COMMON_NAME, Species)) %>%
  dplyr::filter(complete.cases(ebird_COMMON_NAME))

species_joined_2 <- species_lookup_species_list %>%
  dplyr::filter(! Species %in% species_joined_1$Species) %>%
  left_join(., clements %>%
              group_by(ebird_COMMON_NAME) %>%
              slice(1) %>%
              mutate(Species=gsub(" ", "_", ebird_SCIENTIFIC_NAME)) %>%
              dplyr::select(ebird_COMMON_NAME, Species)) %>%
  dplyr::filter(complete.cases(ebird_COMMON_NAME))

species_joined_final <- species_joined_1 %>%
  bind_rows(species_joined_2)

density_dat3 <- density_dat %>%
  left_join(., grids %>%
              st_set_geometry(NULL)) %>%
  left_join(., species_grid_range_lookup %>%
              left_join(., species_joined_final %>%
                          group_by(Species) %>%
                          slice(1))) %>%
  dplyr::filter(complete.cases(Species)) %>%
  mutate(abundance=(10^density)*(area_km2/2.59)*prop_covered) %>%
  dplyr::filter(!abundance==0)

sum_dat2 <- density_dat %>% summarise(max_abund = median[which.max(abundance)],
                                     min_abund = median[which.min(abundance)],
                                     ave_abund = exp(mean(log(abundance))),
                                     wave_abund = exp(weighted.mean(log(abundance), area_km2)),
                                     prop_grids = length(unique(grid_id))/575,
                                     logit_prop = qlogis(prop_grids),
                                     area_sum = sum(area_km2))

# rewriting from above
density_dat2 <- sum_dat2 %>% left_join(., density_dat)

# figure for the abundance occupancy relationship (using all the data)
p0 <- ggplot(density_dat2, aes(x=logit_prop, y = median))+
  geom_point()+
  scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Density per grid") + xlab("Proportion of occupancy (logit)")

p0


p0 <- ggplot(density_dat2, aes(x=log(area_sum), y = median))+
  geom_point()+
  scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Density per grid") + xlab("area_km2")

p0

# getting stats
# non-phylogenetic model
model0 <- lmer(log(median) ~ logit_prop + (1|ebird_COMMON_NAME), data = density_dat2)
summary(model0)

# non-phylogenetic model
model01 <- lmer(log(median) ~ log(area_sum) + (1|ebird_COMMON_NAME), data = density_dat2)
summary(model01)


############
# Density
# looking at the abundance occupancy relationship (data point per species)
p1 <- ggplot(sum_dat, aes(x=logit_prop, y = log10(wave_den)))+
  geom_point()+
  ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Mean density per grid") + xlab("Proportion of occupancy (logit)")

p2 <- ggplot(sum_dat, aes(x=logit_prop, y = log10(max_den)))+
  geom_point(colour = "red")+
  ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Maximum density per grid") + xlab("Proportion of occupancy (logit)")


p3 <- ggplot(sum_dat, aes(x=logit_prop, y = log10(min_den)))+
  geom_point(colour = "blue")+
  ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Minimum density per grid") + xlab("Proportion of occupancy (logit)")

grid.arrange(p1, p2, p3, nrow = 1)

# looking at the abundance occupancy relationship (data point per species)
p1 <- ggplot(sum_dat, aes(x=log(area_sum), y = log10(ave_den)))+
  geom_point()+
  ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Mean density per grid") + xlab("Extent of occurance")

p2 <- ggplot(sum_dat, aes(x=log(area_sum), y = log10(max_den)))+
  geom_point(colour = "red")+
  ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Maximum density per grid") + xlab("Extent of occurance")


p3 <- ggplot(sum_dat, aes(x=log(area_sum), y = log10(min_den)))+
  geom_point(colour = "blue")+
  ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Minimum density per grid") + xlab("Extent of occurance")

grid.arrange(p1, p2, p3, nrow = 1)

###########
# Abundance
##########
# looking at the abundance occupancy relationship (data point per species)
p1 <- ggplot(sum_dat2, aes(x=logit_prop, y = log10(ave_abund)))+
  geom_point()+
  #ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Abundance") + xlab("Proportion of occupancy (logit)")

p2 <- ggplot(sum_dat2, aes(x=logit_prop, y = log10(max_abund)))+
  geom_point(colour = "red")+
  #ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Abundance") + xlab("Proportion of occupancy (logit)")


p3 <- ggplot(sum_dat2, aes(x=logit_prop, y = log10(min_abund)))+
  geom_point(colour = "blue")+
  #ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Abundance") + xlab("Proportion of occupancy (logit)")

grid.arrange(p1, p2, p3, nrow = 1)

# looking at the abundance occupancy relationship (data point per species)
p1 <- ggplot(sum_dat2, aes(x=log(area_sum), y = log10(wave_abund)))+
  geom_point()+
  #ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Abundance") + xlab("Extent of occurance")

p2 <- ggplot(sum_dat2, aes(x=log(area_sum), y = log10(max_abund)))+
  geom_point(colour = "red")+
  #ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Abundance") + xlab("Extent of occurance")


p3 <- ggplot(sum_dat2, aes(x=log(area_sum), y = log10(min_abund)))+
  geom_point(colour = "blue")+
  #ylim(c(-7,4)) +
  #scale_y_log10()+
  geom_smooth(method="lm")+
  ylab("Abundance") + xlab("Extent of occurance")

grid.arrange(p1, p2, p3, nrow = 1)



# Modeling

# non-pylogenetic model
model1 <- lm(log(ave_den) ~ logit_prop , data = sum_dat)
summary(model1)
model2 <- lm(log(max_den) ~ logit_prop , data = sum_dat)
summary(model2)
model3 <- lm(log(min_den) ~ logit_prop , data = sum_dat)
summary(model3)

# this is something I will use - the first one is still slightly positive 
model1b <- lm(log(ave_den) ~ log(area_sum) , data = sum_dat)
summary(model1b) # sig but basically fat
model2b <- lm(log(max_den) ~ log(area_sum) , data = sum_dat)
summary(model2b) # very positive 
model3b <- lm(log(min_den) ~ log(area_sum) , data = sum_dat)
summary(model3b) # very negative


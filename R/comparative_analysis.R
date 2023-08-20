# comparative analysis
# packages
library(gridExtra)
library(lme4)
library(here)
#library(metafor) # this takes way too long
library(mixmeta) # this takes too long
library(tidyverse)
library(magrittr)
# library(sf)
# library(geosphere)
library(ape)
library(norm2)
library(phylolm)
#################
#################

# grids <- st_read(here("Data/grid_5_degree.geojson")) %>%
#   rename(grid_id=ID)

# species with range size....
# 7773 species (575 grids) 
# density in log10 scale (just need to remember this)

density_dat <- readRDS("Data/grid_specific_density_and_abundance.RDS")
tree <- read.tree(here("Data", "phylo_data", "phy.tre"))

ll<-readr::read_csv("Data/clements_clean.csv")
density_dat %>% left_join(ll)-> density_dat

density_dat$density
# density_data is already grouped by bird spp
# some species have NA so reduced to 
aor_dat <- density_dat %>% ungroup() %>% group_by(TipLabel) %>% 
                          summarise(range_size = sum(area_km2*prop_covered), #,
                                    max_den = density[which.max(density)], # max density per species
                                    max_den_se = se[which.max(density)],    
                                    min_den = density[which.min(density)],
                                    min_den_se = se[which.min(density)],
                                    ave_den = mean(density),  # average or mean density
                                    ave_den_se = sqrt(mean(se^2)), # uncertainty (not used in the analysis but unbiased)
                                    #wave_den = weighted.mean(density, area_km2*prop_covered),
                                    #wave_den_se = sqrt(weighted.mean(se^2, area_km2*prop_covered)),
                                    #prop_grids = length(unique(grid_id))/575, # this is not as good as range_size
                                    #logit_prop = qlogis(prop_grids),
                                     ) 

# not matching to the tree names (e.g. subspecies etc, genera)
aor_dat <- filter(aor_dat,!is.na(TipLabel))
dat <- data.frame(aor_dat)
row.names(dat) <- dat$TipLabel
dim(dat) # 7464 species

###################
# Will's analysis 
###################

res_list <- vector(mode = "list", length = length(tree))
#summary(lm(ave_den~log(area_found),data=dat))

# using species range
run_phylo_model <- function(i){
  onetree <- tree[[i]]
  #run phylo model here
  model.phylo <- phylolm(ave_den ~ log10(range_size),
                         data = dat,
                         phy = onetree,
                         model = "BM")
  
  #extract coefficient and SE
  smp<-summary(model.phylo)
  return(smp$coefficients[ ,1:2])
}

phylo_model_out <- map(1:100,run_phylo_model)
int_dat <- map_df(phylo_model_out, ~.x[1, ])
slope_dat <- map_df(phylo_model_out, ~.x[2, ])

#library(norm2)
int_res <- miInference(as.list(int_dat$Estimate),as.list(int_dat$StdErr))
slope_res <- miInference(as.list(slope_dat$Estimate),as.list(slope_dat$StdErr))

int_res
slope_res
#non-significant

# t value
# intercept
int_res$est/int_res$std.err
(1-pt(int_res$est/int_res$std.err, int_res$df))*2

# slope
slope_res$est/slope_res$std.err
(1-pt(slope_res$est/slope_res$std.err, slope_res$df))*2

# get CI
# intercept
int_res$est - qt(0.975, int_res$df)*int_res$std.err
int_res$est + qt(0.975, int_res$df)*int_res$std.err

# slope
slope_res$est - qt(0.975, slope_res$df)*slope_res$std.err
slope_res$est + qt(0.975, slope_res$df)*slope_res$std.err


######################
#
#####################

plot <- ggplot(dat, aes(x=log10(range_size), y = ave_den)) +
  geom_point(color="gray30", alpha = 0.5)+
  #ylim(c(-7,4)) +
  geom_abline(intercept = int_res$est, slope = slope_res$est, col = "blue") +
  #geom_smooth(method="lm")+ 
  ylab("Mean density per grid (log10)") + xlab("Species range size (log10)") +
  theme_bw()

plot # 550 x 450

#ggsave("fig3.png", width = 24, height = 10, units = "cm")
#unlink("fig3.png")

################################################
###############################################

# use for Kim's stuff
density_dat <- readRDS("Data/grid_specific_density_and_abundance.RDS")
tree <- read.tree(here("Data", "phylo_data", "phy.tre"))

ll<-readr::read_csv("Data/clements_clean.csv")
density_dat %>% left_join(ll)-> density_dat

#for(i in 1:num_tree){
trimed_tree <- drop.tip(tree[[1]], which((tree$tip.label %in% dat$Phylogeny) == FALSE))  
cor_tree <- vcv(trimed_tree,corr=T)  
#saveRDS(spp_AOR_data2, "Data/spp_AOR_data2.RDS")
saveRDS(aor_dat, "Data/aor_dat.RDS")

aor_dat <- readRDS(here("Data/aor_dat.RDS"))

model <- rma.mv(yi = ave_den, V = ave_den_se^2, 
       mod = ~ log(area_found),
       random = list(~1|Phylogeny, ~1|Species),
       R= list(phylogeny = cor_tree),
       test = "t",
       control = list(optimizer = "optim", optmethod="Nelder-Mead"),
       data = aor_dat)


############################
# Old
###########################################################################
##########################################################################
# this includes information for how many grides they are found
sum_dat <- density_dat %>% summarise(max_den = median[which.max(median)],
                                     min_den = median[which.min(median)],
                                     ave_den = exp(mean(log(median))),
                                     prop_grids = length(unique(grid_id))/575,
                                     logit_prop = qlogis(prop_grids),
                                     area_sum = sum(area_km2))

# older version
# we need to incorporate error into this
spp_AOR_data_main <- density_dat %>% summarise(range_size = sum(area_km2*prop_covered),
                                               max_den = median[which.max(median)],
                                               min_den = median[which.min(median)],
                                               ave_den = exp(mean(log(median))),  # ignoring weighing
                                               wave_den = exp(weighted.mean(log(median), area_km2*prop_covered)),
                                               prop_grids = length(unique(grid_id))/575,
                                               logit_prop = qlogis(prop_grids),
                                               max_abund = median[which.max(abundance)],
                                               min_abund = median[which.min(abundance)],
                                               ave_abund = exp(mean(log(abundance))),
                                               wave_abund = exp(weighted.mean(log(abundance), area_km2*prop_covered)))

#saveRDS(spp_AOR_data2, "Data/spp_AOR_data2.RDS")
saveRDS(spp_AOR_data_main, "Data/spp_AOR_data_main.RDS")

spp_AOR_data_main <- readRDS(here("Data/spp_AOR_data_main.RDS"))
#####################################################################
####################################################################





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

# get data for look up table
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

######





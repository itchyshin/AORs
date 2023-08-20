# preparing for comparative analysis dataset

library(ape)
library(tidyverse)
library(here)


density_dat <- readRDS("Data/grid_specific_density_and_abundance.RDS")

# tree data
tree <- read.tree(here("Data", "phylo_data", "phy.tre"))

ll<-readr::read_csv("Data/clements_clean.csv")
density_dat %>% left_join(ll)-> density_dat

density_dat$density
# density_data is already grouped by bird spp
# some species have NA so reduced to 
aor_dat <- density_dat %>% ungroup() %>% group_by(TipLabel) %>% 
  summarise(range_size = sum(area_km2*prop_covered), 
            ave_den = mean(density),  # average or mean density
            ave_den_se = sqrt(mean(se^2)), # uncertainty (not used in the analysis but unbiased)
  ) 

# not matching to the tree names (e.g. subspecies etc, genera)
aor_dat <- filter(aor_dat,!is.na(TipLabel))
comp_dat <- data.frame(aor_dat)
row.names(comp_dat) <- comp_dat$TipLabel
dim(comp_dat) # 7464 species

saveRDS(comp_dat, here("Data/comp_dat.RDS"))

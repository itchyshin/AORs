# comparative analysis

##################
# packages & data
##################

library(here)
library(tidyverse)
library(ape)
library(norm2)
library(phylolm)


dat <- readRDS(here("Data/comp_dat.RDS"))
# tree data
tree <- read.tree(here("Data", "phylo_data", "phy.tre"))

#####################################
# Phylogenetic comparative  analysis 
#####################################

res_list <- vector(mode = "list", length = length(tree))
#summary(lm(ave_den~log(area_found),data=dat))

# function to run phylogenetic comparative models
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

# running the model with 100 trees
phylo_model_out <- map(1:100,run_phylo_model)
int_dat <- map_df(phylo_model_out, ~.x[1, ])
slope_dat <- map_df(phylo_model_out, ~.x[2, ])

# combiniong 100 models
int_res <- miInference(as.list(int_dat$Estimate),as.list(int_dat$StdErr))
slope_res <- miInference(as.list(slope_dat$Estimate),as.list(slope_dat$StdErr))

# esimates and se
int_res
slope_res


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


##########
# Figure 3
##########

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

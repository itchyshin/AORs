# Zr hetero analysis

pacman::p_load(readxl, 
               readr, 
               metafor, 
               tidyverse, 
)

# reading files = Zr
Zr_csv <- list.files(path = "./Data/Zr", pattern = "*.csv", 
                     full.names = TRUE) %>% lapply(read_csv)

# get names of each .csv file
Zr_filenames <- list.files(path = "./Data/Zr", 
                           pattern = "*.csv", full.names = FALSE)
names(Zr_csv) <- Zr_filenames
Zr <- Zr_csv

# delete NAs, zero variance
for (i in 1:length(Zr)) {
  Zr[[i]] <- Zr[[i]][!is.na(Zr[[i]]$es) & !is.na(Zr[[i]]$var) & Zr[[i]]$var != 0 & !is.na(Zr[[i]]$year_pub), ]
}
#### delete +-Inf
for (i in 1:length(Zr)) {
  Zr[[i]] <- Zr[[i]] %>% na.omit()
}

# modeling
model_Zr <- NA
for (i in 1:length(Zr)) {
  model_Zr[i] <- rma.mv(yi = es, V = var, 
                        random = list(~1|study_ID/obs_ID), 
                        method = "REML", test = "t", data = Zr[[i]], 
                        sparse=TRUE, control=list(optimizer="optim")) %>% list()
}


length(model_Zr)


tau2 <- sapply(model_Zr, function(x) sum(x$sigma2))

# average tau2
mean(tau2)
median(tau2)
range(tau2)

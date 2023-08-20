library(tidyverse)


setwd("Data/data_summary/")

dat <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)

setwd("../..")

length(unique(dat$COMMON_NAME))

# try with a lot more data
setwd("Data/correlation_results_chunks/")
# 
df_all_correlations <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)

setwd("../..")

length(unique(df_all_correlations$SAMPLING_EVENT_IDENTIFIER))

length(unique(df_all_correlations$LOCALITY_ID))

sum(is.na(df_all_correlations$COUNTRY))
length(unique(df_all_correlations$COUNTRY))


sum(is.na(df_all_correlations$STATE_CODE))

length(unique(df_all_correlations$STATE_CODE))


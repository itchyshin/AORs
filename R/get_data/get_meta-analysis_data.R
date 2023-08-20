# creating data for meta-analysis from correlation files

library(tidyverse)

setwd(here("Data/correlation_results_chunks/"))

# combining all correlation
df <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)

dim(df)

# creating relevant variables
meta_data <- df   %>%  mutate(N = ln_parameter + 2,
                              Zr1 = atanh(ln_estimate), # person's correlation (log = to be used)
                              Zr2 = atanh(ra_estimate), # spearman's correlation (this is not quite correct)
                              VZr = 1/(N - 3), # sampling variance
)

# # | !is.na(COUNTY) should not be used - many countries do not have counties
# each of list has N = > 10 but they could end up still 2 or 3
meta_data <- meta_data %>% filter(N > 3, !is.na(Zr1)) %>%
  mutate(effort_distance_km = if_else(is.na(EFFORT_DISTANCE_KM), 0, EFFORT_DISTANCE_KM),
         time_per_km = DURATION_MINUTES/(effort_distance_km + 1), # note this is 
         effort_time = DURATION_MINUTES, # this will be better (effort_time) - 
         weights = 1/VZr)

# removing potential duplicates
data_full <- distinct(meta_data)

# getting data which is requried for meta-analysis
data_full %>% select(COUNTRY, STATE_CODE, N, Zr1, VZr, time_per_km, effort_time, weights) -> dat

# the final length
dim(dat)
sum(complete.cases(dat))

saveRDS(dat, here("Rdata/dat_MA.RDS"))
#saveRDS(data_full, here("Rdata/data_full.RDS"))

########
# extra
########

# the number of species and the number of observations

setwd(here("Data/data_summary/"))

df2 <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)

setwd("../..")

table_s1 <- df2 %>%
  group_by(COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarize(number_of_times_in_a_model=sum(N),
            number_of_total_individuals=sum(total_individuals))

length(unique(table_s1$COMMON_NAME))
sum(table_s1$number_of_total_individuals)

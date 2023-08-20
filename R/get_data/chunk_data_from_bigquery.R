# This is an R script to get data from bigquery
# into chunks on my local computer...

## packages
library(bigrquery)
library(dbplyr)
library(tidyr)
library(lubridate)
library(readr)
library(dplyr)
library(ggplot2)

# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa_may_2020')

# first I want to see how many observations
# see how many observations per state_code in bigquery
state_code_obs <- ebird %>%
  dplyr::select(STATE_CODE) %>%
  group_by(STATE_CODE) %>%
  summarize(N=n()) %>%
  collect(n=Inf)

state_code_obs %>%
  dplyr::filter(N>1000000) %>%
  nrow()

# looks like 63 state codes have > 1 million observations
# so will split these up into each of their own dataset
chunk_function_v1 <- function(state_code_name){
  
  message(paste0("Downloading data for ", state_code_name))
  
  ebird_data <- ebird %>%
    dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                  LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                  DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                  COUNTRY, STATE_CODE, COUNTY) %>%
    dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
    dplyr::filter(STATE_CODE == state_code_name) %>%
    collect(n=Inf)
  
  saveRDS(ebird_data, paste0("S:/CTC/ebird_may_2020_chunked/", state_code_name, ".RDS"))
  
}

# list of all state codes that > 1 million observations
state_codes <- state_code_obs %>%
  dplyr::filter(N>1000000) %>%
  arrange(N) %>%
  .$STATE_CODE

# apply the function
lapply(state_codes[1:53], chunk_function_v1)


# now I want to group state codes into sizes that make sense.
# and filter by these and save these as 'chunks' to our data_chunked
# folder on my hard drive
states <- state_code_obs %>%
  dplyr::filter(N<=1000000) %>%
  arrange(N) %>%
  mutate(cum=cumsum(N))

# now will manually create 'chunks'
# which will be a list of state_codes
# to extract from eBird
chunk_1 <- states %>%
  dplyr::filter(cum<1000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_1)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_1.RDS")


states_2 <- states %>%
  dplyr::filter(!STATE_CODE %in% chunk_1) %>%
  mutate(cum=cumsum(N))

chunk_2 <- states_2 %>%
  dplyr::filter(cum<2500000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_2)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_2.RDS")

states_3 <- states_2 %>%
  dplyr::filter(!STATE_CODE %in% chunk_2) %>%
  mutate(cum=cumsum(N))

chunk_3 <- states_3 %>%
  dplyr::filter(cum<2500000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_3)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_3.RDS")


states_4 <- states_3 %>%
  dplyr::filter(!STATE_CODE %in% chunk_3) %>%
  mutate(cum=cumsum(N))

chunk_4 <- states_4 %>%
  dplyr::filter(cum<2500000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_4)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_4.RDS")

states_5 <- states_4 %>%
  dplyr::filter(!STATE_CODE %in% chunk_4) %>%
  mutate(cum=cumsum(N))

chunk_5 <- states_5 %>%
  dplyr::filter(cum<2500000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_5)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_5.RDS")

states_6 <- states_5 %>%
  dplyr::filter(!STATE_CODE %in% chunk_5) %>%
  mutate(cum=cumsum(N))

chunk_6 <- states_6 %>%
  dplyr::filter(cum<2500000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_6)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_6.RDS")

states_7 <- states_6 %>%
  dplyr::filter(!STATE_CODE %in% chunk_6) %>%
  mutate(cum=cumsum(N))

chunk_7 <- states_7 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_7)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_7.RDS")

states_8 <- states_7 %>%
  dplyr::filter(!STATE_CODE %in% chunk_7) %>%
  mutate(cum=cumsum(N))

chunk_8 <- states_8 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_8)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_8.RDS")

states_9 <- states_8 %>%
  dplyr::filter(!STATE_CODE %in% chunk_8) %>%
  mutate(cum=cumsum(N))

chunk_9 <- states_9 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_9)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_9.RDS")

states_10 <- states_9 %>%
  dplyr::filter(!STATE_CODE %in% chunk_9) %>%
  mutate(cum=cumsum(N))

chunk_10 <- states_10 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_10)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_10.RDS")

states_11 <- states_10 %>%
  dplyr::filter(!STATE_CODE %in% chunk_10) %>%
  mutate(cum=cumsum(N))

chunk_11 <- states_11 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_11)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_11.RDS")

states_12 <- states_11 %>%
  dplyr::filter(!STATE_CODE %in% chunk_11) %>%
  mutate(cum=cumsum(N))

chunk_12 <- states_12 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_12)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_12.RDS")

states_13 <- states_12 %>%
  dplyr::filter(!STATE_CODE %in% chunk_12) %>%
  mutate(cum=cumsum(N))

chunk_13 <- states_13 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_13)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_13.RDS")

states_14 <- states_13 %>%
  dplyr::filter(!STATE_CODE %in% chunk_13) %>%
  mutate(cum=cumsum(N))

chunk_14 <- states_14 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_14)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_14.RDS")

states_15 <- states_14 %>%
  dplyr::filter(!STATE_CODE %in% chunk_14) %>%
  mutate(cum=cumsum(N))

chunk_15 <- states_15 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_15)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_15.RDS")

states_16 <- states_15 %>%
  dplyr::filter(!STATE_CODE %in% chunk_15) %>%
  mutate(cum=cumsum(N))

chunk_16 <- states_16 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_16)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_16.RDS")

states_17 <- states_16 %>%
  dplyr::filter(!STATE_CODE %in% chunk_16) %>%
  mutate(cum=cumsum(N))

chunk_17 <- states_17 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_17)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_17.RDS")

states_18 <- states_17 %>%
  dplyr::filter(!STATE_CODE %in% chunk_17) %>%
  mutate(cum=cumsum(N))

chunk_18 <- states_18 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_18)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_18.RDS")

states_19 <- states_18 %>%
  dplyr::filter(!STATE_CODE %in% chunk_18) %>%
  mutate(cum=cumsum(N))

chunk_19 <- states_19 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_19)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_19.RDS")

states_20 <- states_19 %>%
  dplyr::filter(!STATE_CODE %in% chunk_19) %>%
  mutate(cum=cumsum(N))

chunk_20 <- states_20 %>%
  dplyr::filter(cum<5000000) %>%
  .$STATE_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(chunk_20)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_20.RDS")

# from above I was able to write out state_codes for 1:53.
state_codes[54:63]


# now I want to manually chunk up the remaining ones by "COUNTY"
# first I will get the number of observations per state/county
# and then will implement the same approach used for states
# using a cumsum
state_county_obs <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(STATE_CODE %in% local(state_codes[54:63])) %>%
  group_by(COUNTY_CODE) %>%
  summarize(N=n()) %>%
  collect(n=Inf)

state_county_obs_cum <- state_county_obs %>%
  ungroup() %>%
  dplyr::filter(complete.cases(COUNTY_CODE)) %>%
  arrange(N) %>%
  mutate(cum=cumsum(N))

chunk_21_codes <- state_county_obs_cum %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE 

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_21_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_21.RDS")

county_codes_22 <- state_county_obs_cum %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_21_codes) %>%
  mutate(cum=cumsum(N))

chunk_22_codes <- county_codes_22 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_22_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_22.RDS")

county_codes_23 <- county_codes_22 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_22_codes) %>%
  mutate(cum=cumsum(N))

chunk_23_codes <- county_codes_23 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_23_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_23.RDS")

county_codes_24 <- county_codes_23 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_23_codes) %>%
  mutate(cum=cumsum(N))

chunk_24_codes <- county_codes_24 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_24_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_24.RDS")


county_codes_25 <- county_codes_24 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_24_codes) %>%
  mutate(cum=cumsum(N))

chunk_25_codes <- county_codes_25 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_25_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_25.RDS")

county_codes_26 <- county_codes_25 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_25_codes) %>%
  mutate(cum=cumsum(N))

chunk_26_codes <- county_codes_26 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_26_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_26.RDS")

county_codes_27 <- county_codes_26 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_26_codes) %>%
  mutate(cum=cumsum(N))

chunk_27_codes <- county_codes_27 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_27_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_27.RDS")

county_codes_28 <- county_codes_27 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_27_codes) %>%
  mutate(cum=cumsum(N))

chunk_28_codes <- county_codes_28 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_28_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_28.RDS")

county_codes_29 <- county_codes_28 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_28_codes) %>%
  mutate(cum=cumsum(N))

chunk_29_codes <- county_codes_29 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_29_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_29.RDS")

county_codes_30 <- county_codes_29 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_29_codes) %>%
  mutate(cum=cumsum(N))

chunk_30_codes <- county_codes_30 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_30_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_30.RDS")

county_codes_31 <- county_codes_30 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_30_codes) %>%
  mutate(cum=cumsum(N))

chunk_31_codes <- county_codes_31 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_31_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_31.RDS")

county_codes_32 <- county_codes_31 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_31_codes) %>%
  mutate(cum=cumsum(N))

chunk_32_codes <- county_codes_32 %>%
  dplyr::filter(cum<5000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_32_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_32.RDS")

county_codes_33 <- county_codes_32 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_32_codes) %>%
  mutate(cum=cumsum(N))

chunk_33_codes <- county_codes_33 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_33_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_33.RDS")

county_codes_34 <- county_codes_33 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_33_codes) %>%
  mutate(cum=cumsum(N))

chunk_34_codes <- county_codes_34 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_34_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_34.RDS")

county_codes_35 <- county_codes_34 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_34_codes) %>%
  mutate(cum=cumsum(N))

chunk_35_codes <- county_codes_35 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_35_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_35.RDS")

county_codes_36 <- county_codes_35 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_35_codes) %>%
  mutate(cum=cumsum(N))

chunk_36_codes <- county_codes_36 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_36_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_36.RDS")

county_codes_37 <- county_codes_36 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_36_codes) %>%
  mutate(cum=cumsum(N))

chunk_37_codes <- county_codes_37 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_37_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_37.RDS")

county_codes_38 <- county_codes_37 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_37_codes) %>%
  mutate(cum=cumsum(N))

chunk_38_codes <- county_codes_38 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_38_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_38.RDS")

county_codes_39 <- county_codes_38 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_38_codes) %>%
  mutate(cum=cumsum(N))

chunk_39_codes <- county_codes_39 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_39_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_39.RDS")

county_codes_40 <- county_codes_39 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_39_codes) %>%
  mutate(cum=cumsum(N))

chunk_40_codes <- county_codes_40 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_40_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_40.RDS")

county_codes_41 <- county_codes_40 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_40_codes) %>%
  mutate(cum=cumsum(N))

chunk_41_codes <- county_codes_41 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_41_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_41.RDS")

county_codes_42 <- county_codes_41 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_41_codes) %>%
  mutate(cum=cumsum(N))

chunk_42_codes <- county_codes_42 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_42_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_42.RDS")

county_codes_43 <- county_codes_42 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_42_codes) %>%
  mutate(cum=cumsum(N))

chunk_43_codes <- county_codes_43 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_43_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_43.RDS")

county_codes_44 <- county_codes_43 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_43_codes) %>%
  mutate(cum=cumsum(N))

chunk_44_codes <- county_codes_44 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_44_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_44.RDS")

county_codes_45 <- county_codes_44 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_44_codes) %>%
  mutate(cum=cumsum(N))

chunk_45_codes <- county_codes_45 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_45_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_45.RDS")

county_codes_46 <- county_codes_45 %>%
  dplyr::filter(!COUNTY_CODE %in% chunk_45_codes) %>%
  mutate(cum=cumsum(N))

chunk_46_codes <- county_codes_46 %>%
  dplyr::filter(cum<10000000) %>%
  .$COUNTY_CODE

ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, PROTOCOL_TYPE, CATEGORY, 
                DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE, COUNTY, COUNTY_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2005-01-01") %>%
  dplyr::filter(COUNTY_CODE %in% local(chunk_46_codes)) %>%
  collect(n=Inf)

saveRDS(ebird_data, "S:/CTC/ebird_may_2020_chunked/chunk_46.RDS")

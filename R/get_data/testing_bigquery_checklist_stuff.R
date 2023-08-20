
## packages
library(bigrquery)
library(dbplyr)
library(tidyr)
library(lubridate)
library(readr)
library(dplyr)
library(ggplot2)


clements <- read_csv("Data/clements_clean.csv") %>%
  group_by(ebird_COMMON_NAME) %>%
  slice(1)

range_sizes <- read_delim("Data/range_sizes/bird.range.season041219.csv", 
                                      ";", escape_double = FALSE, trim_ws = TRUE)

joined_1 <- clements %>%
  left_join(., range_sizes %>%
              rename(TipLabel=Species) %>%
              mutate(TipLabel=gsub(" ", "_", TipLabel))) %>%
  dplyr::filter(complete.cases(range.size.km2))



# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa_may_2020')


## extract data
dat <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2019-01-01") %>%
  dplyr::filter(OBSERVATION_DATE < "2019-12-31") %>%
  group_by(SAMPLING_EVENT_IDENTIFIER) %>%
  summarize(N=n()) %>%
  collect(n=Inf)

# get rid of checklists with < 5 species richness
dat2 <- dat %>%
  dplyr::filter(N>=10)

# sample 100000 checklists
checklists <- dat2 %>%
  sample_n(10000) %>%
  .$SAMPLING_EVENT_IDENTIFIER

# get random data from bigquery for a test case
all_ebird_data <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(OBSERVATION_DATE > "2019-01-01") %>%
  dplyr::filter(OBSERVATION_DATE < "2019-12-31") %>%
  dplyr::filter(SAMPLING_EVENT_IDENTIFIER %in% local(checklists)) %>%
  collect(n=Inf)

saveRDS(all_ebird_data, "Data/temp_play_data.RDS")

all_ebird_data <- readRDS("Data/temp_play_data.RDS")

get_correlation_function <- function(checklist_id){
  
  message(paste0("Analyzing checklist ID: ", checklist_id))
  
  checklist_dat <- all_ebird_data %>%
    dplyr::filter(SAMPLING_EVENT_IDENTIFIER==checklist_id) %>%
    left_join(., joined_1 %>%
                rename(COMMON_NAME=ebird_COMMON_NAME)) %>%
    mutate(OBSERVATION_COUNT=as.numeric(as.character(OBSERVATION_COUNT)))
  
  summary_df <- if (nrow(checklist_dat %>%
      dplyr::filter(complete.cases(range.size.km2))) > 4 & 
      nrow(checklist_dat %>%
           dplyr::filter(complete.cases(OBSERVATION_COUNT))) > 4) {
    
    cor_obj <- cor.test(checklist_dat$OBSERVATION_COUNT, checklist_dat$range.size.km2)
    cor_obj_v2 <- cor.test(log(checklist_dat$OBSERVATION_COUNT), log(checklist_dat$range.size.km2))
    
    summary_df <- checklist_dat %>%
      dplyr::select(SAMPLING_EVENT_IDENTIFIER, LATITUDE, LONGITUDE,
                    DURATION_MINUTES, EFFORT_DISTANCE_KM) %>%
      distinct() %>%
      mutate(total_checklist_richness=nrow(checklist_dat)) %>%
      mutate(correlation=cor_obj$estimate,
             df=cor_obj$parameter,
             ci_lwr=cor_obj$conf.int[1],
             ci_upr=cor_obj$conf.int[2],
             correlation_v2=cor_obj_v2$estimate,
             df_v2=cor_obj_v2$parameter,
             ci_lwr_v2=cor_obj_v2$conf.int[1],
             ci_upr_v2=cor_obj_v2$conf.int[2])
  } else {
    summary_df <- checklist_dat %>%
      dplyr::select(SAMPLING_EVENT_IDENTIFIER, LATITUDE, LONGITUDE,
                    DURATION_MINUTES, EFFORT_DISTANCE_KM) %>%
      distinct() %>%
      mutate(total_checklist_richness=nrow(checklist_dat)) %>%
      mutate(correlation=NA,
             df=NA,
             ci_lwr=NA,
             ci_upr=NA,
             correlation_v2=NA,
             df_v2=NA,
             ci_lwr_v2=NA,
             ci_upr_v2=NA)
  }
  
  return(summary_df)
}

lists <- all_ebird_data %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER) %>%
  distinct() %>%
  .$SAMPLING_EVENT_IDENTIFIER

data_test_1_list <- lapply(lists[1:4000], get_correlation_function)
data_test_1 <- bind_rows(data_test_1_list)
data_test_2_list <- lapply(lists[4001:8000], get_correlation_function)
data_test_2 <- bind_rows(data_test_2_list)
data_test_3_list <- lapply(lists[8001:10000], get_correlation_function)
data_test_3 <- bind_rows(data_test_3_list)

data_test <- data_test_1 %>%
  bind_rows(data_test_2) %>%
  bind_rows(data_test_3) %>%
  distinct()

ggplot(data_test, aes(x=correlation))+
  geom_histogram(color="black", fill="gray80")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Pearson Correlation")+
  ylab("Number of checklist samples")

saveRDS(data_test, "Data/test_meta_test_data.RDS")

# get the number of species richness for all checklists in the eBird dataset
dat <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  group_by(SAMPLING_EVENT_IDENTIFIER) %>%
  summarize(N=n()) %>%
  collect(n=Inf)

all_checklists_trimmed <- dat %>%
  dplyr::filter(N>=10)

saveRDS(all_checklists_trimmed, "Data/bigquery_stuff/checklists_with_more_than_10_species.RDS")
    
# see how many observations per state_code in bigquery
state_code_obs <- ebird %>%
  dplyr::select(STATE_CODE) %>%
  group_by(STATE_CODE) %>%
  summarize(N=n()) %>%
  collect(n=Inf)



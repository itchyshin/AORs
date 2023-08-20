# Get correlation data for ebird abundance and
# the range size of the species on the checklist
# this loops through data chunks stored as RDSs outside of this repository
# and then saves dataframes of the results

# packages
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(broom)
library(purrr)

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

range_data <- joined_1 %>%
  dplyr::select(ebird_COMMON_NAME, range.size.km2) %>%
  rename(COMMON_NAME=ebird_COMMON_NAME)

cor_fun <- function(df) cor.test(log(df$OBSERVATION_COUNT), log(df$range.size.km2), method = "pearson") %>% 
  tidy()

cor_fun2 <- function(df) cor.test(rank(df$OBSERVATION_COUNT, na.last="keep"), 
                                  rank(df$range.size.km2, na.last="keep"), method = "pearson") %>%
  tidy()

skew_fun <- function(df) e1071::skewness(log(df$OBSERVATION_COUNT))

# build one big function
# to apply to a dataset

get_correlation_function <- function(file_name){
  
  message(paste0("Analyzing file: ", file_name))
  
  # read in data
  all_ebird_data <- readRDS(paste0("S:/CTC/ebird_may_2020_chunked/", file_name))
  
  # select only checklists which
  # have greater than 10 species (this still get the number of species to be fewer than 3 as some do not have range data)
  # and which have abundance information for ALL species reported
  # on the checklist
  checklist_info <- all_ebird_data %>%
    dplyr::filter(CATEGORY=="species") %>%
    mutate(OBSERVATION_COUNT=as.numeric(as.character(OBSERVATION_COUNT))) %>%
    group_by(SAMPLING_EVENT_IDENTIFIER) %>%
    summarize(SR=n(),
              include_xs=sum(is.na(OBSERVATION_COUNT)))
  
  lists <- checklist_info %>%
    dplyr::filter(SR>=10) %>%
    dplyr::filter(include_xs==0) %>%
    .$SAMPLING_EVENT_IDENTIFIER
  
  
  all_ebird_data <- all_ebird_data %>%
    dplyr::filter(CATEGORY=="species") %>%
    mutate(OBSERVATION_COUNT=as.numeric(as.character(OBSERVATION_COUNT))) %>%
    dplyr::filter(SAMPLING_EVENT_IDENTIFIER %in% lists) %>%
    group_by(SAMPLING_EVENT_IDENTIFIER) %>%
    left_join(., range_data, by="COMMON_NAME")
  
  list_metadata <- all_ebird_data %>%
    group_by(SAMPLING_EVENT_IDENTIFIER) %>%
    summarize(SR=n(),
              missing_counts=sum(is.na(OBSERVATION_COUNT)),
              missing_range=sum(is.na(range.size.km2))) %>%
    mutate(model_size=SR-missing_range) %>%
    dplyr::filter(model_size>=4) %>%
    left_join(., all_ebird_data %>%
                ungroup() %>%
                dplyr::select(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID, LATITUDE, LONGITUDE,
                              PROTOCOL_TYPE, DURATION_MINUTES, OBSERVATION_DATE,
                              EFFORT_DISTANCE_KM, EFFORT_AREA_HA,
                              NUMBER_OBSERVERS, COUNTRY, STATE_CODE, COUNTY) %>%
                distinct(), by="SAMPLING_EVENT_IDENTIFIER")
  
  # get a data summary where I summarize the number of observations per each species
  # and the total number of individuals for each species across all checklists
  dat_summary <- all_ebird_data %>%
    dplyr::filter(complete.cases(range.size.km2)) %>%
    group_by(COMMON_NAME, SCIENTIFIC_NAME) %>%
    summarize(N=n(),
              total_individuals=sum(OBSERVATION_COUNT)) %>%
    mutate(chunk_name=gsub(".RDS", " ", file_name))
  
  saveRDS(dat_summary, paste0("Data/data_summary/", file_name))
  
  nested_data <- all_ebird_data %>%
    ungroup() %>%
    dplyr::filter(SAMPLING_EVENT_IDENTIFIER %in% list_metadata$SAMPLING_EVENT_IDENTIFIER) %>%
    group_by(SAMPLING_EVENT_IDENTIFIER) %>%
    nest() 
  
  cor_results <- nested_data %>%
    mutate(model=map(data, cor_fun))
  
  corr_pr <- cor_results %>%
    dplyr::select(-data) %>%
    unnest(cols=c(model)) %>%
    dplyr::select(-method, -alternative) %>%
    rename(ln_estimate=estimate,
           ln_statistic=statistic,
           ln_p.value=p.value,
           ln_parameter=parameter,
           ln_conf.low=conf.low,
           ln_conf.high=conf.high)
  
  cor_results2 <- nested_data %>%
    mutate(model=map(data, cor_fun2))

  corr_sp <- cor_results2 %>%
    dplyr::select(-data) %>%
    unnest(cols=c(model)) %>%
    dplyr::select(-method, -alternative) %>%
    rename(ra_estimate=estimate,
           ra_statistic=statistic,
           ra_p.value=p.value,
           ra_parameter=parameter,
           ra_conf.low=conf.low,
           ra_conf.high=conf.high)
  
  skewness_results <- nested_data %>%
    mutate(model=map(data, skew_fun))
  
  skew_out <- skewness_results %>%
    unnest(cols=c(model)) %>%
    dplyr::select(-data) %>%
    rename(skewness=model)
    
  
  corr_results <- corr_pr %>%
    left_join(., corr_sp) %>%
    left_join(., skew_out) %>%
    left_join(., list_metadata)
  
  saveRDS(corr_results, paste0("Data/correlation_results_chunks/", file_name))
  
}

files <- list.files("S:/CTC/ebird_may_2020_chunked/")

file_information <- as.data.frame(file.info(list.files("S:/CTC/ebird_may_2020_chunked/", full.names=TRUE))) %>%
  mutate(file_name=files) %>%
  arrange(size)

lapply(file_information$file_name, get_correlation_function)



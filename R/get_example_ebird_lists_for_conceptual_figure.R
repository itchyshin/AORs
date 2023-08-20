# get example eBird checklists
# for methods figure

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

file_name="S:/CTC/ebird_may_2020_chunked/chunk_14.RDS"


all_ebird_data <- readRDS("S:/CTC/ebird_may_2020_chunked/chunk_14.RDS")

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

#saveRDS(dat_summary, paste0("Data/data_summary/", file_name))

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

corr_results <- corr_pr %>%
  left_join(., corr_sp) %>%
  left_join(., list_metadata)

# get a few lists manually
lists <- c("S10217872", "S27954200", "S34411835")

list_dat <- all_ebird_data %>%
  dplyr::filter(SAMPLING_EVENT_IDENTIFIER %in% lists)

list_1 <- list_dat %>%
  dplyr::filter(SAMPLING_EVENT_IDENTIFIER=="S10217872") %>%
  ggplot(., aes(x=log(range.size.km2), y=log(OBSERVATION_COUNT)))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="blacK"))+
  theme(panel.grid=element_blank())+
  xlab("Range size")+
  ylab("Local abundance")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_blank())+
  theme(panel.border=element_rect(color="orange"))+
  theme(axis.title=element_text(size=7))

list_1

ggsave("fig/example_panel_for_fig1_1.png", width=1.4, height=1.2, units="in")

list_2 <- list_dat %>%
  dplyr::filter(SAMPLING_EVENT_IDENTIFIER=="S27954200") %>%
  ggplot(., aes(x=log(range.size.km2), y=log(OBSERVATION_COUNT)))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="blacK"))+
  theme(panel.grid=element_blank())+
  xlab("Range size")+
  ylab("Local abundance")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_blank())+
  theme(panel.border=element_rect(color="purple"))+
  theme(axis.title=element_text(size=7))

list_2

ggsave("fig/example_panel_for_fig1_2.png", width=1.4, height=1.2, units="in")

list_3 <- list_dat %>%
  dplyr::filter(SAMPLING_EVENT_IDENTIFIER=="S34411835") %>%
  ggplot(., aes(x=log(range.size.km2), y=log(OBSERVATION_COUNT)))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="blacK"))+
  theme(panel.grid=element_blank())+
  xlab("Range size")+
  ylab("Local abundance")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_blank())+
  theme(panel.border=element_rect(color="green"))+
  theme(axis.title=element_text(size=7))

list_3

ggsave("fig/example_panel_for_fig1_3.png", width=1.4, height=1.2, units="in")






# this script will get the unique locality_id
# and save them out as a shapefile

library(tidyverse)
library(sf)

setwd("Data/correlation_results_chunks/")

df <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)

setwd("../..")

unique_points <- df %>%
  ungroup() %>%
  dplyr::select(LOCALITY_ID, LATITUDE, LONGITUDE) %>%
  distinct()

points_shape <- unique_points %>%
  st_as_sf(coords=c("LONGITUDE", "LATITUDE"), crs=4326)


st_write(points_shape, "Data/unique_points_shapefile/unique_localities.shp")

# work on making a map
# get subset of points
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot()+
  geom_sf(data = world)

random_points <- df %>%
  ungroup() %>%
  dplyr::filter(COUNTRY=="United States") %>%
  sample_n(1000000)

plot_points <- df %>%
  dplyr::filter(COUNTRY != "United States") %>%
  bind_rows(random_points)

png("Figure_test.png", width=6, height=4.5, res=300, units="in")

ggplot()+
  geom_sf(data = world)+
  geom_point(data=plot_points, aes(x=LONGITUDE, y=LATITUDE))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  theme(axis.text=element_blank())+
  theme(axis.title=element_blank())+
  theme(axis.ticks=element_blank())

dev.off()

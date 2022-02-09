library(sf)
library(ggplot2)
library(ggforce)
# library(tmap)

# Import the data ---------------------------------------------------------

# read all files of the shapefile
stratum <- st_read("strata/all_strata.shp")

stratum.new <- select(stratum, -c ("DIV", "strtm_t")) %>% 
  # Validating geometries to make operations on them 
  st_make_valid() %>%  
  # grouping by the stratum id
  group_by(stratum) %>% 
  # Combining all polygons that have the same stratum ID
  dplyr::summarize(geometry = st_union(geometry))

# simple plotting
plot(stratum$geometry)

##Lets add the cluster IDs to the abundance data. Seems like it could be useful 

clust.ID <- select (clust.fish, c ("clusterID"))
clust.ID$taxa_name <- row.names(clust.ID)

fish.abun.clust <- merge(x = fish.abun.clean, y= clust.ID, by = "taxa_name", all.x = TRUE)

#keep only the columns of  interest 

fish.abun.clustCL <- select(fish.abun.clust, -c ("year_obs", "season", "vessel", "distance_km", "area_swept_km2"))

# group by year and stratum #

fish.abun.gr <- fish.abun.clustCL %>% 
  group_by(stratum, year_surv) %>% 
  dplyr::summarise(Unique_FE = n_distinct(clusterID))

#join stratum table and unique FE table 
#first I need to filter teh data because I don't have the infor for all stratum

stratum.FE <- merge(x = fish.abun.gr, y= stratum.new, by = "stratum", all.y = TRUE) %>% 
  sf::st_as_sf()

# MAP
map1 <- stratum.FE %>% 
  # filter(year_surv %in% c(2000:2003)) %>% 
  ggplot() +
  geom_sf(aes(fill = Unique_FE)) +
  ggtitle(label = "map shit") +
  scale_fill_viridis_c () +
  theme_light() +
  # facet_wrap(~year_surv)
  ggforce::facet_wrap_paginate(~year_surv,
                               nrow = 2, ncol = 2, page = 2)

ggsave("plot.png", map1, width =  10, height = 10)

library(sf)
library(ggplot2)
library(ggforce)

# library(tmap)

# Import the data ---------------------------------------------------------

# read all files of the shapefile
stratum.shpfile <- st_read("strata/all_strata.shp")

stratum.new <- select(stratum.shpfile, -c ("DIV", "strtm_t")) %>% 
  # Validating geometries to make operations on them 
  st_make_valid() %>%  
  # grouping by the stratum id
  group_by(stratum) %>% 
  # Combining all polygons that have the same stratum ID
  dplyr::summarize(geometry = st_union(geometry))

plot(stratum.new$geometry)
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
  filter(year_surv %in% c(1995)) %>% 
  ggplot() +
  geom_sf(aes(fill = Unique_FE)) +
  ggtitle(label = "Dispersal Richness") +
  scale_fill_viridis_c () +
  theme_light() +
  facet_wrap(~year_surv)
  #ggforce::facet_wrap_paginate(~year_surv,
                              # nrow = 2, ncol = 2, page = 2)

ggsave("plot1995.png", map1, width =  10, height = 10)


#Effective functional diversity with Shannon in clustered groups
effective.FD <-fish.abun.clust %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, vessel, trip, set, clusterID)%>%
  dplyr::summarize(group_biomass = sum(density_kgperkm2))%>%
  #calculate average biomass per stratum for each functional group
  group_by(stratum, year_surv, clusterID)%>%
  dplyr::summarize(group_biomass = mean(group_biomass))%>%
  #calculate the Hill number for each stratum in each year
  group_by(stratum, year_surv)%>%
  dplyr::summarize(effective_species = exp(diversity(group_biomass, "shannon")))


#Effective species Diversity with Shannnon on species                       
effective.species <-fish.abun.clean %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, vessel, trip, set, taxa_name)%>%
  dplyr::summarize(group_biomass = sum(density_kgperkm2))%>%
  #calculate average biomass per stratum for each functional group
  group_by(stratum, year_surv, taxa_name)%>%
  dplyr::summarize(group_biomass = mean(group_biomass))%>%
  #calculate the Hill number for each stratum in each year
  group_by(stratum, year_surv)%>%
  dplyr::summarize(effective_species = exp(diversity(group_biomass, "shannon")))


effective.FD.strata <- merge(x = effective.FD , y= stratum.new, by = "stratum", all.y = TRUE) %>% 
  sf::st_as_sf()

effective.species.strata <- merge(x = effective.species, y= stratum.new, by = "stratum", all.y = TRUE) %>% 
  sf::st_as_sf()

# MAP
map2 <- effective.species.strata %>% 
  filter(year_surv %in% c(1995,2000,2005,2010)) %>% 
  ggplot() +
  geom_sf(aes(fill = effective_species)) +
  ggtitle(label = "Effective Species Diversity") +
  scale_fill_viridis_c () +
  theme_light() +
  facet_wrap(~year_surv)
#ggforce::facet_wrap_paginate(~year_surv,
# nrow = 2, ncol = 2, page = 2)

ggsave("EffectSpecies1995-2010.png", map2, width =  10, height = 10)

# MAP
map3 <- effective.FD.strata %>% 
  filter(year_surv %in% c(1995,2000,2005,2010)) %>% 
  ggplot() +
  geom_sf(aes(fill = effective_species)) +
  ggtitle(label = "Effective Functional Diversity") +
  scale_fill_viridis_c () +
  theme_light() +
  facet_wrap(~year_surv)
#ggforce::facet_wrap_paginate(~year_surv,
# nrow = 2, ncol = 2, page = 2)

ggsave("EffectiveFD1995-2010.png", map3, width =  10, height = 10)

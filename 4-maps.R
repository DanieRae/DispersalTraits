library(sf)
library(ggplot2)
# library(tmap)

# Import the data ---------------------------------------------------------

# read all files of the shapefile
stratum <- st_read("strata/all_strata.shp")

stratum.new <- select(stratum, -c ("DIV", "strtm_t"))

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

strata.filter.fish <- stratum.new$stratum

fish.abun.filter <- fish.abun.gr %>% filter(stratum %in% strata.filter.fish)

fish.abun.filter.strat <- fish.abun.filter$stratum 

stratum.filter <- stratum.new %>% filter(stratum %in% fish.abun.filter.strat)
   
  
stratum.FE <- merge(x = fish.abun.filter, y= stratum.filter, by = "stratum", all.x = FALSE)

#some duplication was happening so i ran this to remove duplicated rows#

stratum.FE.nodup <- stratum.FE %>% distinct(stratum, year_surv, .keep_all = TRUE)

#MAP

map1 <- ggplot(data = stratum.FE.nodup) +
  geom_sf( data = stratum.FE.nodup$geometry)

map2 <- ggplot(data = stratum.FE.nodup) +
  geom_sf( data = stratum.FE.nodup$geometry, 
           aes(fill = stratum.FE.nodup$Unique_FE))+
  ggtitle(label = "map shit")+
  scale_fill_viridis_c ()



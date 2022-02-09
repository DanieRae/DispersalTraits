#Statistics and response variables#

#need to call in the abundance data



# ga.fish.traits.dist <- gawdis(fish.traits)

functional.diversity.FT <- dbFD(fish.traits.dist,
                                corr = "lingoes")

#With the clustering of the species here a a few responses that should be looked at#

#FOR how much over redundancy in each cluster is there? Does the pattern of over redundancy change spatially? How many FEs have just one species? Is there a clear pattern in what characteristics have only 1 species? What are the ecological services provided by FE groups with high FOR vs those with none? Do we see change over time (quick/slow/more or less diversity in different areas), larger home ranges, higher dispersal capabilitys in certian FEs?

#Include abundance distributions (biomass) across FE's, before and during collapse 

#How does the FE diversity compare to hill numbers locally (alpha and gamma) gradients of beta accross the space


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
                                                     f      
  
fish.abun.sum <- fish.abun.gr %>% dplyr::summarise(total =n(),
                                                   cluster = clusterID)


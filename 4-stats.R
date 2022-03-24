#Statistics and response variables#

#need to call in the abundance data



# ga.fish.traits.dist <- gawdis(fish.traits)

#Abundance per community
Abun.fish <-fish.abun.clean %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, vessel, trip, set, taxa_name)%>%
  dplyr::summarize(group_biomass = sum(density_kgperkm2))%>%
  #calculate average biomass per stratum for each functional group
  group_by(stratum, year_surv, taxa_name)%>%
  dplyr::summarize(group_biomass = mean(group_biomass))

Abun.fish.wide <- as.data.table (spread(Abun.fish, taxa_name, group_biomass, fill = 0))

Abun.fish.wide.1995 <- subset(Abun.fish.wide,year_surv=="1995")

Abun.fish.wide.1995 <- select(Abun.fish.wide.1995, -c("year_surv","stratum"))

Abun.fish.wide.1995 <- sort(names(Abun.fish.wide.1995, decreasing = TRUE))


Identity <- functcomp(fish.traits.60percent, Abun.fish.wide.1995 )

functional.diversity.FT <- dbFD(fish.traits.60percent, Abun.fish.wide.1995,
                                corr = "lingoes")

#With the clustering of the species here a a few responses that should be looked at#

#FOR how much over redundancy in each cluster is there? Does the pattern of over redundancy change spatially? How many FEs have just one species? Is there a clear pattern in what characteristics have only 1 species? What are the ecological services provided by FE groups with high FOR vs those with none? Do we see change over time (quick/slow/more or less diversity in different areas), larger home ranges, higher dispersal capabilitys in certian FEs?

#Include abundance distributions (biomass) across FE's, before and during collapse 

#How does the FE diversity compare to hill numbers locally (alpha and gamma) gradients of beta accross the space

#Effective functional diversity
effective.FD.simp <-fish.abun.clust %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
    group_by(stratum, year_surv, vessel, trip, set, clusterID)%>%
    dplyr::summarize(group_biomass = sum(density_kgperkm2))%>%
              #calculate average biomass per stratum for each functional group
    group_by(stratum, year_surv, clusterID)%>%
    dplyr::summarize(group_biomass = mean(group_biomass))%>%
              #calculate the Hill number for each stratum in each year
    group_by(stratum, year_surv)%>%
    dplyr::summarize(effective_species = exp(diversity(group_biomass, "simpson")))


#Effective species Diversity                        
effective.species.simp <-fish.abun.clean %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, vessel, trip, set, taxa_name)%>%
  dplyr::summarize(group_biomass = sum(density_kgperkm2))%>%
  #calculate average biomass per stratum for each functional group
  group_by(stratum, year_surv, taxa_name)%>%
  dplyr::summarize(group_biomass = mean(group_biomass))%>%
  #calculate the Hill number for each stratum in each year
  group_by(stratum, year_surv)%>%
  dplyr::summarize(effective_species = exp(diversity(group_biomass, "simpson")))


effective.FD.strata.simp <- merge(x = effective.FD.simp , y= stratum.new, by = "stratum", all.y = TRUE) %>% 
  sf::st_as_sf()

effective.species.strata.simp <- merge(x = effective.species.simp, y= stratum.new, by = "stratum", all.y = TRUE) %>% 
  sf::st_as_sf()

# MAP
map3 <- effective.species.strata.simp %>% 
  filter(year_surv %in% c(1995,2000,2005,2010)) %>% 
  ggplot() +
  geom_sf(aes(fill = effective_species)) +
  ggtitle(label = "Effective Species Diversity") +
  scale_fill_viridis_c () +
  theme_light() +
  facet_wrap(~year_surv)
#ggforce::facet_wrap_paginate(~year_surv,
# nrow = 2, ncol = 2, page = 2)

ggsave("EffectSpecies1995-2010.png", map1, width =  10, height = 10)

# MAP
map4 <- effective.FD.strata.simp %>% 
  filter(year_surv %in% c(1995,2000,2005,2010)) %>% 
  ggplot() +
  geom_sf(aes(fill = effective_species)) +
  ggtitle(label = "Effective Functional Diversity") +
  scale_fill_viridis_c () +
  theme_light() +
  facet_wrap(~year_surv)
#ggforce::facet_wrap_paginate(~year_surv,
# nrow = 2, ncol = 2, page = 2)

ggsave("EffectiveFD1995-2005-2015.png", map2, width =  10, height = 10)

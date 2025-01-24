#DATA FROM PAGE 2 AND PAGE 3#
#THIS PAGE RUNS ALL THE DIFFERENT MAPS USED TO VISUALIZE SPATIOTEMPOTAL PATTERNS#

# Install and load packages ----
# install.packages("sf")
# install.packages("ggplot2")
# install.packages("ggforce")
# install.packages("patchwork")
# install.packages("dplyr")
# install.packages("stringr")
# install.packages("cmocean")
# install.packages("tmap")

library(sf)
library(ggplot2)
library(ggforce)
library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)
library(cmocean)
library(tmap)
library(vegan)


#EFFECTIVE DISPERSAL DIVERSITY ####
#Effective functional diversity with Shannon in clustered groups
fish.abun.clust.gr <- fish.abun.clust %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, clusterID) %>%
  dplyr::summarize(group_biomass = sum(group_biomass))%>%
#calculate average biomass per stratum for each functional group
  group_by(stratum, year_surv, clusterID) %>%
  dplyr::summarize(group_biomass = mean(group_biomass)) %>%
  ungroup()


#calculate the Hill number for each stratum in each year
effective.dispersal <- fish.abun.clust.gr %>%
  group_by(stratum, year_surv) %>%
  dplyr::summarize(effective_dispersal = exp(diversity(group_biomass, "shannon")))


effective.dispersal.strata <-
  merge(x = effective.dispersal ,
        y = stratum.shpfile,
        by = "stratum",
        all.y = TRUE) %>%
  sf::st_as_sf()


#EFFECTIVE SPECIES DIVERSITY ####
#Effective species Diversity with Shannnon on species

effective.species <- fish.abun.clean  %>%
  #calculate the Hill number for each stratum in each year
  group_by(stratum, year_surv) %>%
  dplyr::summarize(effective_species = exp(diversity(group_biomass, "shannon")))

effective.species.strata <-
  merge(x = effective.species,
        y = stratum.shpfile,
        by = "stratum",
        all.y = TRUE) %>%
  sf::st_as_sf()


## MAPS - EFFECTIVE DIVERSITY -----

map2 <- effective.species.strata %>%
  filter(year_surv %in% c(1995, 2000, 2005, 2010, 2015)) %>%
  ggplot() +
  geom_sf(aes(fill = effective_species),
          color = NA) +
  #coord_sf (xlim = c(-61, -46), ylim = c(42.5, 56)) +
  labs() +
  scale_fill_viridis_c(option = "magma") +
  theme_light() +
  facet_wrap(~ year_surv, ncol = 5) +
  labs(fill = "Taxonomic Diversity", colour = "") + # legend titles
  theme(plot.title = element_text(lineheight = .8, size = 15, hjust = 0.2), # title
        axis.text.x = element_blank(), # remove x axis labels
        #axis.text.y = element_blank(), # remove y axis labels
        axis.ticks = element_blank(), # remove axis ticks
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # remove grid lines
        strip.text.x = element_text(size = 14,vjust = 0, face = "bold",
                                    family = "Arial Narrow"),
        legend.position = "bottom",
        legend.key.size = unit(0.8, 'cm'))
#ggforce::facet_wrap_paginate(~year_surv,
# nrow = 2, ncol = 2, page = 2)

# ggsave(here("analysis", "figures",
#             "Spatiotemporal Variation in Effective Species Diversity.png"),
#        map2,
#        width =  10,
#        height = 10)


map3 <- effective.dispersal.strata %>%
  filter(year_surv %in% c(1995, 2000, 2005,2010,2015)) %>%
  ggplot() +
  geom_sf(aes(fill = effective_dispersal),
          color = NA) + #This removed the boarders for the strata, can use size =0.01 to make boarders as samll as possible
  #ggtitle() +
  scale_fill_viridis_c() +
  theme_light() +
  facet_wrap(~ year_surv, ncol = 5) +
  labs(fill = "Dispersal Group Diversity", colour = "") + # legend titles
  theme(plot.title = element_text(lineheight = .8, size = 15, hjust = 0.2), # title
        axis.text.x = element_blank(), # remove x axis labels
        #axis.text.y = element_blank(), # remove y axis labels
        axis.ticks = element_blank(), # remove axis ticks
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # remove grid lines
        strip.text.x = element_text(size = 14,vjust = 0, face = "bold",
                                    family = "Arial Narrow"),
        legend.position = "bottom",
        legend.key.size = unit(0.8, 'cm'))

# ggsave(here("analysis", "figures",
#             "Spatiotemporal Variation in Effective Dispersal Diversity.png"),
#        map3,
#        width =  10,
#        height = 10)

map.diversity <- map2 / map3

# Uncomment to save figure
# ggsave(here("analysis", "figures", "Taxonomic-DispersalGroup.png"),
#        map.diversity,
#        width =  12,
#        height = 10)


#STRATUM DEPTH MAP----

stratum.depth.geom <-
  merge (x = stratum.shpfile, y = stratum.depth, by = "stratum") %>%
  sf::st_as_sf()

map.depth <- stratum.depth.geom %>%
  ggplot() +
  geom_sf(aes(fill = depth.ave),
          color = "dark grey",
          size = 0.1)+
  coord_sf (xlim = c(-59, -47))+
  scale_fill_cmocean(name = "deep")+
  theme_light() +
  labs(fill = "Depth (m)", colour = "") + # legend titles
  theme(plot.title = element_text(lineheight = .8, size = 15, hjust = 0.2), # title
        axis.text.x = element_text(color = "black", size = 12), # remove x axis labels
        axis.text.y = element_text(color = "black", size = 12), # remove y axis labels
        #axis.ticks = element_blank(), # remove axis ticks
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), # remove grid lines
        legend.position = "bottom",
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(face = "bold"))

map.depth
ggsave("StratumDepth.png", map.depth, width =  10, height = 10)



#END-----

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
library(stringr)
library(cmocean)
library(tmap)
library(vegan)
library(tidyr)


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
#ggsave("StratumDepth.png", map.depth, width =  10, height = 10)
 


#Dispersal groups with hill number ####

#Biomass of clustered dispersal groups run with Shannon diversity
fish.abun.clust.gr <- fish.abun.clust %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, clusterID) %>%
  dplyr::summarize(group_biomass = sum(group_biomass))
#calculate average biomass per stratum for each functional group
# group_by(stratum, year_surv, clusterID) %>%
# dplyr::summarize(group_biomass = mean(group_biomass)) %>%
# ungroup()


#calculate the Hill number for each stratum in each year
effective.dispersal <- fish.abun.clust.gr %>%
  group_by(stratum, year_surv) %>%
  dplyr::summarize(effective_species = exp(diversity(group_biomass, "shannon")))


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
  scale_fill_viridis_c (option = "magma") +
  theme_light() +
  facet_wrap(~ year_surv, ncol = 5)+
  labs(fill = "Taxonomic Diversity", colour = "") + # legend titles
  theme(plot.title = element_text(lineheight = .8, size = 15, hjust = 0.2), # title
        axis.text.x = element_blank(), # remove x axis labels
        #axis.text.y = element_blank(), # remove y axis labels  
        axis.ticks = element_blank(), # remove axis ticks
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # remove grid lines
        strip.text.x = element_text(size=14,vjust=0, face = "bold",
                                    family="Arial Narrow"),
        legend.position="bottom",
        legend.key.size = unit(0.8, 'cm'))
#ggforce::facet_wrap_paginate(~year_surv,
# nrow = 2, ncol = 2, page = 2)

# ggsave("Spatiotemporal Variation in Effective Species Diversity.png",
#        map2,
#        width =  10,
#        height = 10)


map3 <- effective.dispersal.strata %>%
  filter(year_surv %in% c(1995, 2000, 2005,2010,2015)) %>%
  ggplot() +
  geom_sf(aes(fill = effective_species), 
          color = NA) + #This removed the boarders for the strata, can use size =0.01 to make boarders as samll as possible
  #ggtitle() +
  scale_fill_viridis_c () +
  theme_light() +
  facet_wrap(~ year_surv, ncol = 5)+
  labs(fill = "Dispersal Group Diversity", colour = "") + # legend titles
  theme(plot.title = element_text(lineheight = .8, size = 15, hjust = 0.2), # title
        axis.text.x = element_blank(), # remove x axis labels
        #axis.text.y = element_blank(), # remove y axis labels  
        axis.ticks = element_blank(), # remove axis ticks
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # remove grid lines
        strip.text.x = element_text(size=14,vjust=0, face = "bold",
                                    family="Arial Narrow"),
        legend.position = "bottom",
        legend.key.size = unit(0.8, 'cm'))
# 
# ggsave("Spatiotemporal Variation in Effective Dispersal Diversity.png",
#        map3,
#        width =  10,
#        height = 10)

map.diversity <- map2 / map3

# ggsave("Taxonomic-DispersalGroup.png",
#        map.diversity,
#        width =  12,
#        height = 10)

# BETA DIVERSITY ----

##BETA DISPERSAL GROUPS ----
# We split our main dataframe into a list of each dataframe corresponding to
# each stratum

fish.abun.clust.gr.split <-
  fish.abun.clust.gr %>%
  # We exclude stratum for which there is only data for 1 year
  filter(!(
    stratum %in%
      c(
        501,
        502,
        503,
        504,
        505,
        506,
        507,
        508,
        509,
        510,
        511,
        512,
        513,
        514,
        515,
        516,
        517,
        518,
        519,
        710,
        915,
        917,
        918
      )
  )) %>%
  # Split things
  split.data.frame(.$stratum)

# Using lapply, we apply a function to each data frame in the list. Here we
# create a presence absence dataframe for each species and year/stratum
# combinations, by pivoting the taxa column into columns
fish.pa.pivot.FD <-
  lapply(fish.abun.clust.gr.split,
         function(df) {
           df <- df %>%
             mutate(
               presence_absence = as.numeric(group_biomass > 0),
               stratum_year = str_c(stratum, "_", year_surv)
             ) %>%
             pivot_wider(names_from = clusterID,
                         values_from = presence_absence,
                         id_cols = stratum_year)
         })

# Still with lapply, we turn our dataframes into matrices and make sure their
# rownames correspond to the stratum/year combination
fish.pa.pivot.FD.mat <- lapply(fish.pa.pivot.FD, function(df) {
  # remove first column to be casted correctly as integer
  mat <- as.matrix(df[,-1])
  mat[is.na(mat)] <- 0
  rownames(mat) <- df$stratum_year
  return(mat)
})

# Still with lapply, we compute the beta diversity (see ?betadiver, equation 15 which is podoni:Jaccard), cast as a
# matrix and return
fish.beta.FD <- lapply(fish.pa.pivot.FD.mat, function(mat) {
  mat <- as.matrix(betadiver(mat, 15))
  return(mat)
})

# Still with lapply, we turn our matrix back into a usable, joinable dataframe
# - we first collect the rownames, pivot, rename, separate the columns and
# wrangle the data so that the "to" year is later than the "from" year
fish.beta.FD.df <- lapply(fish.beta.FD, function(mat) {
  as.data.frame(mat) %>%
    tibble::rownames_to_column() %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    unique() %>%
    dplyr::rename(from = rowname, to = name) %>%
    tidyr::separate(from, c("from_stratum", "from_year"), "_") %>%
    tidyr::separate(to, c("to_stratum", "to_year"), "_") %>%
    select(-to_stratum) %>%
    dplyr::rename(stratum = from_stratum) %>%
    filter(to_year > from_year)
})

# We turn the list into a fully combined dataframe
fish.beta.FD.df.full <- bind_rows(fish.beta.FD.df)

# We join to the stratum spatial data
fish.beta.FD.df.full.joined <- fish.beta.FD.df.full %>%
  left_join(stratum.shpfile, by = "stratum") %>%
  st_as_sf()

###MAP - BETA DISP.GR----
# We select the years and plot
beta.map.Disp.GR <- fish.beta.FD.df.full.joined %>%
  filter(from_year == 1996, to_year == 2017) %>%
  ggplot() +
  geom_sf(aes(fill = value),
          color = NA) +
  #coord_sf (xlim = c(-61, -46), ylim = c(42.5, 56)) +
  labs(
    title = "Dispersal Group Turn Over",
    subtitle = "1996-2017",
    fill = "Dissimilarity",
    x = "Longitude",
    y = "Latitude"
  ) +
  scale_fill_viridis_c () +
  theme_light()+
  theme(plot.title = element_text(lineheight = .8, size = 15), # title
        #axis.text.x = element_blank(), # remove x axis labels
        #axis.text.y = element_blank(), # remove y axis labels  
        axis.ticks = element_blank(), # remove axis ticks
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
beta.map.Disp.GR
## BETA SPECIES ----
# We split our main dataframe into a list of each dataframe corresponding to
# each stratum

fish.abun.clean.split <- fish.abun.clean %>%
  # We exclude stratum for which there is only data for 1 year
  filter(!(
    stratum %in%
      c(
        501,
        502,
        503,
        504,
        505,
        506,
        507,
        508,
        509,
        510,
        511,
        512,
        513,
        514,
        515,
        516,
        517,
        518,
        519,
        710,
        915,
        917,
        918
      )
  )) %>%
  # Split things
  split.data.frame(.$stratum)

# Using lapply, we apply a function to each data frame in the list. Here we
# create a presence absence dataframe for each species and year/stratum
# combinations, by pivoting the taxa column into columns
fish.pa.pivot <-
  lapply(fish.abun.clean.split,
         function(df) {
           df <- df %>%
             mutate(
               presence_absence = as.numeric(group_biomass > 0),
               stratum_year = str_c(stratum, "_", year_surv)
             ) %>%
             pivot_wider(names_from = taxa_name,
                         values_from = presence_absence,
                         id_cols = stratum_year)
         })

# Still with lapply, we turn our dataframes into matrices and make sure their
# rownames correspond to the stratum/year combination
fish.pa.pivot.mat <- lapply(fish.pa.pivot, function(df) {
  # remove first column to be casted correctly as integer
  mat <- as.matrix(df[,-1])
  mat[is.na(mat)] <- 0
  rownames(mat) <- df$stratum_year
  return(mat)
})

# Still with lapply, we compute the beta diversity (see ?betadiver, equation 15 which is podoni:Jaccard), cast as a
# matrix and return
fish.beta <- lapply(fish.pa.pivot.mat, function(mat) {
  mat <- as.matrix(betadiver(mat, 15))
  return(mat)
})

# Still with lapply, we turn our matrix back into a usable, joinable dataframe
# - we first collect the rownames, pivot, rename, separate the columns and
# wrangle the data so that the "to" year is later than the "from" year
fish.beta.df <- lapply(fish.beta, function(mat) {
  as.data.frame(mat) %>%
    tibble::rownames_to_column() %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    unique() %>%
    dplyr::rename(from = rowname, to = name) %>%
    tidyr::separate(from, c("from_stratum", "from_year"), "_") %>%
    tidyr::separate(to, c("to_stratum", "to_year"), "_") %>%
    select(-to_stratum) %>%
    dplyr::rename(stratum = from_stratum) %>%
    filter(to_year > from_year)
})

# We turn the list into a fully combined dataframe
fish.beta.df.full <- bind_rows(fish.beta.df)

# We join to the stratum spatial data
fish.beta.df.full.joined <- fish.beta.df.full %>%
  left_join(stratum.shpfile, by = "stratum") %>%
  st_as_sf()

###MAP BETA.SPECIES----
# We select the years and plot
beta.map.species <- fish.beta.df.full.joined %>%
  filter(from_year == 1996, to_year == 2017) %>%
  ggplot() +
  geom_sf(aes(fill = value), 
          color = NA) +
  #coord_sf (xlim = c(-61, -46), ylim = c(42.5, 58)) +
  labs(
    title = "Species Turn Over",
    subtitle = "1996-2017",
    fill = "Dissimilarity",
    x = "Longitude",
    y = "Latitude"
  ) +
  scale_fill_viridis_c () +
  theme_light() +
  theme(legend.position = "NULL")+
  theme(plot.title = element_text(lineheight = .8, size = 15), # title
        #axis.text.x = element_blank(), # remove x axis labels
        #axis.text.y = element_blank(), # remove y axis labels  
        axis.ticks = element_blank(), # remove axis ticks
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#ggforce::facet_wrap_paginate(~year_surv,
# nrow = 2, ncol = 2, page = 2)

beta.map <- beta.map.species + beta.map.Disp.GR

# ggsave("BetaDiversity.png",
#        beta.map,
#        width =  10,
#        height = 10)

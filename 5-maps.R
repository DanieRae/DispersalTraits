library(sf)
library(ggplot2)
library(ggforce)
library(patchwork)

# library(tmap)

# Import the strata data ------
#Strata depth
strata.depth <- read.csv("shrimp_current_length_data.csv", fileEncoding ="UTF-8-BOM", stringsAsFactors = TRUE)

strata.depth <- select (strata.depth, c("stratum", "depth")) %>%
  group_by(stratum)%>%
  dplyr::summarize(depth = mean(depth))

# read all files of the shapefile
stratum.shpfile <- st_read("strata/all_strata.shp")

missingdepth<- stratum.new[!stratum.new$stratum %in% strata.depth$stratum,]

#Strata geometry
stratum.new <- select(stratum.shpfile, -c ("DIV", "strtm_t")) %>% 
  # Validating geometries to make operations on them 
  st_make_valid() %>%  
  # grouping by the stratum id
  group_by(stratum) %>% 
  # Combining all polygons that have the same stratum ID
  dplyr::summarize(geometry = st_union(geometry))

# Trying to fix broken/missing stratums
stratum.new.broken <- select(stratum.shpfile, -c ("DIV", "strtm_t")) %>% 
  # Validating geometries to make operations on them 
  st_make_valid() %>%  
  # grouping by the stratum id
  group_by(stratum) %>% 
  # Combining all polygons that have the same stratum ID
  dplyr::summarize(geometry = st_union(geometry)) %>% 
  # trying to get brokenstratums
  dplyr::filter(nchar(stratum) > 3)


plot(stratum.new$geometry)

stratum.new <- merge(x =stratum.new, y=strata.depth,by = "stratum" )

##Lets add the cluster IDs to the abundance data. Seems like it could be useful 

clust.ID <- select (fish.traits.40NA.clust, c ("clusterID"))
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

stratum.FE <- merge(x = fish.abun.gr, y= stratum.new, by = "stratum", all.y = TRUE)%>%
  sf::st_as_sf()

stratum.FE.depth <- merge (x = stratum.FE, y = strata.depth, by = "stratum")%>%
  sf::st_as_sf()

#MAP TESTS ----
map.test <- stratum.FE.depth %>%
  ggplot()+
  geom_sf(aes(fill = depth))+
  scale_fill_gradient(low = 'white', high = "purple")

map1 <- stratum.FE %>% facet_wrap(~year_surv)
  filter(year_surv %in% c(1995)) %>% 
  ggplot() +
  geom_sf(aes(fill = Unique_FE)) +
  ggtitle(label = "Dispersal Richness") +
  scale_fill_viridis_c () +
  theme_light() +
  
  #ggforce::facet_wrap_paginate(~year_surv,
                              # nrow = 2, ncol = 2, page = 2)

ggsave("plot1995.png", map1, width =  10, height = 10)


#Shannon Diversity on Functional Groups ####
#Effective functional diversity with Shannon in clustered groups
fish.abun.clean.summary.FD<-fish.abun.clust %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, vessel, trip, set, clusterID)%>%
  dplyr::summarize(group_biomass = sum(density_kgperkm2))%>%
  #calculate average biomass per stratum for each functional group
  group_by(stratum, year_surv, clusterID)%>%
  dplyr::summarize(group_biomass = mean(group_biomass))%>% 
    ungroup()
  
    
    #calculate the Hill number for each stratum in each year
effective.FD <- fish.abun.clean.summary.FD %>%
  group_by(stratum, year_surv)%>%
  dplyr::summarize(effective_species = exp(diversity(group_biomass, "shannon")))

  
# BETA Clustered------------------
  # We split our main dataframe into a list of each dataframe corresponding to
  # each stratum
  
  fish.abun.clean.summary.FD.split <- fish.abun.clean.summary.FD %>% 
    # We exclude stratum for which there is only data for 1 year
    filter(!(stratum %in% 
               c(501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511,
                 512, 513, 514, 515, 516, 517, 518, 519, 710, 915, 917, 918))) %>%
    # Split things
    split.data.frame(.$stratum)
  
  # Using lapply, we apply a function to each data frame in the list. Here we 
  # create a presence absence dataframe for each species and year/stratum 
  # combinations, by pivoting the taxa column into columns
  fish.pa.pivot.FD <- 
    lapply(fish.abun.clean.summary.FD.split, 
           function(df){
             df <- df %>% 
               mutate(presence_absence = as.numeric(group_biomass > 0), 
                      stratum_year = str_c(stratum, "_", year_surv)) %>% 
               pivot_wider(names_from = clusterID, values_from = presence_absence, 
                           id_cols = stratum_year)
           }
    )
  
  # Still with lapply, we turn our dataframes into matrices and make sure their 
  # rownames correspond to the stratum/year combination
  fish.pa.pivot.FD.mat <- lapply(fish.pa.pivot.FD, function(df) {
    # remove first column to be casted correctly as integer
    mat <- as.matrix(df[,-1])
    mat[is.na(mat)] <- 0
    rownames(mat) <- df$stratum_year
    return(mat)
  }
  )
  
  # Still with lapply, we compute the beta diversity (see ?betadiver, equation 15 which is podoni:Jaccard), cast as a 
  # matrix and return
  fish.beta.FD <- lapply(fish.pa.pivot.FD.mat, function(mat) {
    mat <- as.matrix(betadiver(mat, 15))
    return(mat)
  }
  )
  
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
    left_join(stratum.new, by = "stratum") %>% 
    st_as_sf()

#Beta Clustered Map 1----
# We select the years and plot
beta.map.FD <- fish.beta.FD.df.full.joined %>% 
    filter(from_year == 1996, to_year == 2017) %>% 
    ggplot() +
    geom_sf(aes(fill = value)) +
    coord_sf (xlim = c(-61, -46), ylim = c(42.5, 58))+
    labs(title = "Functional Group turn over", subtitle = "1996-2017", fill = "Dissimilarity", x = "Longitude", y = "Latitude") +
    scale_fill_viridis_c () +
    theme_light()
    
  
#Shannon Diversity on Species ####
#Effective species Diversity with Shannnon on species     
# We first summarize the biomass data to get the man biomass for each year 
# and stratum
fish.abun.clean.summary <- fish.abun.clean %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, vessel, trip, set, taxa_name)%>%
  dplyr::summarize(group_biomass = sum(density_kgperkm2))%>%
  #calculate average biomass per stratum for each functional group
  group_by(stratum, year_surv, taxa_name) %>%
  dplyr::summarize(group_biomass = mean(group_biomass)) %>% 
  ungroup()

# BETA Species------------------
# We split our main dataframe into a list of each dataframe corresponding to
# each stratum
  
fish.abun.clean.summary.split <- fish.abun.clean.summary %>% 
  # We exclude stratum for which there is only data for 1 year
  filter(!(stratum %in% 
             c(501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511,
               512, 513, 514, 515, 516, 517, 518, 519, 710, 915, 917, 918))) %>%
  # Split things
  split.data.frame(.$stratum)

# Using lapply, we apply a function to each data frame in the list. Here we 
# create a presence absence dataframe for each species and year/stratum 
# combinations, by pivoting the taxa column into columns
fish.pa.pivot <- 
  lapply(fish.abun.clean.summary.split, 
         function(df){
           df <- df %>% 
             mutate(presence_absence = as.numeric(group_biomass > 0), 
                    stratum_year = str_c(stratum, "_", year_surv)) %>% 
             pivot_wider(names_from = taxa_name, values_from = presence_absence, 
                         id_cols = stratum_year)
         }
  )

# Still with lapply, we turn our dataframes into matrices and make sure their 
# rownames correspond to the stratum/year combination
fish.pa.pivot.mat <- lapply(fish.pa.pivot, function(df) {
  # remove first column to be casted correctly as integer
  mat <- as.matrix(df[,-1])
  mat[is.na(mat)] <- 0
  rownames(mat) <- df$stratum_year
  return(mat)
}
)

# Still with lapply, we compute the beta diversity (see ?betadiver, equation 15 which is podoni:Jaccard), cast as a 
# matrix and return
fish.beta <- lapply(fish.pa.pivot.mat, function(mat) {
  mat <- as.matrix(betadiver(mat, 15))
  return(mat)
}
)

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
  left_join(stratum.new, by = "stratum") %>% 
  st_as_sf()

#Beta map 1----
# We select the years and plot
beta.map <- fish.beta.df.full.joined %>% 
  filter(from_year == 1996, to_year == 2017) %>% 
  ggplot() +
  geom_sf(aes(fill = value)) +
  coord_sf (xlim = c(-61, -46), ylim = c(42.5, 58))+
  labs(title= "Species turn over",subtitle = "1996-2017", fill = "Dissimilarity", x = "Longitude", y = "Latitude") +
  scale_fill_viridis_c () +
  theme_light()+
  theme(legend.position= "NULL")
#ggforce::facet_wrap_paginate(~year_surv,
# nrow = 2, ncol = 2, page = 2)

beta.map1 <- beta.map+beta.map.FD

ggsave("BetaDiversity1996-2017.png", beta.map1, width =  10, height = 10)

# ALPHA Maps-----

effective.species <- fish.abun.clean.summary  %>%
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
  labs(title = "Effective Species Diversity") +
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
  ggtitle(label = "Effective Functional Group Diversity") +
  scale_fill_viridis_c () +
  theme_light() +
  facet_wrap(~year_surv)
#ggforce::facet_wrap_paginate(~year_surv,
# nrow = 2, ncol = 2, page = 2)

ggsave("EffectiveFD1995-2010.png", map3, width =  10, height = 10)

# MAP
map4 <- FEve.comm.depth %>% 
  filter(year %in% c(1995,2000,2005,2010)) %>% 
  ggplot() +
  geom_sf(aes(fill = V1)) +
  ggtitle(label = "Functional Diversity") +
  scale_fill_viridis_c () +
  theme_light() +
  facet_wrap(~year)

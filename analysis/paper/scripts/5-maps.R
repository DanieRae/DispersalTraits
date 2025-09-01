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
# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")

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
library(rnaturalearth)
library(rnaturalearthdata)


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
          color = "black",
          size = 0.1)+
  coord_sf (xlim = c(-59, -47))+
  scale_fill_cmocean(name = "deep")+
  theme_light() +
  labs(fill = "Depth (m)", colour = "") + # legend titles
  theme(plot.title = element_text(lineheight = .8, size = 15, hjust = 0.2), # title
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        #axis.ticks = element_blank(), # remove axis ticks
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), # remove grid lines
        legend.position = c(0.05, 0.0), # Adjust coordinates to your preference
        legend.justification = c("left","bottom"),# Optional: Adjust justification of the legend box
        legend.direction = "horizontal",
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(face = "bold"))

map.depth
# ggsave(here("analysis", "figures", "StratumDepth.png"),
#             map.depth, width =  10, height = 10)


# NFL Map
#Load NAFO divisions
nafo <- st_read(here ("analysis", "data","raw_data","spatial","nafo division"))


# Get Canada and ocean basemap
canada <- ne_countries(country = "Canada", scale = "medium", returnclass = "sf")
ocean <- ne_download(scale = "medium", type = "ocean", category = "physical", returnclass = "sf")

# Define times zone
crs_proj <- 32622  # UTM Zone 22N


# Transform both layers
nafo_proj <- st_transform(nafo, crs_proj)
canada_proj <- st_transform(canada, crs_proj)


# Create a label point column using projected CRS
nafo_proj <- nafo_proj %>%
  mutate(label_pt = st_centroid(geometry))

#Manually fixing the position of some labels
nafo_proj$label_pt[nafo_proj$Division == "4T"] <- st_sfc(st_point(c(-364927.1, 5305732)), crs = crs_proj)
nafo_proj$label_pt[nafo_proj$Division == "4W"] <- st_sfc(st_point(c(-276927.1, 4950524)), crs = crs_proj)
nafo_proj$label_pt[nafo_proj$line_id == "4R-1"] <- st_sfc(st_point(c(0, 0)), crs = crs_proj)
nafo_proj$label_pt[nafo_proj$Division == "3O"] <- st_sfc(st_point(c(366927.1, 4950524)), crs = crs_proj)
nafo_proj$label_pt[nafo_proj$Division == "3N"] <- st_sfc(st_point(c(656927.1, 4950524)), crs = crs_proj)
nafo_proj$label_pt[nafo_proj$Division == "3M"] <- st_sfc(st_point(c(906927.1, 5150524)), crs = crs_proj)
nafo_proj$label_pt[nafo_proj$line_id == "4Vs"] <- st_sfc(st_point(c(-6927.1, 4950524)), crs = crs_proj)
nafo_proj$label_pt[nafo_proj$line_id == "4S"] <- st_sfc(st_point(c(-206927.1, 5470524)), crs = crs_proj)

# automatically creaking the dimension of the map
bbox_proj <- st_bbox(c(xmin = -603435.4, xmax = 1004576.3 , ymin = 4883720.4, ymax = 6650203.4), crs = st_crs(crs_proj))

#plot the map
NFL.MAP <- ggplot() +
  # Optional: ocean background
  geom_sf(data = st_as_sfc(bbox_proj), fill = "lightblue", color = NA) +
  # Canada landmass
  geom_sf(data = canada_proj, fill = "darkseagreen", color = "black", size = 0.3) +
  #label landmass
      annotate("text", x = 120085.99, y =  5398478, label = "Newfoundland", color = "Black", size =3) +
  # NAFO divisions
  geom_sf(data = nafo_proj, fill = NA, color = "gray10", size = 0.4) +
  # Labels (now safely centered)
  geom_sf_text(
    data = nafo_proj %>%
      filter(!Division %in% c("4R-1", "0-B", "4X")),
    aes(geometry = label_pt, label = Division),
    stat = "sf_coordinates",
    size = 3, color = "black"
  ) +
  coord_sf(xlim = c(bbox_proj["xmin"], bbox_proj["xmax"]),
           ylim = c(bbox_proj["ymin"], bbox_proj["ymax"]),
           expand = FALSE) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "gray90"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),

      )

NFL.MAP
# ggsave(here("analysis", "figures", "NFL.MAP.png"),
#        NFL.MAP, width =  10, height = 10)

# Figure 1
#need to match the latitude limits

figure1<-NFL.MAP+map.depth
# ggsave(here("analysis", "figures", "figure1.png"),
#        figure1, width =  10, height = 10)


#END-----


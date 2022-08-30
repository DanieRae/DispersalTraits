### Make SSPM model
library(sspm)

# The package needs to think we are using points, so we make some
stratum_points <- stratum.shpfile.abun %>% 
  st_point_on_surface() %>% 
  select(stratum, new_bin)

plot <-stratum.shpfile.abun %>% 
  ggplot() + 
  geom_sf() + 
  geom_sf(data = stratum_points)

#Need to take this sf file and pull just the data without the geometries, we are using this rather than the "abun.feve.bin.noNA" data table because this table has more strata than we have shpfiles for.

stratum.abun <- as.data.frame(stratum.shpfile.abun)

stratum.abun<- select(stratum.abun, -c("geometry"))

stratum.shpfile.geom <- select(stratum.shpfile.abun, c("stratum", "geometry"))%>%
  sf::st_as_sf()

stratum.shpfile.geom <- distinct(stratum.shpfile.geom) %>%
  sf::st_as_sf()

# The package needs to know the boundaries of all the polygons, but because of all
# the gaps we have to take convex hull and buffer
stratum.shpfile.abun.combined <- st_convex_hull(stratum.shpfile.abun) %>% 
  st_union() %>% st_make_valid() %>% st_buffer(.1) %>% st_make_valid() %>% 
  st_as_sf() %>% 
  mutate(stratum = 1)  %>% 
  dplyr::rename(geometry = x)

# Take the data, join them to the new points, 
# then create two columns, for Unique ID and for time
stratum.abun.spatial <- stratum.abun %>% 
  left_join(stratum_points) %>% 
  st_as_sf() %>% 
  mutate(uniqueID = paste0(stratum, new_bin), 
         year_basis = as.factor(as.numeric(str_sub(as.character(new_bin), 
                                                   end = 4))))

# Use sspm
# we create the dataset object
dataset <- spm_as_dataset(data = stratum.abun.spatial, 
                          name = "abundance_model",
                          time = "year_basis", 
                          uniqueID = "uniqueID")

# We create the boundary object
bound <- spm_as_boundary(stratum.shpfile.abun.combined, 
                         boundary = "stratum", 
                         patches = select(stratum.shpfile, 
                                          "stratum"), points = NULL)
 

# dataset_smoothed <- dataset %>%  
#   spm_smooth(bin_mean_biomass ~ stratum + new_bin + smooth_space(bs = "mrf"),
#              boundaries = bound, 
#              family = tw)

# helper function 
stratum_mrf <- sspm:::ICAR(data_frame = dataset@data, boundaries = dataset@boundaries, dimension = "space", time = dataset@time, k=10, bs = "mrf", xt = NULL )

stratum_mrf_pen <-sspm:::ICAR_space(patches = stratum.shpfile.geom, space = "stratum")


#Exploring gams with eric 
gam_test <- 
  gam(
    log(bin_sd_biomass) ~ bin_mean_FEve + new_bin + s(stratum, bs = "mrf", xt = list(penalty = stratum_mrf_pen)),
    data = stratum.abun,
    method = "REML",
    family = gaussian()
  )


dataset_smoothed <- dataset %>%  
  spm_smooth(bin_sd_biomass ~ bin_mean_FEve + smooth_space(bs = "mrf"),
             boundaries = bound, 
             family = Gamma(link="log"))

fit <- spm_smoothed_fit(dataset_smoothed)[[1]]
View(fit$smooth)

View(fit$smooth[[1]][["xt"]][["penalty"]])


dataset_smoothed_plot <-
  ggplot(data = abun.feve.bin.noNA.spatial, aes(y = log(bin_sd_biomass), x = bin_mean_FEve)) +
  geom_point() +
  geom_line(aes(y = fitted(fit)),
            colour = "red", size = 1.2) +
  theme_bw()+
  labs(title = "Marine Community Stability in Response to Dispersal Diversity") +
  xlab("Dispersal Diversity")+
  ylab ("Biomass")

veca <- unique(sort(stratum.shpfile.geom$stratum))
vecb <- unique(sort(abun.feve.bin.noNA$stratum))
vecc <- unique(sort(stratum.abun$stratum))

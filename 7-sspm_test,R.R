### Make SSPM model
library(sspm)

# The package needs to think we are using points, so we make some
stratum_points <- stratum.shpfile.abun %>% 
  st_point_on_surface() %>% 
  select(stratum, new_bin)

stratum.shpfile.abun %>% ggplot() + geom_sf() + geom_sf(data = stratum_points)

# The package needs to know the boundaries of all the polygons, but because of all
# the gaps we have to take convex hull and buffer
stratum.shpfile.abun.combined <- st_convex_hull(stratum.shpfile.abun) %>% 
  st_union() %>% st_make_valid() %>% st_buffer(.1) %>% st_make_valid() %>% 
  st_as_sf() %>% 
  mutate(stratum = 1) %>% 
  rename(geometry = x)

# Take the data, join them to the new points, 
# then create two columns, for Unique ID and for time
abun.feve.bin.noNA.spatial <- abun.feve.bin.noNA %>% 
  left_join(stratum_points) %>% 
  st_as_sf() %>% 
  mutate(uniqueID = paste0(stratum, new_bin), 
         year_basis = as.factor(as.numeric(str_sub(as.character(new_bin), 
                                                   end = 4))))

# Use sspm
# we create the dataset object
dataset <- spm_as_dataset(data = abun.feve.bin.noNA.spatial, 
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

dataset_smoothed <- dataset %>%  
  spm_smooth(bin_mean_biomass ~ smooth_space(bs = "mrf"),
             boundaries = bound, 
             family = tw)

fit <- spm_smoothed_fit(dataset_smoothed)[[1]]
View(fit$smooth)

View(fit$smooth[[1]][["xt"]][["penalty"]])
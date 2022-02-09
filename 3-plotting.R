library(sf)
library(ggplot2)
# library(tmap)

# Import the data ---------------------------------------------------------

# read all files of the shapefile
stratum <- st_read("data/all_strata.shp")

# simple plotting
plot(stratum$geometry)
#LIBRARIES ----
library(sspm)


stratum.shpfile.spdiv <-
  merge(stratum.shpfile.adjusted, effective.species.stab, by = "stratum")

# The package needs to think we are using points, so we make some
stratum_points_spdiv <- stratum.shpfile.spdiv %>%
  st_point_on_surface() %>%
  select(stratum, new_bin)


#Need to take this sf file and pull just the data without the geometries, we are using this rather than the "abun.feve.bin.noNA" data table because this table has more strata than we have shpfiles for.

stratum.spdiv <- as.data.frame(stratum.shpfile.spdiv)

stratum.spdiv <- select(stratum.spdiv, -c("geometry"))

stratum.shpfile.spdiv.geom <- select(stratum.shpfile.spdiv, c("stratum", "geometry")) %>%
  sf::st_as_sf()

stratum.shpfile.spdiv.geom <- distinct(stratum.shpfile.spdiv.geom) %>%
  sf::st_as_sf()


# The package needs to know the boundaries of all the polygons, but because of all
# the gaps we have to take convex hull and buffer
stratum.shpfile.spdiv.combined <- st_convex_hull(stratum.shpfile.spdiv) %>%
  st_union() %>% st_make_valid() %>% st_buffer(.1) %>% st_make_valid() %>%
  st_as_sf() %>%
  mutate(stratum = 1)  %>%
  dplyr::rename(geometry = x)

# Take the data, join them to the new points,
# then create two columns, for Unique ID and for time
stratum.spdiv.spatial <- stratum.spdiv %>%
  left_join(stratum_points_spdiv) %>%
  st_as_sf() %>%
  mutate(uniqueID = paste0(stratum, new_bin),
         year_basis = as.factor(as.numeric(str_sub(as.character(new_bin),
                                                   end = 4))))

# Use sspm
# we create the dataset object
dataset_spdiv <- spm_as_dataset(data = stratum.spdiv.spatial,
                                name = "abundance_model",
                                time = "year_basis",
                                uniqueID = "uniqueID")

# We create the boundary object
bound <- spm_as_boundary(stratum.shpfile.spdiv.combined,
                         boundary = "stratum",
                         patches = select(stratum.shpfile,
                                          "stratum"), points = NULL)


# dataset_smoothed <- dataset %>%
#   spm_smooth(bin_mean_biomass ~ stratum + new_bin + smooth_space(bs = "mrf"),
#              boundaries = bound,
#              family = tw)


#Val written model with spm, using predctor variable that is not easily accepted by spm
dataset_smoothed <- dataset_spdiv %>%
  spm_smooth(bin_sd_biomass ~ Z_effectivesp + smooth_space(bs = "mrf"),
             boundaries = bound,
             family = Gamma(link = "log"),
             predict = FALSE)

fit <- spm_smoothed_fit(dataset_smoothed)[[1]]

#Exploring gams with eric ----

# helper function
stratum_mrf <- sspm:::ICAR(data_frame = dataset_spdiv@data, boundaries = bound, dimension = "space", time = dataset_spdiv@time, k = 10, bs = "mrf", xt = NULL )

#penalty
stratum_mrf_pen <- sspm:::ICAR_space(patches = stratum.shpfile.spdiv.geom, space = "stratum")

stratum.spdiv.mut <- stratum.spdiv %>%
  mutate(stratum = as.factor(stratum))

#MODEL ----
gam.stability.spdiv <-
  gam(
    log((bin_sd_biomass)^-1) ~ Z_effectivesp + s(new_bin, bs = "re", k = 6) + s(stratum, bs = "mrf", xt = list(penalty = stratum_mrf_pen), k= 200) ,
    data = stratum.spdiv.mut,
    method = "REML",
    type = "terms",
    family = gaussian()
  )

gam.check(gam.stability.spdiv, rep = 500)
plot(gam.stability.spdiv)


##ERICS CODE ----
summary(gam.stability.spdiv)
term_predictors.spdiv <- predict(gam.stability.spdiv, type = "terms", se.fit = TRUE)
disperse_div_fit_spdiv <- as.vector(term_predictors.spdiv$fit[,"Z_effectivesp"]) #the as.vector part is just to make sure this is a vector, and not a 1D matrix, for plotting.

disperse_intercept_spdiv <- as.vector(attr(term_predictors.spdiv,"const")) #like this because the intercept is a single number, so it's just stored with the predicted values as a constant
disperse_div_fit_spdiv <- disperse_div_fit_spdiv + disperse_intercept_spdiv

#Confidence Int
disp_div_se_spdiv <- as.vector(term_predictors.spdiv$se.fit[,"Z_effectivesp"])

##PLOT GAM ----
dataset_smoothed_plot_spdiv <-
  ggplot(data = stratum.spdiv.mut, aes(y = log((bin_sd_biomass)^-1), x = Z_effectivesp)) +
  geom_point() +
  geom_line(aes(y = disperse_div_fit_spdiv),
            colour = "blue", size = 1.2) +
  geom_ribbon(aes(ymin = disperse_div_fit_spdiv - 1.96*disp_div_se_spdiv,
                  ymax = disperse_div_fit_spdiv + 1.96*disp_div_se_spdiv),
              fill = "blue", alpha = 0.2) +
  theme_bw() +
  xlab("Taxonomic Diversity (z-score scaled)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

dataset_smoothed_plot_spdiv

stabilitymaps <- dataset_smoothed_plot + dataset_smoothed_plot_spdiv

ggsave(here("analysis", "figures", "stabilitymaps.png"),
       stabilitymaps,
       height = 10,
       width = 20)

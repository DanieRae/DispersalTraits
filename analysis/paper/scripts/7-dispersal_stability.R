#LIBRARIES ----
# install.packages("devtools")
# devtools::install_github("pedersen-fisheries-lab/sspm")
# if you wish to build the vignettes:
# devtools::install_github("pedersen-fisheries-lab/sspm", build_vignettes = TRUE)

library(sspm)
library(patchwork)


# Make SSPM model ----

stratum.shpfile.abun <-
  merge(stratum.shpfile.adjusted, abun.feve.bin.noNA, by = "stratum")

# The package needs to think we are using points, so we make some
stratum_points <- stratum.shpfile.abun %>%
  st_point_on_surface() %>%
  select(stratum, new_bin)


#Need to take this sf file and pull just the data without the geometries, we are using this rather than the "abun.feve.bin.noNA" data table because this table has more strata than we have shpfiles for.

stratum.abun <- as.data.frame(stratum.shpfile.abun)

stratum.abun <- select(stratum.abun, -c("geometry"))

stratum.shpfile.geom <- select(stratum.shpfile.abun, c("stratum", "geometry")) %>%
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


#Val written model with spm, using predctor variable that is not easily accepted by spm
dataset_smoothed <- dataset %>%
  spm_smooth(inverse_sd ~ Z_FEve + smooth_space(bs = "mrf"),
             boundaries = bound,
             family = Gamma(link = "log"),
             predict = FALSE)

fit <- spm_smoothed_fit(dataset_smoothed)[[1]]

#Exploring gams with eric ----

# helper function
stratum_mrf <- sspm:::ICAR(data_frame = dataset@data, boundaries = bound,
                           dimension = "space", time = dataset@time, k = 10, bs = "mrf", xt = NULL )

#penalty
stratum_mrf_pen <- sspm:::ICAR_space(patches = stratum.shpfile.geom, space = "stratum")

stratum.abun.mut <- stratum.abun %>%
  mutate(stratum = as.factor(stratum))

#MODEL ----
gam.stability <-
  gam(
    log((bin_sd_biomass)^-1) ~ Z_FEve +Z_effectivesp + s(new_bin, bs = "re", k = 6) + s(stratum, bs = "mrf", xt = list(penalty = stratum_mrf_pen), k= 200),
    data = stratum.abun.mut,
    method = "REML",
    type = "terms",
    family = gaussian()
  )

gam.check(gam.stability, rep = 500)
plot(gam.stability)


##ERICS CODE ----
summary(gam.stability)
term_predictors <- predict(gam.stability, type = "terms", se.fit = TRUE)

disperse_div_predictions <- stratum.abun.mut %>%
  mutate(
    #fitted values for each term for confidence intervals
    FEve_fit = as.vector(term_predictors$fit[,"Z_FEve"]),
    effectivesp_fit = as.vector(term_predictors$fit[,"Z_effectivesp"]),
    # standard errors for each term for confidence intervals
    FEve_se = as.vector(term_predictors$se.fit[,"Z_FEve"]),
    effectivesp_se = as.vector(term_predictors$se.fit[,"Z_effectivesp"]),
    #residual values
    resid =  residuals(gam.stability))


##PLOT GAM ----
FEve_partial_resid_plot <-
  ggplot(data = disperse_div_predictions, aes(y = FEve_fit + resid, x = Z_FEve)) +
  geom_point() +
  geom_line(aes(y = FEve_fit),
            colour = "red", size = 1.2) +
  geom_ribbon(aes(ymin = FEve_fit - 1.96*FEve_se,
                  ymax = FEve_fit + 1.96*FEve_se),
              fill = "red", alpha = 0.2) +
  theme_bw() +
  xlab("Dispersal Diversity (z-score scaled)") +
  ylab("Effect of dispersal diversity") +
  theme(axis.title.x = element_text(size = 25,face = "bold"),
        axis.title.y = element_text(size = 25, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

FEve_partial_resid_plot

effectivesp_partial_resid_plot <-
  ggplot(data = disperse_div_predictions, aes(y = effectivesp_fit + resid, x = Z_effectivesp)) +
  geom_point() +
  geom_line(aes(y = effectivesp_fit),
            colour = "blue", size = 1.2) +
  geom_ribbon(aes(ymin = effectivesp_fit - 1.96*effectivesp_se,
                  ymax = effectivesp_fit + 1.96*effectivesp_se),
              fill = "blue", alpha = 0.2) +
  theme_bw() +
  xlab("Taxonomic diversity (z-score scaled)") +
  ylab("Effect of taxonomic diversity\non stability") +
  theme(axis.title.x = element_text(size = 25,face = "bold"),
        axis.title.y = element_text(size = 25, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

effectivesp_partial_resid_plot

diversity_effects_plot <- FEve_partial_resid_plot + effectivesp_partial_resid_plot

ggsave(here("analysis", "figures", "diversity_effects_gam.png"),
       diversity_effects_plot ,
       height = 10,
       width = 30)


#stratum effect plot

  #the as.vector part is just to make sure this is a vector, and not a 1D matrix, for plotting.

stratum_effect_plot <- gam.stability %>%
  mutate(stratum_effect = as.vector(term_predictors$fit[,"s(stratum)"])) %>%
  group_by(stratum) %>%
  summarize(stratum_effect = stratum_effect[1]) %>%
  left_join(stratum.shpfile.geom) %>%
  st_as_sf() %>%
  ggplot(aes(fill = stratum_effect)) +
  geom_sf()

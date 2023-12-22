### Running the functional diversity indices
library(FD)

##ABUNDANCE DATA WRANGLING ----
#abundance by site (strata,year)
abun.fish.wide <-
  as.data.table(spread(fish.abun.complete, taxa_name, group_biomass))

abun.fish.biomass <-
  select(abun.fish.wide,-c("stratum", "year_surv"))

abun.fish.names <- sort(names(abun.fish.biomass), decreasing = FALSE)

abun.fish.biomass <- abun.fish.biomass[, ..abun.fish.names]

abun.fish.biomass <- as.matrix(abun.fish.biomass)

# FD ----

functional.diversity <- dbFD(fish.traits.40NA, abun.fish.biomass, corr = "lingoes")

##FD SAVED FILE ----
saveRDS(functional.diversity,
        file = here("analysis", "data", "derived_data", "FunctionalDiv.rds"))

functional.diversity <-
  readRDS(here("analysis", "data", "derived_data", "FunctionalDiv.rds"))

#Run script to check difference in Dispersal diversity and effective species using MIN BIOMASS versus 0 BIOMASS#


# PART 1 - MIN BIOMASS ----

#FILLING ABUNDANCE DATA WITH MISSING SPECIES ABUNDANCES ----
#table of minimun biomass values
min_biomass <- fish.abun.filtered %>%
  group_by(taxa_name) %>%
  dplyr::summarise(min_biomass = min(density_kgperkm2))


fish.abun.complete.NA <- fish.abun.clean %>%
  #filter(year_surv == 1996)%>%
  #select(taxa_name, stratum, year_surv)%>%
  tidyr::complete(nesting(stratum, year_surv), taxa_name)

#Fill NA with min biomass
fish.abun.complete.minbio <- fish.abun.complete.NA %>%
  left_join(min_biomass) %>%
  mutate(final_biomass = ifelse(is.na(group_biomass), min_biomass, group_biomass)) %>%
  mutate(weight = ifelse(is.na(group_biomass), 0, 1)) %>%
  select(-c("group_biomass", "min_biomass"))

fish.abun.complete.noweight <- select(fish.abun.complete.minbio, -c(weight))

##WIDE----
#abundance by site (strata,year)
abun.fish.min.wide <-
  as.data.table(spread(fish.abun.complete.noweight, taxa_name, final_biomass))

abun.fish.min.biomass <-
  select(abun.fish.min.wide,-c("stratum", "year_surv"))

abun.fish.min.names <- sort(names(abun.fish.min.biomass), decreasing = FALSE)

abun.fish.min.biomass <- abun.fish.min.biomass[, ..abun.fish.min.names]


abun.fish.min.biomass <- as.matrix(abun.fish.min.biomass)


#FD --------
#If you need to re-run the diversity calculation you can use this code, and run it
functional.diversity.minbio <- dbFD(fish.traits.40NA, abun.fish.min.biomass, corr = "lingoes")
saveRDS(functional.diversity.minbio, file = "FunctionalDivMinBIO.rds")

#load functional diversity metrics that used min-biomass filled data
functional.diversity.minbio <-
  readRDS(here("analysis", "data", "derived_data", "FunctionalDivMinBIO.rds"))

year.strat <- select(abun.fish.min.wide, c("stratum", "year_surv"))

##FEVE ----
FEve.comm.min.biomass <-
  cbind(
    functional.diversity.minbio$FEve,
    stratum = year.strat$stratum,
    year = year.strat$year_surv
  )

FEve.comm.min.biomass <- as.data.frame(FEve.comm.min.biomass)

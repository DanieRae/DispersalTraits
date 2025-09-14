#Statistics and response variables#

# Install and load packages ----
#
# install.packages("grid")
# install.packages("gridExtra")
# install.packages("lme4")
# install.packages("sjPlot")
# install.packages("FD")
# install.packages("sf")
# install.packages("ggplot2")
# install.packages("ggforce")
# install.packages("patchwork")
# install.packages("dplyr")
# install.packages("stringr")
# install.packages("terra")
# install.packages("proj4")
# install.packages("spdep")
# install.packages("mgcv")
# install.packages("viridis")
# install.packages("AICcmodavg")


library(grid)
library(gridExtra)
library(lme4)
library(sjPlot)
library(FD)
library(sf)
library(ggplot2)
library(ggforce)
library(patchwork)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(terra)
library(proj4)
library(spdep)
library(mgcv)
library(viridis)
library(AICcmodavg)
library(here)

# --- Load fonctional diversity data --------

# Run script 6.1 first if "FunctionalDiv.rds" is missing
functional.diversity <- readRDS(
  here("analysis", "data", "derived_data","FunctionalDiv.rds"))

## functional evenness (FEve) ----
abun.fish.wide <- fish.abun.complete %>%
  pivot_wider(names_from = taxa_name, values_from = group_biomass) %>%
  as.data.table()

year.strat <- abun.fish.wide %>% select(stratum, year_surv)

FEve.comm <- data.frame(
  FEve = functional.diversity$FEve,
  stratum = year.strat$stratum,
  year = year.strat$year_surv
)


# ---- Define year bins (5-year averages) ----
#5 year averages for biomass, FEve, and effective species
tags <-
  c("1995-1999",
    "2000-2004",
    "2005-2009",
    "2010-2014",
    "2015-2017")

# ---- Biomass stability ----
# Total biomass per stratum/year

fish.abun.strata <- fish.abun.complete %>%
  group_by(stratum, year_surv) %>%
  dplyr::summarize(group_biomass = sum(group_biomass))

# Assign bins
abun.bin <- fish.abun.strata %>%
  mutate(new_bin = case_when(
      year_surv <= 1999 ~ tags[1],
      year_surv > 1999 & year_surv < 2005 ~ tags[2],
      year_surv >= 2005 & year_surv < 2010 ~ tags[3],
      year_surv >= 2010 & year_surv < 2015 ~ tags[4],
      year_surv >= 2015 & year_surv <= 2017 ~ tags[5]
    )
  )

abun.bin.mean <- abun.bin %>%
  #mean biomass for the binned years
  group_by(stratum, new_bin) %>%
  dplyr::summarize(bin_mean_biomass = mean(group_biomass))

# Bin-level biomass mean & sd
abun.bin.sd <- abun.bin %>%
  group_by(stratum, new_bin) %>%
  dplyr::summarise(bin_sd_biomass = sd(group_biomass))

abun.bin.biomass <- merge(abun.bin.mean, abun.bin.sd)


# ---- Functional Evenness (FEve) binned ----

FEve.bin <- FEve.comm %>%
  mutate(
    new_bin = case_when(
      year <= 1999 ~ tags[1],
      year > 1999 & year < 2005 ~ tags[2],
      year >= 2005 & year < 2010 ~ tags[3],
      year >= 2010 & year < 2015 ~ tags[4],
      year >= 2015 & year <= 2017 ~ tags[5]
    )
  )

FEve.bin.mean <- FEve.bin %>%
  #mean biomass for the binned years
  group_by(stratum, new_bin) %>%
  dplyr::summarize(bin_mean_FEve = mean(FEve))


# ---- Effective species diversity ----
effective.species.filled <- fish.abun.complete  %>%
  #calculate the Hill number for each stratum in each year
  group_by(stratum, year_surv) %>%
  dplyr::summarize(effective_species = exp(diversity(group_biomass, "shannon")))

effective.species.filled.bin <- effective.species.filled %>%
  mutate(
    new_bin = case_when(
      year_surv <= 1999 ~ tags[1],
      year_surv > 1999 & year_surv < 2005 ~ tags[2],
      year_surv >= 2005 & year_surv < 2010 ~ tags[3],
      year_surv >= 2010 & year_surv < 2015 ~ tags[4],
      year_surv >= 2015 & year_surv <= 2017 ~ tags[5]
    )
  )
effective.species.filled.bin.mean <- effective.species.filled.bin %>%
  #mean biomass for the binned years
  group_by(stratum, new_bin) %>%
  dplyr::summarize(bin_mean_effectivesp = mean(effective_species))

# ---- Merge FEve & biomass stability ----
abun.feve.bin <- inner_join(FEve.bin.mean, abun.bin.sd, by = c("stratum", "new_bin")) %>%
  mutate(
    new_bin = factor(new_bin),
    stratum = factor(stratum)
  )

abun.feve.bin.noNA <- abun.feve.bin %>%
  drop_na()%>%
  mutate(stratum = as.character(stratum))

# Scale variables
abun.feve.bin.noNA <- abun.feve.bin.noNA %>%
  mutate(
    Z_FEve = scale(bin_mean_FEve),
    inverse_sd = bin_sd_biomass^-1
  )


# ---- Merge effective species & stability ----

effective.species.stab <- inner_join(
  effective.species.filled.bin.mean, abun.bin.sd,
  by = c("stratum", "new_bin")
) %>%
  drop_na() %>%
  mutate(
    new_bin = factor(new_bin),
    Z_effectivesp = scale(bin_mean_effectivesp),
    stratum = as.character(stratum)
    )

# ---- Final combined dataset ----
stab.feve.species <- inner_join(abun.feve.bin.noNA, effective.species.stab,
                                by = c("stratum", "new_bin"))

# ---- Correlation: Taxonomic vs Functional Evenness ----
cor.test(~ Z_FEve + Z_effectivesp, data = stab.feve.species)

# ---- Plot: correlation ----

plot.binned.effective.species.FEve <- stab.feve.species %>%
  ggplot(aes(x = Z_effectivesp, y = Z_FEve)) +
  xlab("Taxonomic Diversity (effective species number)") +
  ylab("Dispersal Diversity (functional evenness)") +
  theme_light() +
  geom_point() +
  stat_smooth(method = glm) +
  # geom_line(aes(y = fitted(lm)),
  #           colour = "red", size = 1.2) +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 20,vjust = 0, face = "bold",
                                    family = "Arial Narrow"),
        plot.title = element_text(lineheight = .8, size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) # title

# Uncomment to save figure
# ggsave(here("analysis","figures",
#             "Influence of Effective Species on Dispersal Diversity.png"),
#        plot.binned.effective.species.FEve,
#        width =  15,
#        height = 10)


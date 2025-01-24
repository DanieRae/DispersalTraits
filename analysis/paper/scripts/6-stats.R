#Statistics and response variables#

# Install and load packages ----

install.packages("grid")
install.packages("gridExtra")
install.packages("lme4")
install.packages("sjPlot")
install.packages("FD")
install.packages("sf")
install.packages("ggplot2")
install.packages("ggforce")
install.packages("patchwork")
install.packages("dplyr")
install.packages("stringr")
install.packages("rgdal")
install.packages("rgeos")
install.packages("proj4")
install.packages("spdep")
install.packages("mgcv")
install.packages("viridis")
install.packages("AICcmodavg")


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
library(rgdal)
library(rgeos)
library(proj4)
library(spdep)
library(mgcv)
library(viridis)
library(AICcmodavg)
library(here)

#FD LOAD --------

# Run page 6.1 first if you don't have the RDS saved
functional.diversity <- readRDS(here("analysis", "data", "derived_data",
                                     "FunctionalDiv.rds"))

##FEVE ----
abun.fish.wide <-
  as.data.table(spread(fish.abun.complete, taxa_name, group_biomass))

year.strat <- select(abun.fish.wide, c("stratum", "year_surv"))

FEve.comm <-
  cbind(
    functional.diversity$FEve,
    stratum = year.strat$stratum,
    year = year.strat$year_surv
  )

FEve.comm <- as.data.frame(FEve.comm)


#YEAR BINS ----
#5 year averages for biomass, FEve, and effective species
tags <-
  c("1995-1999",
    "2000-2004",
    "2005-2009",
    "2010-2014",
    "2015-2017")

##STABILITY ----
#mean and sd of biomass per bin

##TOTAL BIOMASS per strata/year
fish.abun.strata <- fish.abun.complete %>%
  group_by(stratum, year_surv) %>%
  dplyr::summarize(group_biomass = sum(group_biomass))

abun.bin <- fish.abun.strata %>%
  mutate(
    new_bin = case_when(
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


abun.bin.sd <- abun.bin %>%
  #sd of the binned years
  group_by(stratum, new_bin) %>%
  dplyr::summarise(bin_sd_biomass = sd(group_biomass))

abun.bin.biomass <- merge(abun.bin.mean, abun.bin.sd)


##DISPERSAL DIVERSITY ----
#mean

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
  dplyr::summarize(bin_mean_FEve = mean(V1))


## EFFECTIVE SPECIES DIVERSITY ----
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

#DISPERSAL AND STAB ----
abun.feve.bin <- merge(FEve.bin.mean,  abun.bin.sd)
abun.feve.bin$new_bin <- as.factor(abun.feve.bin$new_bin)
abun.feve.bin$stratum <- as.factor(abun.feve.bin$stratum)

#too many rows with NAs to model
abun.feve.bin.noNA <- abun.feve.bin[complete.cases(abun.feve.bin), ]


#have to scale the variables as they have different scales
abun.feve.bin.noNA$Z_FEve <- scale(abun.feve.bin.noNA$bin_mean_FEve)
abun.feve.bin.noNA$inverse_sd <- abun.feve.bin.noNA$bin_sd_biomass^-1


# EFFECTIVE SPECIES STAB ----

effective.species.stab <- merge(effective.species.filled.bin.mean, abun.bin.sd)

effective.species.stab <- effective.species.stab[complete.cases(effective.species.stab), ]

effective.species.stab$new_bin <- as.factor(effective.species.stab$new_bin)

effective.species.stab$Z_effectivesp <- scale(effective.species.stab$bin_mean_effectivesp)

#TAXONOMIC VS DISPERSAL DIV MERGE----
stab.feve.species <- merge(abun.feve.bin.noNA, effective.species.stab)


### CORR - SPDIV/FEVE BINNED change this to corr ----
binned.effective.species.FEve <-
  merge(x = FEve.bin.mean, y = effective.species.filled.bin.mean)

lm <- cor.test( ~ Z_FEve + Z_effectivesp,
                data = stab.feve.species)
lm

##PLOT - SPDIV/FEVE BINNED ----
#Effective species/FEve pearson correlation plot

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



##PLOT - STAB/FEVE LOCAL ***ask eric----
#functional evenness and stab by community over time
abun.feve.bin %>%
  ggplot(aes(x = bin_mean_FEve, y = log(bin_sd_biomass))) +
  geom_point(size = 2) +
  #geom_text()+
  theme_bw() +
  xlab("DD") +
  ylab("Stab") +
  ggtitle("Change in Local Community Dispersal Diversity over a 20 Year Period") +
  geom_path() +
  scale_color_grey() +
  facet_wrap_paginate( ~ stratum,
                       nrow = 5,
                       ncol = 5,
                       page = 1) +
  theme(plot.title = element_text(lineheight = .8, size = 20, hjust = 0.5), # title
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #axis.text.x = element_blank(), # remove x axis labels
        #axis.text.y = element_blank(), # remove y axis labels
        axis.ticks = element_blank(), # remove axis ticks
        panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), # remove grid lines
        strip.text.x = element_text(size = 14,vjust = 0, face = "bold",
                                    family = "Arial Narrow"),
        legend.position = "none")

####NEED TO CLEAN PAST HERE####

#RATE OF BIOMASS CHANGE ----

#what is the rate of change of biomass in a strata from year1-yearn
Abun.fish.rate <- abun.feve.bin %>%
  group_by(stratum) %>%
  #because we don't have individual taxa we can;t look at sd
  arrange(stratum, new_bin) %>%
  dplyr::mutate(rate.change = log(bin_mean_FEve / lag(bin_mean_FEve)))

#this could be used for the rate of compositional change but no with the current state of summed biomass accross the strata
#Abun.fish.rate.SD <-Abun.fish.rate %>%
#  group_by(stratum,year_surv)%>%
#  dplyr:: summarise(sd2.rate.change = (weighted.sd(rate.change, weight, na.rm = TRUE))^2)


##PLOT - RATE OF CHANGE----

plot.ratechange <- Abun.fish.rate %>%
  ggplot(aes(x = bin_mean_FEve, y = rate.change)) +
  labs(title = "Temporal varation in biomass for NFL as a response to community functional evenness") +
  xlab("Functional Evenness") +
  ylab("Rate of biomass change") +
  geom_point() +
  geom_smooth(method = lm, aes(color = new_bin)) +
  theme_light()

plot.ratechange

#ylim(0,2000)
#facet_wrap(~year)



#have to scale the variables as they have different scales
abun.feve.bin.noNA$Z_FEve <- scale(abun.feve.bin.noNA$bin_mean_FEve)
abun.feve.bin.noNA$Z_sdbio <-
  scale(abun.feve.bin.noNA$bin_sd_biomass)



#Trying data out in a linear model
lm.test <- lm(Z_sdbio ~ Z_FEve, data = abun.feve.bin.noNA)
lm.test.resid <- rstandard(lm.test)

#testing model residuals against the random factors to determine if they need to be included
plot(
  lm.test.resid ~ as.factor(abun.feve.bin.noNA$stratum),
  xlab = "Stratum",
  ylab = "Standardized residuals"
)
abline(0, 0, lty = 2)

plot(
  lm.test.resid ~ as.factor(abun.feve.bin.noNA$new_bin),
  xlab = "Year",
  ylab = "Standardized residuals"
)
abline(0, 0, lty = 2)
#there is a lot of variation around the means, this indicates that these factors should be included

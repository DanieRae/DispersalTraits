#Statistics and response variables#

#Packages----
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
library(stringr)
library(rgdal)
library(rgeos)
library(proj4)
library(spdep)
library(mgcv)
library(viridis)

#NOTES ----
#With the clustering of the species here a a few responses that should be looked at#

#FOR how much over redundancy in each cluster is there? Does the pattern of over redundancy change spatially? How many FEs have just one species? Is there a clear pattern in what characteristics have only 1 species? What are the ecological services provided by FE groups with high FOR vs those with none? Do we see change over time (quick/slow/more or less diversity in different areas), larger home ranges, higher dispersal capabilitys in certian FEs?

#Include abundance distributions (biomass) across FE's, before and during collapse

#How does the FE diversity compare to hill numbers locally (alpha and gamma) gradients of beta accross the space

#need to call in the abundance data

#EFFECTIVE SPECIES FILLED ----
effective.species.filled <- fish.abun.complete  %>%
  #calculate the Hill number for each stratum in each year
  group_by(stratum, year_surv) %>%
  dplyr::summarize(effective_species = exp(diversity(final_biomass, "shannon")))

#ABUNDANCE DATA WRANGLING ----

##TOTAL BIOMASS---- 
#per strata/year
Abun.fish.strata <- fish.abun.complete %>%
  group_by(stratum, year_surv) %>%
  dplyr::summarize(group_biomass = sum(final_biomass))

fish.abun.complete.noweight <- select(fish.abun.complete, -c(weight))

##WIDE----
#abundance by site (strata,year)
Abun.fish.wide <-
  as.data.table (spread(fish.abun.complete.noweight,taxa_name, final_biomass))

Abun.fish.biomass <-
  select(Abun.fish.wide,-c ("stratum", "year_surv"))

Abun.fish.names <- sort(names(Abun.fish.biomass), decreasing = FALSE)

Abun.fish.biomass <- Abun.fish.biomass[, ..Abun.fish.names]


Abun.fish.biomass <- as.matrix(Abun.fish.biomass)


#FD --------

#Identity <- functcomp(fish.traits.40NA, Abun.fish.biomass)

# functional.diversity.FT <- dbFD(fish.traits.40NA, Abun.fish.biomass, corr = "lingoes")

##FD SAVED FILE/LOAD ----
#saveRDS(functional.diversity.FT, file = "FunctionalDiv.rds")

#functional.diversity.FT <- readRDS("FunctionalDiv.rds")
functional.diversity.FT.new <- readRDS("FunctionalDivNEW.rds")

year.strat <- select(Abun.fish.wide, c("stratum", "year_surv"))

##FEVE ----
FEve.comm <-
  cbind(
    functional.diversity.FT.new$FEve,
    stratum = year.strat$stratum,
    year = year.strat$year_surv
  )

FEve.comm <- as.data.frame(FEve.comm)

CWM.comm <-
  cbind(
    stratum = year.strat$stratum,
    functional.diversity.FT.new$CWM,
    year = year.strat$year_surv
  )

#BINS ----
#5 year averages for biomass, FEve, and effective species
tags <-
  c("1995-1999",
    "2000-2004",
    "2005-2009",
    "2010-2014",
    "2015-2017")

##ABUNDANCE ----
#mean and sd of biomass per bin

abun.bin <- Abun.fish.strata %>%
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


##FEVE ---- 
#mean and sd per bin

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

FEve.bin.sd <- FEve.bin %>%
  #sd of the binned years
  group_by(stratum, new_bin) %>%
  dplyr::summarise(bin_sd_FEve = sd(V1))

FEve.bin.FD <- merge(FEve.bin.mean, FEve.bin.sd)

## EFFECTIVE SPECIES ----
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

#SPECIES-DISPERSAL DIV ----
effective.species <-
  dplyr::rename(effective.species, year = year_surv)

effective.species.FEve <-
  merge (x = FEve.comm, y = effective.species)

effective.species.FEve$year <-
  as.factor(effective.species.FEve$year)

##PLOT - SPDIV/FEVE ----
#Effective species/FEve plot

plot.effective.FEve <- effective.species.FEve %>%
  filter (year %in% c(1995,2000,2005,2010)) %>%
  ggplot(aes(x = effective_species, y = V1)) +
  labs(title = "Relationship between FEve and Effective species") +
  xlab ("Effective Species Diveristy") +
  ylab ("Dispersal Diversity (Functional Evenness)") +
  theme_light() +
  geom_point() +
  geom_smooth(method = lm) +
  #facet_wrap( ~ year)+
  theme(legend.position = "none", 
        strip.text.x = element_text(size=14,vjust=0, face = "bold",
                                    family="Arial Narrow"),
        plot.title = element_text(lineheight = .8, size = 20),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) # title
       
#BINNED SPECIES-DISPERSAL DIV ----
binned.effective.species.FEve <-
  merge (x = FEve.bin.mean, y = effective.species.filled.bin.mean)

# GAM - SPDIV/FEVE BINNED ----

gam.spdiv.feve.binned <- gam(bin_mean_FEve ~ bin_mean_effectivesp, 
                             data = binned.effective.species.FEve,
                             method = "REML",
                             family = gaussian)
summary(gam.spdiv.feve.binned)

##PLOT - SPDIV/FEVE BINNED ----
#Effective species/FEve plot

plot.binned.effective.species.FEve <- binned.effective.species.FEve %>%
  #filter (year %in% c(1995,2000,2005,2010)) %>%
  ggplot(aes(x = bin_mean_effectivesp, y = bin_mean_FEve)) +
  labs(title = "Influence of Effective Species Diversity on Dispersal Diversity") +
  xlab ("Effective Species Diveristy") +
  ylab ("Dispersal Diversity (Functional Evenness)") +
  theme_light() +
  geom_point() +
  geom_smooth(method = lm)+
  #facet_wrap( ~ year)+
  theme(legend.position = "none", 
        strip.text.x = element_text(size=14,vjust=0, face = "bold",
                                    family="Arial Narrow"),
        plot.title = element_text(lineheight = .8, size = 20),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) # title

##PLOT - FEVE/TIME ----
#functional evenness by community over time

plot.FEve.time <- FEve.comm %>%
  #filter(stratum %in% c(201, 202, 203, 204, 205, 206, 207, 208, 209, 210)) %>%
  ggplot (aes(x = year, y = V1)) +
  geom_point(size = 3) +
  #geom_text()+
  theme_bw() +
  xlab("Year")+
  ylab("Dispersal Diversity")+
  ggtitle("Change in Local Community Dispersal Diversity over a 20 Year Period") +
  geom_smooth() +
  scale_color_grey() +
  facet_wrap_paginate( ~ stratum,
              nrow = 6,
              ncol = 6,
              page = 2)+
  theme(plot.title = element_text(lineheight = .8, size = 20, hjust = 0.5), # title
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #axis.text.x = element_blank(), # remove x axis labels
        #axis.text.y = element_blank(), # remove y axis labels  
        axis.ticks = element_blank(), # remove axis ticks
        panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(), # remove grid lines
        strip.text.x = element_text(size=14,vjust=0, face = "bold",
                                    family="Arial Narrow"),
        legend.position="none")


# should I do this all again in the clusters?


##FD/depth plot----

# FD.Depth.Year <- stratum.FE.depth %>%
#   filter(year_surv %in% c(1995, 2000, 2005, 2010, 2015, 2017)) %>%
#   ggplot(aes(x = depth, y = Unique_FE)) +
#   geom_point() +
#   geom_smooth() +
#   facet_wrap( ~ year_surv) +
#   labs(x = "Stratum Depth", y = "Fungtional Groups")
# 
# FEve.comm.depth <-
#   merge (x = FEve.comm, y = stratum.new, by = "stratum") %>%
#   sf::st_as_sf()
# 
# FEve.Depth.Year <- FEve.comm.depth %>%
#   filter(year %in% c(1995, 2000, 2005, 2010, 2015, 2017)) %>%
#   ggplot(aes(x = depth, y = V1)) +
#   geom_point() +
#   geom_smooth() +
#   facet_wrap( ~ year) +
#   labs(x = "Stratum Depth", y = "Fungtional Evenness")
# 

#MERGE FEVE AND STAB ----
abun.feve.bin <- merge (FEve.bin.FD,  abun.bin.biomass)
abun.feve.bin$new_bin <- as.factor(abun.feve.bin$new_bin)
abun.feve.bin$stratum <- as.factor(abun.feve.bin$stratum)

#Need to figure out how to select the first and last year in every bin
FEve.lag <- FEve.bin %>%
  #sd of the binned years
  filter(year %in% c(1995, 2000, 2005, 2010, 2015))

##PLOT - logSDbiomass----
plot.logSDbiomass <- abun.feve.bin %>%
  ggplot(aes(x = bin_sd_FEve, y = log(bin_sd_biomass))) +
  labs(title = "Temporal varation in biomass for NFL as a response to community functional evenness") +
  xlab ("Functional Evenness [sd(FEve)]") +
  ylab ("log(sd(biomass))") +
  geom_point() +
  geom_smooth(method = lm, aes(color = new_bin))

ggsave(
  "Biomass stability (sd) with FEve (sd).png",
  plot.logSDbiomass,
  width =  10,
  height = 10
)

##PLOT - Mean FEve with logsd Biomass ----
plot.meanfeve <- abun.feve.bin %>%
  filter(stratum %in% c(201:230)) %>%
  ggplot(aes(x = bin_mean_FEve, y = log(bin_sd_biomass))) +
  labs(title = "Temporal varation in biomass for NFL as a response to community functional evenness") +
  xlab ("Mean Functional Evenness") +
  ylab ("log(sd(biomass))") +
  geom_point(aes(color = new_bin)) +
  geom_path() +
  theme() +
  facet_wrap( ~ stratum)

ggsave(
  "Biomass stability (sd) with FEve (mean).png",
  plot.meanfeve,
  width =  10,
  height = 10
)


#RATE OF BIOMASS CHANGE ----

#what is the rate of change of biomass in a strata from year1-yearn
Abun.fish.rate <- abun.feve.bin %>%
  group_by(stratum) %>%
  #because we don't have individual taxa we can;t look at sd
  arrange(stratum, new_bin) %>%
  dplyr::mutate(rate.change = log(bin_mean_biomass / lag(bin_mean_biomass)))

#this could be used for the rate of compositional change but no with the current state of summed biomass accross the strata
#Abun.fish.rate.SD <-Abun.fish.rate %>%
#  group_by(stratum,year_surv)%>%
#  dplyr:: summarise(sd2.rate.change = (weighted.sd(rate.change, weight, na.rm = TRUE))^2)


##PLOT - RATE OF CHANGE----

plot.ratechange <- Abun.fish.rate %>%
  ggplot(aes(x = bin_mean_FEve, y = rate.change)) +
  labs(title = "Temporal varation in biomass for NFL as a response to community functional evenness") +
  xlab ("Functional Evenness") +
  ylab ("Rate of biomass change") +
  geom_point() +
  geom_smooth(method = lm, aes(color = new_bin)) +
  theme_light()

plot.ratechange

#ylim(0,2000)
#facet_wrap(~year)


#LM MODELS----
#too many rows with NAs to model
abun.feve.bin.noNA <- abun.feve.bin[complete.cases(abun.feve.bin), ]

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

###Lets try with glm #WHY YOU NO WORK
glm1 <- glmer(
  Z_sdbio ~ Z_FEve + (1 | stratum) + (1 | new_bin),
  data = abun.feve.bin.noNA,
  family = gaussian(link = "log"),
  start = 0
)

source(file = "functions/glmm_funs.R")

if (!require("coefplot2"))
  remotes::install_github(
    "palday/coefplot2",
    subdir = "pkg",
    upgrade = "always",
    quiet = TRUE
  )

library(coefplot2)
overdisp_fun(glm1)

summary(glm1)
tab_model(glm1)
sjPlot::plot_model(glm1)

effects.V1 <- effect(term = "V1", mod = glm1)
effects.V1 <- as.data.frame(effects.V1)

abun.feve.bin.noNA.fit <-
  filter(abun.feve.bin.noNA, stratum %in% c(201:202))


ggplot() +
  geom_point(data = rate.feve.noNA.fit, aes(x = V1, y = sd2.rate.change)) +
  geom_point(data = effects.V1, aes(x = V1, y = fit), color = "red") +
  geom_line(data = effects.V1, aes(x = V1, y = fit), color = "red") +
  geom_ribbon(
    data = effects.V1,
    aes (x = V1, ymin = lower, ymax = upper),
    alpha = 0.3,
    fill = "red"
  ) +
  facet_wrap( ~ stratum)

#GAM ----

##GAM1 scaled variables with the NA's removed from the data ----
gam1 <- gam(Z_sdbio ~ s(Z_FEve) ,
            data = abun.feve.bin.noNA, method = 'REML')
summary(gam1)
gam.check(gam1)

data_plotgam1 <-
  ggplot(data = abun.feve.bin.noNA, aes(y = Z_sdbio, x = Z_FEve)) +
  geom_point() +
  geom_line(aes(y = fitted(gam1)),
            colour = "red", size = 1.2) +
  #geom_line(aes(y = fitted(gam2)),
  #         colour = "blue", size = 1.2)+
  theme_bw()
data_plotgam1

##GAM2 unscaled with a log of the dependent variable  ----
gam2 <- gam(log(bin_sd_biomass) ~ s(bin_mean_FEve) ,
            data = abun.feve.bin.noNA,
            method = 'REML')
summary(gam2)

data_plotgam2 <-
  ggplot(data = abun.feve.bin.noNA, aes(y = log(bin_sd_biomass), x = bin_mean_FEve)) +
  geom_point() +
  geom_line(aes(y = fitted(gam2)),
            colour = "red", size = 1.2) +
  #geom_line(aes(y = fitted(gam2)),
  #         colour = "blue", size = 1.2)+
  theme_bw()+
  labs(title = "Marine Community Stability in Response to Dispersal Diversity") +
  xlab("Dispersal Diversity (Functional Evenness)")+
  ylab ("Stability (SD of Biomass) ")

data_plotgam2
## GAM3 noNA and no log ----
#By loging the dependent varable I achieve two 'normal' distributions for my data. this gives a linear regression. Would I go back to lmms instead?

gam3 <- gam(bin_sd_biomass ~ s(bin_mean_FEve) ,
            data = abun.feve.bin.noNA,
            method = 'REML')
summary(gam3)
#no different than the scaled version
data_plotgam3 <-
  ggplot(data = abun.feve.bin.noNA, aes(y = bin_sd_biomass, x = bin_mean_FEve)) +
  geom_point() +
  geom_line(aes(y = fitted(gam3)),
            colour = "red", size = 1.2) +
  #geom_line(aes(y = fitted(gam2)),
  #         colour = "blue", size = 1.2)+
  theme_bw()
data_plotgam3

AIC(gam1, gam2, gam3, gamm1)
#visually gam1 and gam3 look very similar but with the AIC gam 1 comes out as the most predictive. Still there isnt a great fit to the data. I think this will need the mixed model or a higher k -> changing K din't work


#GAMM ----
gamm1 <- gam(Z_sdbio ~ s(Z_FEve) +
               s(stratum, bs = "re") ,
             data = abun.feve.bin.noNA,
             method = 'REML')

gamm1_summary <- summary(gamm1)

gamm2 <- gam(
  log(bin_sd_biomass) ~ bin_mean_FEve +
    s(stratum, bs = "re", k = 12) +
    s(new_bin, bs ="re", k = 6),
  data = abun.feve.bin.noNA,
  method = 'REML'
)
summary(gamm2)

data_plotgamm2 <-
  ggplot(data = abun.feve.bin.noNA, aes(y = log(bin_sd_biomass), x = bin_mean_FEve)) +
  geom_point() +
  geom_line(aes(y = fitted(gamm2)),
            colour = "red", size = 1.2) +
  #geom_line(aes(y = fitted(gamm1)),
   #        colour = "blue", size = 1.2)+
  theme_bw()
data_plotgamm2


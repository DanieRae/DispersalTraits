#Statistics and response variables#

library(grid)
library(gridExtra)
library(lme4)
library(sjPlot)


#NOTES ----
#With the clustering of the species here a a few responses that should be looked at#

#FOR how much over redundancy in each cluster is there? Does the pattern of over redundancy change spatially? How many FEs have just one species? Is there a clear pattern in what characteristics have only 1 species? What are the ecological services provided by FE groups with high FOR vs those with none? Do we see change over time (quick/slow/more or less diversity in different areas), larger home ranges, higher dispersal capabilitys in certian FEs?

#Include abundance distributions (biomass) across FE's, before and during collapse 

#How does the FE diversity compare to hill numbers locally (alpha and gamma) gradients of beta accross the space

#need to call in the abundance data


#Abundance ----
##Per taxa by year within strata----
Abun.fish <-fish.abun.clean %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, vessel, trip, set, taxa_name)%>%
  dplyr::summarize(group_biomass = sum(density_kgperkm2))%>%
  #calculate average biomass per stratum for each functional group
  group_by(stratum, year_surv, taxa_name)%>%
  dplyr::summarize(group_biomass = mean(group_biomass))

##Total biomass per strata/year----
Abun.fish.strata <-Abun.fish %>%
  group_by(stratum, year_surv)%>%
  dplyr::summarize(group_biomass = mean(group_biomass))

#Wide version of the abundance by site (strata,year)----
Abun.fish.wide <- as.data.table (spread(Abun.fish, taxa_name, group_biomass, fill = 0))

Abun.fish.biomass <- select(Abun.fish.wide, -c ("stratum", "year_surv"))

Abun.fish.names <- sort(names(Abun.fish.biomass),decreasing = FALSE)

Abun.fish.biomass <- Abun.fish.biomass[, ..Abun.fish.names]


Abun.fish.biomass <- as.matrix(Abun.fish.biomass)

#Betadiversity----

beta.div.j<- lapply(fish.pa.pivot.mat, function(mat) {
  mat <- beta.div.comp(mat, coef= "J", quant = T)
  return(mat)
  }
)
  

#Function diversity measures --------

#Identity <- functcomp(fish.traits.40NA, Abun.fish.biomass) 

#functional.diversity.FT <- dbFD(fish.traits.40NA, 
                                #Abun.fish.biomass, corr = "lingoes")

saveRDS(functional.diversity.FT, file = "FunctionalDiv.rds")

functional.diversity.FT <-readRDS("FunctionalDiv.rds")

year.strat <- select(Abun.fish.wide, c("stratum", "year_surv"))

FEve.comm <- cbind(functional.diversity.FT$FEve, stratum = year.strat$stratum, year = year.strat$year_surv)

FEve.comm <- as.data.frame(FEve.comm)

CWM.comm <- cbind(stratum = year.strat$stratum, functional.diversity.FT$CWM, year = year.strat$year_surv)



##FD/FEve plot ----
#functional evenness by community over time 

FEve.comm %>%
  filter(stratum %in% c(201,202,203,204,205,206,207,208,209,210)) %>%
    ggplot (aes(x = year, y= V1))+  
    geom_point(size =3)+
    #geom_text()+
    theme_bw()+
    #xlab(paste("PCoA 1 -",fish.cdm.eig[1], "%", sep ="" ))+
    #ylab(paste("PCoA 2 -",fish.cdm.eig[2], "%", sep ="" ))+
    ggtitle("FEve community 201") +
    geom_smooth() + 
    theme(legend.position="none", asp=1) +
    scale_color_grey()+
    facet_wrap(~stratum)

# should I do this all again in the clusters?


##FD/depth plot----

FD.Depth.Year <-stratum.FE.depth %>%
  filter(year_surv %in% c(1995,2000,2005,2010,2015,2017)) %>% 
  ggplot(aes(x= depth, y= Unique_FE))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~year_surv)+
  labs(x= "Stratum Depth",y= "Fungtional Groups")

FEve.comm.depth <- merge (x = FEve.comm, y = stratum.new, by = "stratum")%>%
  sf::st_as_sf()

FEve.Depth.Year <-FEve.comm.depth %>%
  filter(year %in% c(1995,2000,2005,2010,2015,2017)) %>% 
  ggplot(aes(x= depth, y= V1))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~year)+
  labs(x= "Stratum Depth",y= "Fungtional Evenness")

#Abundance over time, sd by FEve----

#Abundance per community by year
Abun.fish.SD <-Abun.fish %>%
  #what is the SD of the biomass 
  group_by(stratum, year_surv)%>%
  dplyr::summarize(group_SD = sd(group_biomass)) %>%
  #log transform this to have a less scew?
  group_by(stratum, year_surv)%>%
  dplyr::mutate(log_SD = log(group_SD)) %>%
  #Need to take 1/SD so then plot with the FD 
  group_by(stratum, year_surv)%>%
  dplyr::mutate(one.group_SD = 1/(group_SD)) %>%
  #log of the 1/sd
  group_by(stratum, year_surv)%>%
  dplyr::mutate(log.one.group_SD = log(one.group_SD)) 
 
 
 
##sd plot----
Abun.fish.SD <- Abun.fish.SD %>% dplyr::rename(year = year_surv)
Abun.fish.FEve <- merge (x = FEve.comm, y = Abun.fish.SD)
Abun.fish.FEve$year <- as.factor(Abun.fish.FEve$year)

plot.logSD <- Abun.fish.FEve %>%
  filter (year %in% c(1995,2000,2005,2010,2015,2017)) %>%
  ggplot(aes(x= V1, y= log.one.group_SD))+
  labs(title = "Temporal varation in biomass for NFL as a response to community functional evenness") +
  xlab ("Functional Evenness") +
  ylab ("Log(SD(biomass)^-1)")+
  geom_point()+
  geom_smooth(method = lm, aes(color = year))
  #ylim(0,2000)
  #facet_wrap(~year)

#Rate of biomass change ----

Abun.fish.rate <-Abun.fish %>%
  #what is the rate of change of biomass in a strata from year1-yearn
  #Rate of biomass change from one year to another
  group_by(stratum,taxa_name)%>%
  arrange(stratum, year_surv)%>%
  dplyr::mutate(rate.change = log(group_biomass / lag(group_biomass))) 

Abun.fish.rate.SD <-Abun.fish.rate %>%
  group_by(stratum,year_surv)%>%
  dplyr:: summarise(sd2.rate.change = (sd(rate.change, na.rm = TRUE))^2)

#Trying by hand/fail  
#Sum of the rate of change, if 1 NA is present it is summing to NA  
#Abund.fish.mean.rate <-Abun.fish.rate%>% 
  ##dplyr::mutate (mean.rate.change.strata.year = mean(rate.change, na.rm = TRUE)) %>%
 # group_by(stratum, year_surv) %>%
  #dplyr::summarise(sum.rate.change.strata.year = sum(rate.change, na.rm = TRUE)) %>%
 # ungroup ()

##Rate of change plot----

Abun.fish.rate.SD <- Abun.fish.rate.SD %>% dplyr::rename(year = year_surv)
Abun.fish.rate.SD.FEve <- merge (x = FEve.comm, y = Abun.fish.rate.SD)
Abun.fish.rate.SD.FEve$year <- as.factor(Abun.fish.rate.SD.FEve$year)
Abun.fish.rate.SD.FEve$stratum <- as.factor(Abun.fish.rate.SD.FEve$stratum)

plot.SDratechange <- Abun.fish.rate.SD.FEve %>%
  filter (stratum %in% c(201:210)) %>%
  #filter (year %in% c(1996,2000,2005,2010,2015,2017)) %>%
  ggplot(aes(x= V1, y= sd2.rate.change))+
  labs(title = "Temporal varation in biomass for NFL as a response to community functional evenness") +
  xlab ("Functional Evenness") +
  ylab ("Rate of biomass change")+
  geom_point()+
  geom_smooth(method = lm, aes(color = stratum))+
  theme_light()
  #facet_wrap(~year)
  theme(legend.position = 'none')


plot.SDratechange
#ylim(0,2000)
#facet_wrap(~year)

##Model rate of biomass change----

###Trying with lm models, checking residuals----
#have to scale the variables as they have different scales
Abun.fish.rate.SD.FEve$Z_V1 <- scale(Abun.fish.rate.SD.FEve$V1)
Abun.fish.rate.SD.FEve$Z_sd2 <- scale(Abun.fish.rate.SD.FEve$sd2.rate.change )

#too many rows with NAs to model
rate.feve.noNA<-Abun.fish.rate.SD.FEve[complete.cases(Abun.fish.rate.SD.FEve),]

#Trying data out in a linear model
lm.test<-lm(Z_sd2 ~ Z_V1, data=rate.feve.noNA)
lm.test.resid <- rstandard(lm.test)

#testing model residuals against the random factors to determine if they need to be included
plot(lm.test.resid ~ as.factor(rate.feve.noNA$stratum),
     xlab = "Stratum", ylab = "Standardized residuals")
abline(0, 0, lty = 2)

plot(lm.test.resid ~ as.factor(rate.feve.noNA$year),
     xlab = "Year", ylab = "Standardized residuals")
abline(0, 0, lty = 2)
#there is a lot of variation around the means, this indicates that these factors should be included

###Lets try with glm
glm1 <- glmer(sd2.rate.change ~ V1 + stratum + (1|year),
     data=rate.feve.noNA, family = gaussian(link = "log"))

source(file = "functions/glmm_funs.R")

if (!require("coefplot2"))
  remotes::install_github("palday/coefplot2",
                          subdir = "pkg",
                          upgrade = "always",
                          quiet = TRUE)

library(coefplot2)
overdisp_fun(glm1)

summary(glm1)
tab_model(glm1)
sjPlot::plot_model(glm1)

effects.V1 <- effect(term = "V1", mod = glm1)
effects.V1 <- as.data.frame(effects.V1)

rate.feve.noNA.fit <- filter(rate.feve.noNA, stratum %in% c(201:202))


ggplot()+
  geom_point(data = rate.feve.noNA.fit, aes(x= V1, y= sd2.rate.change))+
  geom_point(data = effects.V1, aes(x= V1, y= fit), color = "red")+
  geom_line(data = effects.V1, aes(x= V1, y= fit), color = "red")+
  geom_ribbon(data = effects.V1, aes (x= V1, ymin= lower, ymax= upper), alpha= 0.3, fill = "red")+
  facet_wrap(~stratum)

#Effective species vs FEve----
effective.species <- effective.species %>% rename(year = year_surv)
effective.species.FEve <- merge (x = FEve.comm, y = effective.species)

##Effective species/FEve plot----

plot.effective.FEve <- effective.species.FEve %>%
  filter (year %in% c(1995,2000,2005,2010,2015,2017)) %>%
  ggplot(aes(x= effective_species, y= V1))+
  labs(title = "Relationship between FEve and Effective species") +
  xlab ("Effective species") +
  ylab ("Functional Evenness)")+
  geom_point()+
  geom_smooth(method = lm, aes(color = year))
#ylim(0,2000)



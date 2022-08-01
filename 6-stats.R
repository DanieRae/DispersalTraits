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


#Abundance ----
#
##Total biomass per strata/year----
Abun.fish.strata <-fish.abun.complete %>%
  group_by(stratum, year_surv)%>%
  dplyr::summarize(group_biomass = sum(final_biomass))

##Wide version of the abundance by site (strata,year)----
Abun.fish.wide <- as.data.table (spread(fish.abun.complete, taxa_name, final_biomass, fill = 0))

Abun.fish.biomass <- select(Abun.fish.wide, -c ("stratum", "year_surv","weight"))

Abun.fish.names <- sort(names(Abun.fish.biomass),decreasing = FALSE)

Abun.fish.biomass <- Abun.fish.biomass[, ..Abun.fish.names]


Abun.fish.biomass <- as.matrix(Abun.fish.biomass)

#Betadiversity----

beta.div.j<- lapply(fish.pa.pivot.mat, function(mat) {
  mat <- beta.div.comp(mat, coef= "J", quant = T)
  return(mat)
  }
)


#Effective species vs FEve----
effective.species <- dplyr::rename(effective.species, year = year_surv)
effective.species.FEve <- merge (x = FEve.comm, y = effective.species)
effective.species.FEve$year <- as.factor(effective.species.FEve$year)

##Effective species/FEve plot----

plot.effective.FEve <- effective.species.FEve %>%
  filter (year %in% c(1995,2000,2005,2010,2015,2017)) %>%
  ggplot(aes(x= effective_species, y= V1))+
  labs(title = "Relationship between FEve and Effective species") +
  xlab ("Effective species") +
  ylab ("Functional Evenness")+
  geom_point()+
  geom_smooth(method = lm, aes(color = year))+
  facet_wrap(~year)
#ylim(0,2000)  

#Function diversity measures --------

Identity <- functcomp(fish.traits.40NA, Abun.fish.biomass) 

#functional.diversity.FT <- dbFD(fish.traits.40NA, Abun.fish.biomass, corr = "lingoes")

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

##Bins for 5 year biomass table and FEve table----
tags <- c("1995-1999","2000-2004","2005-2009","2010-2014","2015-2017")

abun.bin <-Abun.fish.strata %>%
  mutate(new_bin = case_when(
    year_surv <=1999 ~ tags[1],
    year_surv >1999 & year_surv <2005 ~ tags[2],
    year_surv >=2005 & year_surv <2010 ~ tags[3],
    year_surv >=2010 & year_surv <2015 ~ tags[4],
    year_surv >=2015 & year_surv <=2017 ~ tags[5]))


FEve.bin <-FEve.comm %>%
  mutate(new_bin = case_when(
    year <=1999 ~ tags[1],
    year >1999 & year <2005 ~ tags[2],
    year >=2005 & year <2010 ~ tags[3],
    year >=2010 & year <2015 ~ tags[4],
    year >=2015 & year <=2017 ~ tags[5]))


##Abundance mean and sd of biomass per bin----
abun.bin.mean <-abun.bin %>%
  #mean biomass for the binned years
  group_by(stratum, new_bin)%>%
  dplyr::summarize(bin_mean_biomass = mean(group_biomass))

abun.bin.sd <-abun.bin %>%
  #sd of the binned years
  group_by(stratum, new_bin)%>%
  dplyr::summarise(bin_sd_biomass = sd(group_biomass))
  
abun.bin.biomass <- merge(abun.bin.mean, abun.bin.sd) 

##FEve mean and sd per bin----
FEve.bin.mean <-FEve.bin %>%
  #mean biomass for the binned years
  group_by(stratum, new_bin)%>%
  dplyr::summarize(bin_mean_FEve = mean(V1))

FEve.bin.sd <-FEve.bin %>%
  #sd of the binned years
  group_by(stratum, new_bin)%>%
  dplyr::summarise(bin_sd_FEve = sd(V1))

FEve.bin.FD <- merge(FEve.bin.mean,FEve.bin.sd)

abun.feve.bin <- merge ( FEve.bin.FD,  abun.bin.biomass)
abun.feve.bin$new_bin <- as.factor(abun.feve.bin$new_bin)
abun.feve.bin$stratum <- as.factor(abun.feve.bin$stratum)


#Need to figure out how to select the first and last year in every bin 
FEve.lag <-FEve.bin %>%
  #sd of the binned years
  filter(year %in% c(1995,2000,2005,2010,2015))
 
##sd plot logSDbiomass----
plot.logSDbiomass <- abun.feve.bin %>%
  ggplot(aes(x= bin_sd_FEve, y= log(bin_sd_biomass)))+
  labs(title = "Temporal varation in biomass for NFL as a response to community functional evenness") +
  xlab ("Functional Evenness [sd(FEve)]") +
  ylab ("log(sd(biomass))")+
  geom_point()+
  geom_smooth(method = lm, aes(color = new_bin))
  
ggsave("Biomass stability (sd) with FEve (sd).png", plot.logSDbiomass, width =  10, height = 10)

##Mean FEve plot with logsd Biomass ----
plot.meanfeve <- abun.feve.bin %>%
  filter(stratum %in% c(201:230))%>%
  ggplot(aes(x= bin_mean_FEve, y= log(bin_sd_biomass)))+
  labs(title = "Temporal varation in biomass for NFL as a response to community functional evenness") +
  xlab ("Mean Functional Evenness") +
  ylab ("log(sd(biomass))")+
  geom_point(aes(color = new_bin))+
  geom_path()+
  theme()+
  facet_wrap(~stratum)

ggsave("Biomass stability (sd) with FEve (mean).png", plot.meanfeve, width =  10, height = 10)


#Rate of biomass change ----

#what is the rate of change of biomass in a strata from year1-yearn
Abun.fish.rate <- abun.feve.bin %>%  
  group_by(stratum)%>%
  #because we don't have individual taxa we can;t look at sd
  arrange(stratum, new_bin)%>%
  dplyr::mutate(rate.change = log(bin_mean_biomass / lag(bin_mean_biomass))) 

#this could be used for the rate of compositional change but no with the current state of summed biomass accross the strata
#Abun.fish.rate.SD <-Abun.fish.rate %>%
#  group_by(stratum,year_surv)%>%
#  dplyr:: summarise(sd2.rate.change = (weighted.sd(rate.change, weight, na.rm = TRUE))^2)


##Rate of change plot----

plot.ratechange <- Abun.fish.rate %>%
  ggplot(aes(x= bin_mean_FEve, y= rate.change))+
  labs(title = "Temporal varation in biomass for NFL as a response to community functional evenness") +
  xlab ("Functional Evenness") +
  ylab ("Rate of biomass change")+
  geom_point()+
  geom_smooth(method = lm, aes(color = new_bin))+
  theme_light()

plot.ratechange

#ylim(0,2000)
#facet_wrap(~year)

##Model rate of biomass change----

###Trying with lm models, checking residuals----
#too many rows with NAs to model
abun.feve.bin.noNA<-abun.feve.bin[complete.cases(abun.feve.bin),]

#have to scale the variables as they have different scales
abun.feve.bin.noNA$Z_FEve <- scale(abun.feve.bin.noNA$bin_mean_FEve)
abun.feve.bin.noNA$Z_sdbio <- scale(abun.feve.bin.noNA$bin_sd_biomass)



#Trying data out in a linear model
lm.test<-lm(Z_sdbio ~ Z_FEve, data=abun.feve.bin.noNA)
lm.test.resid <- rstandard(lm.test)

#testing model residuals against the random factors to determine if they need to be included
plot(lm.test.resid ~ as.factor(abun.feve.bin.noNA$stratum),
     xlab = "Stratum", ylab = "Standardized residuals")
abline(0, 0, lty = 2)

plot(lm.test.resid ~ as.factor(abun.feve.bin.noNA$new_bin),
     xlab = "Year", ylab = "Standardized residuals")
abline(0, 0, lty = 2)
#there is a lot of variation around the means, this indicates that these factors should be included

###Lets try with glm #WHY YOU NO WORK
glm1 <- glmer(Z_sdbio ~ Z_FEve + (1|stratum) + (1|new_bin),
     data=abun.feve.bin.noNA, family = gaussian(link = "log"), start = 0)

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

##Trying gams time ----

#gam1 will used the scaled variables with the NA's removed from the data
gam1 <- gam(Z_sdbio ~ s(Z_FEve) ,
            data=abun.feve.bin.noNA, method = 'REML')
summary(gam1)
gam.check(gam1)

data_plotgam1 <- ggplot(data = abun.feve.bin.noNA, aes(y = Z_sdbio, x = Z_FEve)) + 
  geom_point() +
  geom_line(aes(y = fitted(gam1)),
            colour = "red", size = 1.2) + 
  #geom_line(aes(y = fitted(gam2)),
   #         colour = "blue", size = 1.2)+
  theme_bw()
data_plotgam1

#gam2 is using the data unscaled but with a log of the dependant varable to visualize it better
gam2 <- gam(log(bin_sd_biomass) ~ s(bin_mean_FEve) ,
            data=abun.feve.bin.noNA, method = 'REML')
summary(gam2)

data_plotgam2 <- ggplot(data = abun.feve.bin.noNA, aes(y = log(bin_sd_biomass), x = bin_mean_FEve)) + 
  geom_point() +
  geom_line(aes(y = fitted(gam2)),
            colour = "red", size = 1.2) + 
  #geom_line(aes(y = fitted(gam2)),
  #         colour = "blue", size = 1.2)+
  theme_bw()
data_plotgam2
#By loging the dependent varable I achieve two 'normal' distributions for my data. this gives a linear regression. Would I go back to lmms instead?

gam3 <- gam(bin_sd_biomass ~ s(bin_mean_FEve) ,
            data=abun.feve.bin.noNA, method = 'REML')
summary(gam3)
#no different than the scaled version
data_plotgam3 <- ggplot(data = abun.feve.bin.noNA, aes(y = bin_sd_biomass, x = bin_mean_FEve)) + 
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
               s(stratum, bs= "re") ,
            data=abun.feve.bin.noNA, method = 'REML')

gamm1_summary <-summary(gamm1)

gamm2 <- gam(log(bin_sd_biomass) ~ bin_mean_FEve +
               s(stratum, bs= "re")+
               new_bin,
            data=abun.feve.bin.noNA, method = 'REML')
summary(gamm2)

data_plotgamm2 <- ggplot(data = abun.feve.bin.noNA, aes(y = log(bin_sd_biomass), x = bin_mean_FEve)) + 
  geom_point() +
  geom_line(aes(y = fitted(gamm2)),
            colour = "red", size = 1.2) + 
  #geom_line(aes(y = fitted(gam2)),
  #         colour = "blue", size = 1.2)+
  theme_bw()
data_plotgamm2


#MRF----
# stratum.shpfile <- st_read("all_strata.shp") %>% 
#   # Validating geometries to make operations on them 
#   st_make_valid() %>%  
#   dplyr::mutate(stratum_id_base = str_sub(stratum, 1, 3)) %>% 
#   group_by(stratum_id_base) %>% 
#   dplyr::summarise(geometry = st_union(geometry))%>%
#   dplyr::rename(stratum = stratum_id_base)

##### Depth Data ####
# Created depth data to be added to the shapefile
# {
#   strata.depth <- as.data.frame(stratum.shpfile$stratum)
#   strata.depth$depth <- seq(1, 345, 1)
#   colnames(strata.depth) <- c("stratum", "depth")
# }

##### Combine Data ####
stratum.shpfile.abun <- merge(stratum.shpfile, abun.feve.bin.noNA,by = "stratum" )
# 
# class(stratum.shpfile)

##### SpatialPolygonsDataFrame ####
shp <- as(stratum.shpfile.abun, 'Spatial')

# Create dataframe from data layer in shapefile 
df <- droplevels(as(shp, 'data.frame'))

## project data
utm.proj <- "+proj=utm +zone=20 +south ellps=WGS84 +datum=WGS84"
shp <- spTransform(shp, CRS(utm.proj))  
shpf <- fortify(shp, region = "stratum")

#Make sure stratum is set to factor or else!
df$stratum <- as.factor(df$stratum)

#Create neighbours list
nb <- poly2nb(shp, row.names = df$stratum)
names(nb) <- attr(nb, "region.id")
str(nb[1:6]) 

#Test using mgcv example
ctrl <- gam.control(nthreads = 6) 

b <-
  gam(
    log(bin_sd_biomass) ~ s(stratum, bs = "mrf", xt = list(nb = nb)),
    data = abun.feve.bin.noNA,
    method = "REML",
    control = ctrl,
    family = gaussian()
  )

plot(b, scheme = 1)


#Fonyas example using gavin simpsons tutorial
ctrl <- gam.control(nthreads = 6) # use 6 parallel threads, reduce if fewer physical CPU cores
m1 <- gam(depth ~ s(stratum, bs = 'mrf', xt = list(nb = nb)), # define MRF smooth
          data = df,
          method = 'REML', # fast version of REML smoothness selection
          control = ctrl,
          family = gaussian()) # fit gaussian


df <- transform(df,
                mrfFull= predict(m1, type = 'response'))

## merge data with fortified shapefile
mdata <- left_join(shpf, df, by = c('id' = 'stratum'))

theme_map <- function(...) {
  theme_minimal() +
    theme(...,
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank())
}

myTheme <- theme_map(legend.position = 'bottom')
myScale <- scale_fill_viridis(name = '%', option = 'plasma',
                              limits = c(1, 345),
                              labels = function(x) x,
                              guide = guide_colorbar(direction = "horizontal",
                                                     barheight = unit(2, units = "mm"),
                                                     barwidth = unit(75, units = "mm"),
                                                     title.position = 'left',
                                                     title.hjust = 0.5,
                                                     label.hjust = 0.5))
myLabs <- labs(x = NULL, y = NULL, title = 'Strata MRF')

ggplot(mdata, aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = mrfFull)) +
  geom_path(col = 'black', alpha = 0.5, size = 0.1) +
  coord_equal() +
  myTheme + myScale + myLabs
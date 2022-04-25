#Statistics and response variables#

#need to call in the abundance data

#Abundance per community by year
Abun.fish <-fish.abun.clean %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, vessel, trip, set, taxa_name)%>%
  dplyr::summarize(group_biomass = sum(density_kgperkm2))%>%
  #calculate average biomass per stratum for each functional group
  group_by(stratum, year_surv, taxa_name)%>%
  dplyr::summarize(group_biomass = mean(group_biomass))

Abun.fish.wide <- as.data.table (spread(Abun.fish, taxa_name, group_biomass, fill = 0))

Abun.fish.biomass <- select(Abun.fish.wide, -c ("stratum", "year_surv"))

Abun.fish.names <- sort(names(Abun.fish.biomass),decreasing = FALSE)

Abun.fish.biomass <- Abun.fish.biomass[, ..Abun.fish.names]


Abun.fish.biomass <- as.matrix(Abun.fish.biomass)

#------ Function diversity measures --------

Identity <- functcomp(fish.traits.40NA, 
                      Abun.fish.biomass) 

functional.diversity.FT <- dbFD(fish.traits.40NA, 
                                Abun.fish.biomass, corr = "lingoes")

year.strat <- select(Abun.fish.wide, c("stratum", "year_surv"))

FEve.comm <- cbind(functional.diversity.FT$FEve, stratum = year.strat$stratum, year = year.strat$year_surv)

FEve.comm <- as.data.frame(FEve.comm)

CWM.comm <- cbind(functional.diversity.FT$CWM, stratum = year.strat$stratum, year = year.strat$year_surv)

#Some plots to visualize the functional evenness by community over time. 

FEve.comm %>%
  filter(stratum %in% c(201,202,203,204,205,206,207,208,209,210)) %>%
    ggplot (aes(x = year, y= V1))+  
    geom_point(size =3)+
    #geom_text()+
    theme_bw()+
    #xlab(paste("PCoA 1 -",fish.cdm.eig[1], "%", sep ="" ))+
    #ylab(paste("PCoA 2 -",fish.cdm.eig[2], "%", sep ="" ))+
    ggtitle("FEve community 201") +
    geom_path() + 
    theme(legend.position="none", asp=1) +
    scale_color_grey()+
    facet_wrap(~stratum)

# should I do this all again in the clusters?




#With the clustering of the species here a a few responses that should be looked at#

#FOR how much over redundancy in each cluster is there? Does the pattern of over redundancy change spatially? How many FEs have just one species? Is there a clear pattern in what characteristics have only 1 species? What are the ecological services provided by FE groups with high FOR vs those with none? Do we see change over time (quick/slow/more or less diversity in different areas), larger home ranges, higher dispersal capabilitys in certian FEs?

#Include abundance distributions (biomass) across FE's, before and during collapse 

#How does the FE diversity compare to hill numbers locally (alpha and gamma) gradients of beta accross the space


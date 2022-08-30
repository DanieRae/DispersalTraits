#Install and load packages, istallation shouldn't be necessary again----

# install.packages("devtools")
# install_github("vqv/ggbiplot")
# install.packages("reshape")

library(dplyr)
library(tidyr)
library(stringr)
library(devtools)
library(ggbiplot)
library(mltools)
library(data.table)
library(vegan)
library(ade4)
library(adespatial)
library(reshape)
library(skimr)

# PART ONE - TRAIT DATA ----

##Load/clean fish trait data----
raw.fish.traits <-
  read.csv("AtlanticFishTraits.csv",
           fileEncoding = "UTF-8-BOM",
           stringsAsFactors = TRUE)


#need to cut out all the dead space after the last species, not sure why that was there but bye bye#
fish.traits <- raw.fish.traits %>% droplevels()

#Making species name the row names and shortening them - might want to remove this as it is causing problems with the visuals on the plots#
rownames(fish.traits) <- fish.traits[, 1]
fish.traits[, 1] <- NULL

old_rownames <- rownames(fish.traits)
new_rownames <- unlist(lapply(lapply(str_split(old_rownames, " "),
                                     FUN = str_sub, 1, 3), paste0, collapse = "_"))
new_rownames[new_rownames == "Myo_sco"] <-
  c("Myo_sco_1", "Myo_sco_2")
rownames(fish.traits) <- sort(new_rownames, decreasing = FALSE)

#quick mutate to remove negative values in these columns#
fish.traits <- fish.traits %>%
  mutate(
    Fresh = abs(Fresh),
    Brack = abs(Brack),
    Saltwater = abs(Saltwater)
  )

##Data cleaning of unnecessary columns, such as references----
fish.traits <-
  select(
    fish.traits,
    -c (
      "Genus",
      "Subfamily",
      "Family",
      "Order",
      "Class",
      "Phylum",
      "Atlantic..Northwest",
      "Common",
      "DemersPelagRef",
      "MigratRef",
      "DepthRangeRef",
      "DepthComRef",
      "DepthRangeComShallow",
      "DepthRangeComDeep",
      "VerticalMove",
      "LTypeMaxM",
      "LTypeComM",
      "LTypeComF",
      "VerticalMoveRef",
      "LatMin",
      "LatMax",
      "LatRef",
      "LongevityWildRef",
      "MaxLengthRef",
      "MaxWeightRef",
      "CommonLengthRef",
      "MainCatchingMethod",
      "II",
      "MSeines",
      "MGillnets",
      "MCastnets",
      "MTraps",
      "MTrawls",
      "MDredges",
      "MLiftnets",
      "MHooksLines",
      "MSpears",
      "MOther",
      "MobilityRef",
      "FecundityRef",
      "SpawnRef",
      "LarvalRef",
      "LLRef",
      "LVRef",
      "PLDRef",
      "ADRef",
      "RaftingRef",
      "SchoolingRef",
      "Notes",
      "Comments",
      "Website.Link"
    )
  )

#Removing rows with missing larval strategy, this is considered an important trait that needs to be available for all species
fish.traits <- fish.traits[!is.na(fish.traits$LarvalStrategy), ]

#skimming data for structure
skim.tot <- skim(fish.traits)

# PART TWO - FISH ABUNDANCE DATA ----

##Load/clean fish abundance data----
raw.fish.abun <-
  read.csv("NL_Biomass.csv", fileEncoding = "UTF-8-BOM")

#I need to make the taxa names match with the other data set
## this doesn't appear to be used
rownames.fish.traits <- rownames(fish.traits)

#This takes the species list and puts it into a character#
species_name <- raw.fish.abun$taxa_name

#We are taking the character to then change the species names from all CAPS to just having the first letter in CAP, returns better data#
species_name <- str_to_sentence(species_name)

#Shortening the names to match
species_name[species_name == "Myoxocephalus scorpioides"] <-
  c("Myoxocephalus scorpioides 1")
species_name[species_name == "Myoxocephalus scorpius"] <-
  c("Myoxocephalus scorpius 2")


new_species_name <-
  unlist(lapply(lapply(str_split(species_name, " "),
                       FUN = str_sub, 1, 3), paste0, collapse = "_"))

fish.abun <- select(raw.fish.abun,-c ("taxa_name"))

fish.abun <- cbind(fish.abun, taxa_name = new_species_name)

#Need to cut out all the species not included in the main data
# %in% use for comparing two vectors of unequal length#

fish.abun.clean <-
  fish.abun %>% filter(taxa_name %in% rownames.fish.traits)


#table of minimun biomass values
min_biomass <- fish.abun.clean %>%
  group_by(taxa_name) %>%
  dplyr::summarise(min_biomass = min(density_kgperkm2))

fish.abun.clean.subset <- fish.abun.clean %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, vessel, trip, set, taxa_name) %>%
  dplyr::summarize(group_biomass = sum(density_kgperkm2)) %>%
  group_by(stratum, year_surv, taxa_name) %>%
  dplyr::summarize(group_biomass = mean(group_biomass)) %>%
  ungroup()

fish.abun.clean.complete <- fish.abun.clean.subset %>%
  #filter(year_surv == 1996)%>%
  #select(taxa_name, stratum, year_surv)%>%
  tidyr::complete(nesting(stratum, year_surv), taxa_name)


fish.abun.complete <- fish.abun.clean.complete %>%
  left_join(min_biomass) %>%
  mutate(final_biomass = ifelse(is.na(group_biomass), min_biomass, group_biomass)) %>%
  mutate(weight = ifelse(is.na(group_biomass), 0, 1)) %>%
  select(-c("group_biomass", "min_biomass"))


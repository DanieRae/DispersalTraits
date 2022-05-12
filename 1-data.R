#Install and load packages, istallation shouldn't be necessary again# 

# install.packages("devtools")
# install_github("vqv/ggbiplot")
#install.packages("reshape")


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
#library(reshape)

#PART ONE#
#load full data and view it#

raw.fish.traits <- read.csv("AtlanticFishTraits.csv", fileEncoding ="UTF-8-BOM", stringsAsFactors = TRUE)

#this was recommended to keep data in an easier format for analysis#
#atlantic.fish.traits <- as_tibble(atlantic.fish.traits)# this is messing with something so im removing it for now

#need to cut out all the dead space after the last species, not sure why that was there but bye bye#
atlantic.fish.traits <- raw.fish.traits %>% droplevels()

#Making species name the row names and shortening them - might want to remove this as it is causing problems with the visuals on the plots#

rownames(atlantic.fish.traits) <- atlantic.fish.traits[,1]
atlantic.fish.traits[,1] <- NULL


old_rownames <- rownames(atlantic.fish.traits)
new_rownames <- unlist(lapply(lapply(str_split(old_rownames, " "), 
                                      FUN = str_sub, 1, 3), paste0, collapse = "_"))
new_rownames[new_rownames == "Myo_sco"] <- c("Myo_sco_1", "Myo_sco_2")
rownames(atlantic.fish.traits) <- sort(new_rownames, decreasing = FALSE)
 

#quick mutate to remove negative values in these columns#
atlantic.fish.traits <- atlantic.fish.traits %>%
  mutate(Fresh = abs(Fresh),
         Brack = abs(Brack),
         Saltwater = abs(Saltwater))


#Data cleaning of unnecessary columns, such as references#

fish.traits <- select(atlantic.fish.traits, -c ("Genus","Subfamily", "Family", "Order", "Class", "Phylum","Atlantic..Northwest","Common", "DemersPelagRef", "MigratRef", "DepthRangeRef", "DepthComRef","LTypeMaxM","LTypeComM","LTypeComF", "VerticalMoveRef", "LatMin", "LatMax", "LatRef", "LongevityWildRef", "MaxLengthRef","MaxWeightRef", "CommonLengthRef", "MainCatchingMethod", "II", "MSeines", "MGillnets", "MCastnets", "MTraps", "MTrawls", "MDredges", "MLiftnets", "MHooksLines", "MSpears", "MOther", "MobilityRef", "FecundityRef", "SpawnRef", "LarvalRef", "LLRef", "LVRef", "PLDRef", "ADRef", "RaftingRef", "SchoolingRef", "Notes", "Comments", "Website.Link"))

#View current data format#

fish.traits <- fish.traits[!is.na(fish.traits$LarvalStrategy),]


skim.tot <- skim(fish.traits)


#PART TWO#
#Load abundance data#

fish.abun <- read.csv("NL_Biomass.csv", fileEncoding ="UTF-8-BOM")


#I need to make the taxa names match with the other data set

rownames.fish.traits <- rownames(fish.traits)

#This takes the species list and puts it into a character#
species_name <- fish.abun$taxa_name

#We are taking the character to then change the species names from all CAPS to just having the first letter in CAP, returns better data#
species_name <- str_to_sentence(species_name)

#Shortening the names to math
species_name[species_name == "Myoxocephalus scorpioides"] <- c("Myoxocephalus scorpioides 1")
species_name[species_name == "Myoxocephalus scorpius"] <- c("Myoxocephalus scorpius 2")


new_species_name <- unlist(lapply(lapply(str_split(species_name, " "), 
                                     FUN = str_sub, 1, 3), paste0, collapse = "_"))

fish.abun <- select(fish.abun, -c ("taxa_name"))

fish.abun <- cbind(fish.abun, taxa_name = new_species_name)

#Need to cut out all the species not included in the main data
# %in% use for comparing two vectors of unequal length#

fish.abun.clean <- fish.abun %>% filter(taxa_name %in% rownames.fish.traits)

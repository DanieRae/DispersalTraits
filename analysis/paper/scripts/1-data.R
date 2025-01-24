# Install and load packages ----

install.packages("dplyr")
install.packages("tidyr")
install.packages("stringr")
install.packages("devtools")
install.packages("mltools")
install.packages("data.table")
install.packages("vegan")
install.packages("ade4")
install.packages("adespatial")
install.packages("reshape")
install.packages("skimr")
install.packages("sf")
install.packages("here")


library(dplyr)
library(tidyr)
library(stringr)
library(devtools)
library(mltools)
library(data.table)
library(vegan)
library(ade4)
library(adespatial)
library(reshape)
library(skimr)
library(sf)
library(here)

# PART ONE - Trait data ----

##Load and clean fish trait data----
raw.fish.traits <-
  read.csv(here("analysis", "data", "raw_data", "tabular",
                "AtlanticFishTraits.csv"),
           fileEncoding = "UTF-8-BOM",
           stringsAsFactors = TRUE)


##Remove spaces ----
fish.traits <- raw.fish.traits %>% droplevels()

#making species name the row names and shortening them #
rownames(fish.traits) <- fish.traits[, 1]
fish.traits[, 1] <- NULL

old_rownames <- rownames(fish.traits)
new_rownames <- unlist(lapply(lapply(str_split(old_rownames, " "),
                                     FUN = str_sub, 1, 3), paste0, collapse = "_"))
new_rownames[new_rownames == "Myo_sco"] <-
  c("Myo_sco_1", "Myo_sco_2")
rownames(fish.traits) <- sort(new_rownames, decreasing = FALSE)

#Mutate to remove negative values in these columns#
fish.traits <- fish.traits %>%
  mutate(
    Fresh = abs(Fresh),
    Brack = abs(Brack),
    Saltwater = abs(Saltwater)
  )

#data cleaning of unnecessary columns, ex. references#
fish.traits <-
  select(
    fish.traits,
    -c(
      "Genus",
      "Subfamily",
      "Family",
      "Order",
      "Class",
      "Phylum",
      "Atlantic..Northwest",
      "Common",
      "DemersPelagRef",
      "Mobility",
      "Saltwater",
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

#removing rows(species) with missing larval strategy#
fish.traits <- fish.traits[!is.na(fish.traits$LarvalStrategy), ]


# PART TWO - Fish abundance data ----

##Load and clean fish abundance data----
raw.fish.abun <-
  read.csv(here("analysis", "data", "raw_data", "tabular",
                "NL_Biomass.csv"), fileEncoding = "UTF-8-BOM")

#This takes the species list and puts it into a character#
species_name <- raw.fish.abun$taxa_name

#Change the species names from all CAPS to just having the first letter in CAP#
species_name <- str_to_sentence(species_name)

#Shortening the names to match#
species_name[species_name == "Myoxocephalus scorpioides"] <-
  c("Myoxocephalus scorpioides 1")
species_name[species_name == "Myoxocephalus scorpius"] <-
  c("Myoxocephalus scorpius 2")


new_species_name <-
  unlist(lapply(lapply(str_split(species_name, " "),
                       FUN = str_sub, 1, 3), paste0, collapse = "_"))

fish.abun <- select(raw.fish.abun, -c("taxa_name"))

fish.abun <- cbind(fish.abun, taxa_name = new_species_name)

#Need to cut out all the species not included in the main data

#make the taxa names match with the other data set
rownames.fish.traits <- rownames(fish.traits)

# %in% use for comparing two vectors of unequal length#

fish.abun.filtered <-
  fish.abun %>% filter(taxa_name %in% rownames.fish.traits)

## abundance data not filled ----
# We first summarize the biomass data to get the mean biomass for each year and stratum
fish.abun.clean <- fish.abun.filtered %>%
  #first find, for each trawl and each functional group, the total biomass of that group in that trawl
  group_by(stratum, year_surv, vessel, trip, set, taxa_name) %>%
  dplyr::summarize(group_biomass = sum(density_kgperkm2)) %>%
  #calculate average biomass per stratum for each functional group
  group_by(stratum, year_surv, taxa_name) %>%
  dplyr::summarize(group_biomass = mean(group_biomass)) %>%
  ungroup()

## filling with zero ----
#Needs to be filled with 0 to determine the change in biomass from year to year, when species were not detected in a year they are given an abundance of 0
fish.abun.complete <- fish.abun.clean %>%
  #filter(year_surv == 1996)%>%
  #select(taxa_name, stratum, year_surv)%>%
  tidyr::complete(nesting(stratum, year_surv), taxa_name)

#Fill data frame NA with 0
fish.abun.complete[is.na(fish.abun.complete)] <- 0

# PART THREE - Stratum shapefile ------

# read all files of the shapefile
stratum.shpfile.raw <- st_read(here("analysis", "data", "raw_data", "spatial",
                                    "strata"))

#strata geometry, this joins the multiple polygons
stratum.shpfile <- stratum.shpfile.raw %>%
  # Validating geometries to make operations on them
  st_make_valid() %>%
  dplyr::mutate(stratum_id_base = str_sub(stratum, 1, 3)) %>%
  group_by(DIV,stratum_id_base) %>%
  dplyr::summarise(geometry = st_union(geometry)) %>%
  dplyr::rename(stratum = stratum_id_base) %>%
  filter(DIV != "2H", DIV != "2G") %>%
  filter(stratum != 613, stratum != 776, stratum != 777, stratum != 778)

#need to filter the adjusted strata, this will be used for the MRF on page 7
strata.keep <- stratum.shpfile$stratum

#this adjusted strata has been altered so that the margins of neighboring strata properly connect, to be used for the MRF
stratum.shpfile.adjusted <- st_read(here("analysis", "data", "derived_data", "spatial",
                                         "strata_adjusted")) %>%
  # Validating geometries to make operations on them
  st_make_valid() %>%
  dplyr::mutate(stratum_id_base = str_sub(stratum, 1, 3)) %>%
  group_by(stratum_id_base) %>%
  dplyr::summarise(geometry = st_union(geometry)) %>%
  dplyr::rename(stratum = stratum_id_base) %>%
  filter(stratum %in% strata.keep)

##Stratum Depth  ----
stratum.depth <- st_read("stratum_depth.csv")

stratum.depth$depth.ave <- as.numeric(stratum.depth$depth.ave)

#END#

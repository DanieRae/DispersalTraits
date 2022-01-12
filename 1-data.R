#Install and load packages, istallation shouldn't be necessary again# 

# install.packages("devtools")
# install_github("vqv/ggbiplot")


library(dplyr)
library(tidyr)
library(stringr)
library(FD)
library(gawdis)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(cluster)
library(mltools)
library(data.table)


#load full data and view it#

atlantic.fish.traits <- read.csv("AtlanticFishTraits.csv", fileEncoding ="UTF-8-BOM", stringsAsFactors = TRUE)

#this was recommended to keep data in an easier format for analysis#
#atlantic.fish.traits <- as_tibble(atlantic.fish.traits)# this is messing with something so im removing it for now

#need to cut out all the dead space after the last species, not sure why that was there but bye bye#
atlantic.fish.traits <- atlantic.fish.traits %>% slice(1:117)

#Making species name the row names and shortening them - might want to remove this as it is causing problems with the visuals on the plots#

rownames(atlantic.fish.traits) <- atlantic.fish.traits[,1]
atlantic.fish.traits[,1] <- NULL


old_rownames <- rownames(atlantic.fish.traits)
new_rownames <- unlist(lapply(lapply(str_split(old_rownames, " "), 
                                      FUN = str_sub, 1, 3), paste0, collapse = "_"))
new_rownames[new_rownames == "Myo_sco"] <- c("Myo_sco_1", "Myo_sco_2")
rownames(atlantic.fish.traits) <- new_rownames
 
norm <- decostand(atlantic.fish.traits, "normalize")

#quick mutate to remove negative values in these columns#
atlantic.fish.traits <- atlantic.fish.traits %>%
  mutate(Fresh = abs(Fresh),
         Brack = abs(Brack),
         Saltwater = abs(Saltwater))

#Data cleaning of unnecessary columns, such as references#

fish.traits <- select(atlantic.fish.traits, -c ("Genus","Subfamily", "Family", "Order", "Class", "Phylum","Atlantic..Northwest","Common", "DemersPelagRef", "MigratRef", "DepthRangeRef", "DepthComRef", "VerticalMoveRef", "LatMin", "LatMax", "LatRef", "LongevityWildRef", "MaxLengthRef", "CommonLengthRef", "MainCatchingMethod", "II", "MSeines", "MGillnets", "MCastnets", "MTraps", "MTrawls", "MDredges", "MLiftnets", "MHooksLines", "MSpears", "MOther", "MobilityRef", "FecundityRef", "SpawnRef", "LarvalRef", "LLRef", "LVRef", "PLDRef", "ADRef", "RaftingRef", "SchoolingRef", "Notes", "Comments", "Website.Link"))

#Veiw current data format#
View(atlantic.fish.traits)

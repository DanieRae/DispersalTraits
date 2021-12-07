library("dplyr")
library("tidyr")

#try to load large data and view it#

atlantic.fish.traits <- read.csv("AtlanticFishTraits.csv", fileEncoding ="UTF-8-BOM", stringsAsFactors = TRUE)

#this was recommended to keep data in an easier format for analysis#
#atlantic.fish.traits <- as_tibble(atlantic.fish.traits)# this is messing with something so im removing it for now

#need to cut out all the dead space after the last species, not sure why that was there but bye bye#
atlantic.fish.traits <- atlantic.fish.traits %>% slice(1:117)

#Making species name the row names

rownames(atlantic.fish.traits) <- atlantic.fish.traits[,1]
atlantic.fish.traits[,1] <- NULL

atlantic.fish.traits <- atlantic.fish.traits %>%
  mutate(Fresh = abs(Fresh),
         Brack = abs(Brack),
         Saltwater = abs(Saltwater))

View(atlantic.fish.traits)

#data looks like it loaded properly, has all necessary NA's for missing info.#
#lets try to look at some of the structure#

str(atlantic.fish.traits)

hist(atlantic.fish.traits$LarvalLength, breaks = 20)

hist(atlantic.fish.traits$Length, breaks = 20)

hist(atlantic.fish.traits$AvePLD, breaks = 20)

count(atlantic.fish.traits$LarvalStrategy)

# Simple Bar Plot
counts <- table(atlantic.fish.traits$LarvalStrategy)
barplot(counts, main="Larval strategy",
        xlab="Stragety type")

plot(atlantic.fish.traits$Length ~ atlantic.fish.traits$LarvalLength)


#Removing unnecessary columns 

fish.traits <- select(atlantic.fish.traits, -c ("Genus","Subfamily", "Family", "Order", "Class", "Phylum","Atlantic..Northwest","Common", "DemersPelagRef", "MigratRef", "DepthRangeRef", "DepthComRef", "VerticalMoveRef", "LatMin", "LatMax", "LatRef", "LongevityWildRef", "MaxLengthRef", "CommonLengthRef", "MainCatchingMethod", "II", "MSeines", "MGillnets", "MCastnets", "MTraps", "MTrawls", "MDredges", "MLiftnets", "MHooksLines", "MSpears", "MOther", "MobilityRef", "FecundityRef", "SpawnRef", "LarvalRef", "LLRef", "LVRef", "PLDRef", "ADRef", "RaftingRef", "SchoolingRef", "Notes", "Comments", "Website.Link"))




#Shall we analysis?

head(fish.traits)


# Apply function FUN to each column in the data.frame, returns
# a vector containing how many NAs are in a given column
fish.traits.NAs <- apply(fish.traits, FUN=function(x){sum(is.na(x))}, MARGIN = 2)

# Identify which columns have 0 NAs
fish.traits.noNas.columns <- names(fish.traits.NAs[fish.traits.NAs == 0])

fish.traits.noNA <- fish.traits[, fish.traits.noNas.columns]


#FD package is suppose to be able to use many different data types
library("FD")
library("gawdis")

# Calculate the distance matrix using Gower coefficient (S15 in numerical ecology book
# page 297)
fish.traits.dist <- gowdis(fish.traits, asym.bin = NULL)
fish.traits.dist

fish.traits.dist.noNa <- gowdis(fish.traits.noNA, asym.bin = NULL)
fish.traits.dist.noNa

image(as.matrix(fish.traits.dist) - as.matrix(fish.traits.dist.noNa), asp=1)

fish.traits.noDispersal <- dplyr::select(fish.traits, -AverageDispersal)
fish.traits.dist.noDispersal <- gowdis(fish.traits.noDispersal, asym.bin = NULL)
fish.traits.dist.noDispersal

image(as.matrix(fish.traits.dist) - as.matrix(fish.traits.dist.noDispersal), asp=1)

fish.traits.60percent <- dplyr::select(fish.traits, -c ("AverageDispersal", "CommonLength", "LTypeComM", "CommonLengthF", "FecundityMin","FecundityMax", "FecundityMean", "FecundityType","SpawningSeason", "LarvalLength","MinPLD","MaxPLD","AvePLD","MinDispersal", "MaxDispersal", "Rafting","Schooling", "WeightFemale", "LTypeComF"))
fish.traits.dist.60percent <- gowdis(fish.traits.60percent, asym.bin = NULL)
fish.traits.dist.60percent

image(as.matrix(fish.traits.dist.noDispersal) - as.matrix(fish.traits.dist.60percent), asp=1)

fish.traits.90percent <- dplyr::select(fish.traits, -c ("AverageDispersal", "CommonLengthF", "FecundityMean", "FecundityType","MinPLD","MaxPLD","AvePLD" , "MinDispersal","MaxDispersal", "WeightFemale", "LTypeComF"))
fish.traits.dist.90percent <- gowdis(fish.traits.90percent, asym.bin = NULL)
fish.traits.dist.90percent

image(as.matrix(fish.traits.dist) - as.matrix(fish.traits.dist.90percent), asp=1)

# Attempting to cluster the distance metrics using UPGMA
fish.traits.UPGMA <- hclust(fish.traits.dist, method = "average")
plot(fish.traits.UPGMA)

fish.traits.UPGMA.noDispersal <- hclust(fish.traits.dist.noDispersal, method = "average")
plot(fish.traits.UPGMA.noDispersal)

fish.traits.UPGMA.60percent <- hclust(fish.traits.dist.60percent, method = "average")
plot(fish.traits.UPGMA.60percent)

fish.traits.UPGMA.90percent <- hclust(fish.traits.dist.90percent, method = "average")
plot(fish.traits.UPGMA.90percent)

heatmap(as.matrix(fish.traits.dist))

# Cut the tree at k groups
k <- 11
#fish.traits.10groups <- cutree(fish.traits.UPGMA.60percent, k)

#fish.traits.10groups <- cutree(fish.traits.UPGMA.noDispersal, k)
# Shameless code steal from numerical ecology in R

# Reorder clusters
fish.traits.noDispersal.UPGMA.reordered <- reorder.hclust(fish.traits.UPGMA.noDispersal, fish.traits.dist.noDispersal)

fish.traits.60percent.UPGMA.reordered <- reorder.hclust(fish.traits.UPGMA.60percent, fish.traits.dist.60percent)

fish.traits.90percent.UPGMA.reordered <- reorder.hclust(fish.traits.UPGMA.90percent, fish.traits.dist.90percent)
# Plot reordered dendrogram with group labels
plot(
  fish.traits.90percent.UPGMA.reordered,
  hang = -1,
   xlab = "11 groups",
  # sub = "",
  # ylab = "Height",
   main = "Atlantic fish dispersal clusters with 90% traits",
  labels = cutree(fish.traits.90percent.UPGMA.reordered, k = k)
)
rect.hclust(fish.traits.90percent.UPGMA.reordered, k = k)
hcoplot(fish.traits.90percent.UPGMA.reordered,fish.traits.dist.90percent, k = k)

plot(
  fish.traits.60percent.UPGMA.reordered,
  hang = -1,
  # xlab = "4 groups",
  # sub = "",
  # ylab = "Height",
  # main = "Chord - Ward (reordered)",
  labels = cutree(fish.traits.60percent.UPGMA.reordered, k = k)
)
rect.hclust(fish.traits.60percent.UPGMA.reordered, k = k)
# Plot the final dendrogram with group colors (RGBCMY...)
# Fast method using the additional hcoplot() function:
# !!! Sourcing the function first
source("functions/hcoplot.R")
hcoplot(fish.traits.UPGMA.90percent, fish.traits.dist.90percent, k = k)

kmeans(fish.traits.dist)


# test2 <- gowdis(fish.traits.noNA, asym.bin = NULL)
# test2
#
# heatmap(as.matrix(test2))

# ga.fish.traits.dist <- gawdis(fish.traits)

functional.diversity.FT <- dbFD(fish.traits.dist,
              corr = "lingoes")

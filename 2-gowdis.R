#DATA REQUIRED PAGE 1 - LINE 1 TO 111#
#THIS SHEET RUNS THE TRAIT FILTERING FOR REMOVINGS DATA WITH TOO MANY NA'S --- USED IN PAGE 3-CLUSTERING, PAGE 4-PCOA, AND PAGE 5-MAPS#

#Libraries---- 
library(FD)
library(skimr)

#DATA FILTERING----
# Apply function FUN to each column in the data.frame, returns
# a vector containing how many NAs are in a given column
fish.traits.NAs <- apply(fish.traits, FUN=function(x){sum(is.na(x))}, MARGIN = 2)

#40NA ----
# Identify which columns have less than 40 NAs
fish.traits.40Nas.columns <- names(fish.traits.NAs[fish.traits.NAs <= 40])

fish.traits.40NA <- fish.traits[, fish.traits.40Nas.columns]

##SKIM DATA STRUC----
skim.40 <- skim(fish.traits.40NA)

#60NA ----
# Identify which columns have less than 60 NAs
fish.traits.60NA.column<- names(fish.traits.NAs[fish.traits.NAs <= 58])

fish.traits.60NA <- fish.traits[, fish.traits.60NA.column]

skim.60 <- skim(fish.traits.60NA)
#GOWERS MATRIX ----
#FD package is suppose to be able to use many different data types
# Calculate the distance matrix using Gower coefficient (S15 in numerical ecology book
# page 297)
#
fish.traits.dist <- gowdis(fish.traits, asym.bin = NULL)

fish.traits.40NA.dist <- gowdis(fish.traits.40NA, asym.bin = NULL)

fish.traits.60NA.dist <- gowdis(fish.traits.60NA, asym.bin = NULL)

#ONEHOT----

fish.traits.40NA.hot <- one_hot(as.data.table(fish.traits.40NA))
fish.traits.40NA.hot.dist <-gowdis(fish.traits.40NA.hot, asym.bin = NULL)

fish.traits.60NA.hot <- one_hot(as.data.table(fish.traits.60NA))
fish.traits.60NA.hot.dist <-gowdis(fish.traits.60NA.hot, asym.bin = NULL)


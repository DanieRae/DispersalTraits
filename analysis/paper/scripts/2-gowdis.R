#DATA REQUIRED PAGE 1 - LINE 1 TO 186#
#THIS SHEET RUNS THE TRAIT FILTERING FOR REMOVINGS DATA WITH TOO MANY NA'S
# USED IN PAGE 3-CLUSTERING, PAGE 4-PCOA, AND PAGE 5-MAPS#

# Install and load packages ----

# install.packages("FD")
# install.packages("skimr")
# install.packages("mltools")
# install.packages("data.table")

library(FD)
library(skimr)
library(mltools)
library(data.table)

#Data filtering ----
# Apply function 'FUN' to each column in the data.frame, returns
# a vector containing how many NAs are in a given column

fish.traits.NAs <- apply(fish.traits, FUN = function(x){sum(is.na(x))}, MARGIN = 2)

#40NA ----
# Identify which columns have less than 40 NAs
fish.traits.40Nas.columns <- names(fish.traits.NAs[fish.traits.NAs <= 40])

fish.traits.40NA <- fish.traits[, fish.traits.40Nas.columns]

#GOWERS MATRIX ----
#FD package is suppose to be able to use many different data types
# Calculate the distance matrix using Gower coefficient (S15 in numerical ecology book
# page 297)
#

fish.traits.40NA.dist <- gowdis(fish.traits.40NA, asym.bin = NULL)

#ONEHOT---- might remove this section

fish.traits.40NA.hot <- one_hot(as.data.table(fish.traits.40NA))


#END ---

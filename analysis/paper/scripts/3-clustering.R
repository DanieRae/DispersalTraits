#DATA REQUIRED FROM PAGE 2#
#THIS SHEET RUNS THE HEIARCHICAL CLUSTERING OF THE FISH TRAITS FOR THE THREE DIFFERENT DATA FILTERING OPTIONS ---- USED IN PAGE 5-MAPS#

# Install and load packages ----

# install.packages("cluster")
# install.packages("gclus")

library(cluster)
library(gclus)
library(here)

#Cluster using UPGMA----
fish.traits.40NA.UPGMA <- hclust(fish.traits.40NA.dist, method = "average")
plot(fish.traits.40NA.UPGMA, hang = -1)

## Reorder clusters----
# fish.traits.UPGMA.reordered <- reorder(fish.traits.UPGMA, fish.traits.dist)
# plot(fish.traits.UPGMA.reordered, hang =-1, main = "Fish Dispersal Clusters", xlab = "Dispersal Trait Distances", sub = "Method:UPGMA")


fish.traits.40NA.UPGMA.reordered <- reorder(fish.traits.40NA.UPGMA,
                                            fish.traits.40NA.dist)

## Silhouette Width -----

#How to determine how many clusters? Silhouette widths shows the optimum number of clusters#

fish_40clust <- fish.traits.40NA.UPGMA.reordered

Si <- numeric(nrow(fish.traits.40NA))
for (k in 2:(nrow(fish.traits.40NA) - 1)) {
  sil <- silhouette(cutree(fish_40clust, k = k), fish.traits.40NA.dist)
  Si[k] <- summary(sil)$avg.width
}

k.best <- which.max(Si)
plot(
  1:nrow(fish.traits.40NA),
  Si,
  type = "h",
  main = "Silhouette-optimal number of clusters",
  xlab = "K(number of clusters)",
  ylab = "Average silhouette width"
)
axis(
  1,
  k.best,
  paste("optiman", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best,
       max(Si),
       pch = 16,
       col = "red",
       cex = 1.5)

#returns an optimum of k=37, however there is very little variation between k=25-45#

##PLOT - Dendogram ----

plot(
  fish.traits.40NA.UPGMA.reordered,
  hang = -1,
  xlab = NULL,
  sub = NULL,
  ylab = "Height",
  main = NULL)

dendogramPlot <- rect.hclust(fish.traits.40NA.UPGMA.reordered, k = 25)

# uncomment this path to save figure
#ggsave(path = "figures", "dendogramPlot.png", dendogramPlot, width =  10, height = 10)

# MIGHT REMOVE THIS SECTION
#Plot the final dendrogram with group colors (RGBCMY...)
# Fast method using the additional hcoplot() function:
# !!! Sourcing the function first
# source("R/hcoplot.R")
# hcoplot(fish.traits.40NA.UPGMA.reordered, fish.traits.40NA.dist, k = 25)


###ADD CLUSTER ID TO TRAIT DATA----

cutree <- cutree(fish.traits.40NA.UPGMA.reordered, k = 25)

fish.traits.40NA.clust <- cbind(fish.traits.40NA, clusterID = cutree)

fish.traits.40NA.clust$clusterID <- as.factor(fish.traits.40NA.clust$clusterID)

# uncomment to write to CSV, used to determine what traits define each group
# write.csv(fish.traits.40NA.clust,
#           here("analysis","data","derived_data", "FishTraitsClustered25GR.csv"),
#           row.names = TRUE)


###ADD CLUSTER ID TO ABUND DATA ----
##Lets add the cluster IDs to the abundance data.

clust.ID <- select(fish.traits.40NA.clust, c("clusterID"))
clust.ID$taxa_name <- row.names(clust.ID)

fish.abun.clust <-
  merge(x = fish.abun.clean,
        y = clust.ID,
        by = "taxa_name",
        all.x = TRUE)

#Redundant?  group by year and stratum see Page 5 line 78#
#
# fish.abun.clust.gr <- fish.abun.clust %>%
#   group_by(stratum, year_surv, clusterID) %>%
#   dplyr::summarise(clust_biomass = sum(final_biomass))
#
#   dplyr::summarise(Unique_FE = n_distinct(clusterID))

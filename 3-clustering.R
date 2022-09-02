#DATA REQUIRED FROM PAGE 2#
#THIS SHEET RUNS THE HEIARCHICAL CLUSTERING OF THE FISH TRAITS FOR THE THREE DIFFERENT DATA FILTERING OPTIONS ---- USED IN PAGE 5-MAPS#

#Libraries----
library(cluster)
library(gclus)

#Cluster using UPGMA----
fish.traits.UPGMA <- hclust(fish.traits.dist, method = "average")
plot(fish.traits.UPGMA,hang= -1)

fish.traits.40NA.UPGMA <- hclust(fish.traits.40NA.dist, method = "average")
plot(fish.traits.40NA.UPGMA, hang= -1)

fish.traits.60NA.UPGMA <- hclust(fish.traits.60NA.dist, method = "average")
plot(fish.traits.60NA.UPGMA, hang= -1)


#comparing heat map#

dend <- as.dendrogram(fish.traits.40NA.UPGMA)
heatmap(
  as.matrix(fish.traits.60NA.dist),
  Rowv = dend,
  symm = T,
  margin=c(3,3))

plot(dend)

# Shameless code steal from numerical ecology in R

## Reorder clusters----
fish.traits.UPGMA.reordered <- reorder(fish.traits.UPGMA, fish.traits.dist)
plot(fish.traits.UPGMA.reordered, hang =-1, main = "Fish Dispersal Clusters", xlab = "Dispersal Trait Distances", sub = "Method:UPGMA")


fish.traits.40NA.UPGMA.reordered <- reorder(fish.traits.40NA.UPGMA, 
                                            fish.traits.40NA.dist)
plot(fish.traits.40NA.UPGMA.reordered, hang =-1)


fish.traits.60NA.UPGMA.reordered <- reorder(fish.traits.60NA.UPGMA, 
                                            fish.traits.60NA.dist)
plot(fish.traits.60NA.UPGMA.reordered, hang =-1)


# SILHOUETTE WIDTH -------------------------------------------------------------------

#How to determine how many clusters? Silhouette widths, I don't want to cluster based on width anymore, id prefer to cluster based on height. Will keep this to show that it

#fish_clust <- fish.traits.UPGMA.reordered 
fish_40clust <- fish.traits.40NA.UPGMA.reordered
#fish_60clust <- fish.traits.60NA.UPGMA.reordered

Si <- numeric(nrow(fish.traits.40NA))
for (k in 2:(nrow(fish.traits.40NA)-1)) {
  sil <- silhouette(cutree(fish_40clust, k=k), fish.traits.40NA.dist)
  Si[k] <- summary(sil)$avg.width
}

k.best <- which.max(Si)
plot(
  1:nrow(fish.traits.40NA),
  Si,
  type="h",
  main="Silhouette-optimal number of clusters",
  xlab= "K(number of clusters)",
  ylab="Average silhouette width"
)
axis(
  1,
  k.best,
  paste("optiman", k.best, sep = "\n"),
  col= "red",
  font = 2,
  col.axis="red"
)
points(k.best,
       max(Si),
       pch=16,
       col="red",
       cex=1.5)


#DENDOGRAMS ----

# Plot reordered dendrogram with group labels
# plot(
#   fish.traits.60percent.UPGMC.reordered,
#   hang = -1,
#   xlab = "11 groups",
#   # sub = "",
#   # ylab = "Height",
#   main = "Atlantic fish dispersal clusters with 90% traits",
#   labels = cutree(fish.traits.60percent.UPGMC.reordered, k = k)
# )
# rect.hclust(fish.traits.60percent.UPGMC.reordered, k = k)
# hcoplot(fish.traits.60percent.UPGMC.reordered,fish.traits.dist.60percent, k = k)

plot(
  fish.traits.40NA.UPGMA.reordered,
  hang = -1,
  xlab = "Groups",
  sub = "(Functional Entities)",
  ylab = "Height",
  main = "Fish Species Clustered by Dispersal Traits")
#labels = cutree(fish.traits.UPGMA.60percent, k =33)

rect.hclust(fish.traits.40NA.UPGMA.reordered, k=25)

# Plot the final dendrogram with group colors (RGBCMY...)
# Fast method using the additional hcoplot() function:
# !!! Sourcing the function first
source("functions/hcoplot.R")
hcoplot(fish.traits.40NA.UPGMA.reordered, fish.traits.40NA.dist, k = 25)

#DATA MANIPULATION ----
##ADD CLUSTER ID TO TRAIT DATA----

cutree <- cutree(fish.traits.40NA.UPGMA.reordered, k =25)

fish.traits.40NA.clust <- cbind(fish.traits.40NA, clusterID = cutree)

fish.traits.40NA.clust$clusterID <- as.factor(fish.traits.40NA.clust $clusterID)

#Write to CSV to determine what traits define each group
write.csv(fish.traits.40NA.clust,"C:\\Users\\Danielle\\Documents\\Graduate school\\Master Thesis\\FishTraitsClustered25GR.csv", row.names = TRUE)

###CLUST W 5 groups ----
cutree5 <- cutree(fish.traits.40NA.UPGMA.reordered, k =8)

fish.traits.40NA.clust5 <- cbind(fish.traits.40NA, clusterID = cutree5)

fish.traits.40NA.clust5$clusterID <- as.factor(fish.traits.40NA.clust5 $clusterID)

#Write to CSV to determine what traits define each group
write.csv(fish.traits.40NA.clust5,"C:\\Users\\Danielle\\Documents\\Graduate school\\Master Thesis\\FishTraitsClustered8GR.csv", row.names = TRUE)


cutree <- cutree(fish.traits.60NA.UPGMA.reordered, k =25)

fish.traits.60NA.clust <- cbind(fish.traits.60NA, clusterID = cutree)

fish.traits.60NA.clust $clusterID <- as.factor(fish.traits.60NA.clust $clusterID)

##ADD CLUSTER ID TO ABUND DATA ----
##Lets add the cluster IDs to the abundance data. Seems like it could be useful

clust.ID <- select (fish.traits.40NA.clust, c ("clusterID"))
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

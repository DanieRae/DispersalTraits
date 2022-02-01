#Shall we analysis - removing data with too many NA for use, distance matrix and clustering#

head(fish.traits)

# Apply function FUN to each column in the data.frame, returns
# a vector containing how many NAs are in a given column
fish.traits.NAs <- apply(fish.traits, FUN=function(x){sum(is.na(x))}, MARGIN = 2)

# Identify which columns have 0 NAs
fish.traits.noNas.columns <- names(fish.traits.NAs[fish.traits.NAs == 0])

fish.traits.noNA <- fish.traits[, fish.traits.noNas.columns]


#FD package is suppose to be able to use many different data types
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

#keeping traits with up to 60% missing data#
fish.traits.60percent <- dplyr::select(fish.traits, -c ("AverageDispersal", "CommonLength", "LTypeComM", "CommonLengthF","LTypeMaxM", "MaxWeightRef", "FecundityMin","FecundityMax", "FecundityMean", "FecundityType","SpawningSeason", "LarvalLength","MinPLD","MaxPLD","AvePLD","MinDispersal", "MaxDispersal", "Rafting","Schooling", "WeightFemale", "LTypeComF"))

fish.traits.dist.60percent <- gowdis(fish.traits.60percent, asym.bin = NULL)
fish.traits.dist.60percent

image(as.matrix(fish.traits.dist) - as.matrix(fish.traits.dist.60percent), asp=1)


# Attempting to cluster the distance metrics using UPGMA
fish.traits.UPGMA <- hclust(fish.traits.dist, method = "average")
plot(fish.traits.UPGMA)

fish.traits.UPGMA.noDispersal <- hclust(fish.traits.dist.noDispersal, method = "average")
plot(fish.traits.UPGMA.noDispersal)

fish.traits.UPGMA.60percent <- hclust(fish.traits.dist.60percent, method = "average")
plot(fish.traits.UPGMA.60percent, hang= -1)

#Lets compare the clustering method to centroid, what will work better
fish.traits.UPGMC.60percent <- hclust(fish.traits.dist.60percent, method = "centroid")
plot(fish.traits.UPGMC.60percent)

#comparing heat map#

dend <- as.dendrogram(fish.traits.UPGMA.60percent)
heatmap(
  as.matrix(fish.traits.dist.60percent),
  Rowv = dend,
  symm = T,
  margin=c(3,3))



# Shameless code steal from numerical ecology in R

# Reorder clusters doesn't seem important -------------------------------

fish.traits.noDispersal.UPGMA.reordered <- reorder(fish.traits.UPGMA.noDispersal, 
                                                   fish.traits.dist.noDispersal)


fish.traits.60percent.UPGMA.reordered <- reorder(fish.traits.UPGMA.60percent, 
                                                 fish.traits.dist.60percent)

plot(fish.traits.60percent.UPGMA.reordered, hang =-1)


# Plots -------------------------------------------------------------------

#How to determine how many clusters? Silhouette widths#

fish_60clust <- fish.traits.UPGMA.60percent

Si <- numeric(nrow(fish.traits.60percent))
for (k in 2:(nrow(fish.traits.60percent)-1)) {
  sil <- silhouette(cutree(fish_60clust, k=k), fish.traits.dist.60percent)
  Si[k] <- summary(sil)$avg.width
}

k.best <- which.max(Si)
plot(
  1:nrow(fish.traits.60percent),
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



k <- 33

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
  fish.traits.UPGMA.60percent,
  hang = -1,
  xlab = "33 groups",
  sub = "(Dispersal trait clusters)",
  ylab = "Height",
  main = "Fish Species (reordered)")
,
  #labels = cutree(fish.traits.UPGMA.60percent, k =33)
)
rect.hclust(fish.traits.UPGMA.60percent, k = 33)
# Plot the final dendrogram with group colors (RGBCMY...)
# Fast method using the additional hcoplot() function:
# !!! Sourcing the function first
source("functions/hcoplot.R")
hcoplot(fish.traits.60percent.UPGMA.reordered, fish.traits.dist.60percent, k = 33)

#adding cluster ID to data

cutree <- cutree(fish.traits.UPGMA.60percent, k =33)

clust.fish <- cbind(fish.traits.60percent, clusterID = cutree)

clust.fish$clusterID <- as.factor(clust.fish$clusterID)


##ONEHOT

fish.traits.60percent.hot <- one_hot(as.data.table(fish.traits.60percent))

fish.traits.60percent.hot.dist <-gowdis(fish.traits.60percent.hot, asym.bin = NULL)

# PCOA two function produce similar plots, cmdscale has open circles and pcoa uses the species names 

# -------------------------------------------------------------------------

## CMDSCALE pcoa method with ggplot 

fish.cmd <- cmdscale(fish.traits.dist.60percent, eig = TRUE, x.ret = TRUE, k =2)

fish.cdm.eig <- round (fish.cmd$eig/sum(fish.cmd$eig)*100,1)

fish.cmd.values <- fish.cmd$points

fish.cmd.data <- data.frame(Sample = rownames(fish.cmd.values), 
                            X= fish.cmd.values [,1], 
                            Y= fish.cmd.values [,2])
fish.cmd.data$clusterID <- clust.fish$clusterID

plot1 <- ggplot(data =fish.cmd.data, aes (x=X, y=Y, 
                                          color = clusterID)
                ) +
  geom_point(size =3)+
  #geom_text()+
  theme_bw()+
  xlab(paste("PCoA 1 -",fish.cdm.eig[1], "%", sep ="" ))+
  ylab(paste("PCoA 2 -",fish.cdm.eig[2], "%", sep ="" ))+
  ggtitle("PCOA") +
  geom_path() + 
  theme(legend.position="none", asp=1) +
  scale_color_grey()

# plot1

# Code taken from ape::biplot.pcoa
get_U_matrix <- function(Y, Eig, vectors, plot.axes=c(1,2), 
                         dir.axis1=1, dir.axis2=1){
  # browser()
  diag.dir <- diag(c(dir.axis1,dir.axis2))
  vectors[,plot.axes] <- vectors[,plot.axes] %*% diag.dir
  n <- nrow(Y)
  points.stand <- scale(vectors[,plot.axes])
  S <- cov(Y, points.stand)
  U <- S %*% diag((Eig[plot.axes]/(n-1))^(-0.5))
  colnames(U) <- colnames(vectors[,plot.axes])
  return(U)
}

U <-  tidyr::drop_na(as.data.frame(
  get_U_matrix(scale(fish.traits.60percent.hot), fish.cdm.eig, fish.cmd.values)))/3 
colnames(U) <- c("x", "y")
U$var_name <- rownames(U)

#current row names do not add up to the row names from the 60percent.hot, why is this happening?#

old_var_name <- U$var_name
new_var_name <- c("eel-like", "elongated", "normal", "other","short","fresh", "brack","bathydemersal", "bathypelagic","benthopelagic", "demersal","neritic","oceanic", "shallow", "vertical", "non-vertical", "length")

U$var_name <- new_var_name

plot1 + geom_segment(data = U,
                     aes(x = 0, xend = x, 
                         y = 0, yend = y),
                     arrow = arrow(length = unit(0.25, "cm")), colour = "red") +
  geom_text(data = U,
             aes(x = x*1.1, y = y*1.1, 
                 label = stringr::str_wrap(var_name,20)), colour = "red")
  
# -------------------------------------------------------------------------

# Biplot with pcoa method using lingoes correction
fish.pcoa <- pcoa(fish.traits.dist.60percent, correction = "lingoes")

numeric_columns_id <- attributes(fish.traits.dist.60percent)$Types == "C"
fish.traits.60percent.numeric <- fish.traits.60percent[,numeric_columns_id]

# biplot.pcoa(fish.pcoa, scale(fish.traits.60percent))
# abline(h=0, lty=3)
# abline(v=0, lty=3)

hot.dist.pcoa <- pcoa(fish.traits.60percent.hot.dist, correction = "lingoes")

biplot.pcoa(hot.dist.pcoa,scale(fish.traits.60percent.hot))
abline(h=0, lty=3)
abline(v=0, lty=3)

test <- ggplot(data.frame(clust.fish),aes(x=clusterID)) 

test + geom_histogram(fill="lightgreen")

test + geom_histogram(fill="blue", aes(y=..count../sum(..count..)))
test

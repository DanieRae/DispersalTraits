#Shall we analysis - removing data with too many NA for use, distance matrix and clustering#

head(fish.traits)

# Apply function FUN to each column in the data.frame, returns
# a vector containing how many NAs are in a given column
fish.traits.NAs <- apply(fish.traits, FUN=function(x){sum(is.na(x))}, MARGIN = 2)

# Identify which columns have less than 40 NAs
fish.traits.40Nas.columns <- names(fish.traits.NAs[fish.traits.NAs <= 40])

fish.traits.40NA <- fish.traits[, fish.traits.40Nas.columns]


#FD package is suppose to be able to use many different data types
# Calculate the distance matrix using Gower coefficient (S15 in numerical ecology book
# page 297)
#
fish.traits.dist <- gowdis(fish.traits, asym.bin = NULL)
fish.traits.dist

fish.traits.40NA.dist <- gowdis(fish.traits.40NA, asym.bin = NULL)
fish.traits.40NA.dist

image(as.matrix(fish.traits.dist) - as.matrix(fish.traits.40NA.dist), asp=1)


#keeping traits with up to 60% missing data#
fish.traits.60NA.column<- names(fish.traits.NAs[fish.traits.NAs <= 58])

fish.traits.60NA <- fish.traits[, fish.traits.60NA.column]

fish.traits.60NA.dist <- gowdis(fish.traits.60NA, asym.bin = NULL)
fish.traits.60NA.dist

image(as.matrix(fish.traits.40NA.dist) - as.matrix(fish.traits.60NA.dist), asp=1)


# Attempting to cluster the distance metrics using UPGMA
fish.traits.UPGMA <- hclust(fish.traits.dist, method = "average")
plot(fish.traits.UPGMA,hang= -1)

fish.traits.40NA.UPGMA <- hclust(fish.traits.40NA.dist, method = "average")
plot(fish.traits.40NA.UPGMA, hang= -1)

fish.traits.60NA.UPGMA <- hclust(fish.traits.60NA.dist, method = "average")
plot(fish.traits.60NA.UPGMA, hang= -1)


#comparing heat map#

dend <- as.dendrogram(fish.traits.60NA.UPGMA)
heatmap(
  as.matrix(fish.traits.60NA.dist),
  Rowv = dend,
  symm = T,
  margin=c(3,3))



# Shameless code steal from numerical ecology in R

# Reorder clusters -------------------------------
fish.traits.UPGMA.reordered <- reorder(fish.traits.UPGMA, fish.traits.dist)
plot(fish.traits.UPGMA.reordered, hang =-1)


fish.traits.40NA.UPGMA.reordered <- reorder(fish.traits.40NA.UPGMA, 
                                            fish.traits.40NA.dist)
plot(fish.traits.40NA.UPGMA.reordered, hang =-1)


fish.traits.60NA.UPGMA.reordered <- reorder(fish.traits.60NA.UPGMA, 
                                            fish.traits.60NA.dist)
plot(fish.traits.60NA.UPGMA.reordered, hang =-1)


# Plots -------------------------------------------------------------------

#How to determine how many clusters? Silhouette widths, I don't want to cluster based on width anymore, id prefer to cluster based on height.

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

rect.hclust(fish.traits.40NA.UPGMA.reordered, k = 30)
# Plot the final dendrogram with group colors (RGBCMY...)
# Fast method using the additional hcoplot() function:
# !!! Sourcing the function first
source("functions/hcoplot.R")
hcoplot(fish.traits.40NA.UPGMA.reordered, fish.traits.40NA.dist, k = 30)

#adding cluster ID to data

cutree <- cutree(fish.traits.40NA.UPGMA.reordered, k =30)

clust.fish <- cbind(fish.traits.40NA, clusterID = cutree)

clust.fish$clusterID <- as.factor(clust.fish$clusterID)

write.csv(clust.fish,"C:\\Users\\Danielle\\Documents\\Grad school\\FishTraitsClustered.csv", row.names = TRUE)

##ONEHOT

fish.traits.40NA.hot <- one_hot(as.data.table(fish.traits.40NA))

fish.traits.40NA.hot.dist <-gowdis(fish.traits.40NA.hot, asym.bin = NULL)

# PCOA two function produce similar plots, cmdscale has open circles and pcoa uses the species names 

# -------------------------------------------------------------------------

## CMDSCALE pcoa method with ggplot 

fish.cmd <- cmdscale(fish.traits.40NA.dist, eig = TRUE, x.ret = TRUE, k =2)

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
  get_U_matrix(scale(fish.traits.40NA.hot), fish.cdm.eig, fish.cmd.values)))/5 
colnames(U) <- c("x", "y")
U$var_name <- rownames(U)

#current row names do not add up to the row names from the 60percent.hot, why is this happening?#

old_var_name <- U$var_name
new_var_name <- c("Eel-like", "Elongated", "Normal", "Other","Short","Fresh", "Brack","Bathydemersal", "Bathypelagic","Benthopelagic", "Demersal","Neritic","Oceanic", "Shallow", "Vertical N", "Vertical Y", "Length", "Larval demersal", "Larval direct", "Larval pelagic")

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

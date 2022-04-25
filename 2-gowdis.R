#Shall we analysis - removing data with too many NA for use, distance matrix and clustering#

# Apply function FUN to each column in the data.frame, returns
# a vector containing how many NAs are in a given column
fish.traits.NAs <- apply(fish.traits, FUN=function(x){sum(is.na(x))}, MARGIN = 2)

# Identify which columns have less than 40 NAs
fish.traits.40Nas.columns <- names(fish.traits.NAs[fish.traits.NAs <= 40])

fish.traits.40NA <- fish.traits[, fish.traits.40Nas.columns]


# Identify which columns have less than 40 NAs
fish.traits.60NA.column<- names(fish.traits.NAs[fish.traits.NAs <= 58])

fish.traits.60NA <- fish.traits[, fish.traits.60NA.column]

#FD package is suppose to be able to use many different data types
# Calculate the distance matrix using Gower coefficient (S15 in numerical ecology book
# page 297)
#
fish.traits.dist <- gowdis(fish.traits, asym.bin = NULL)
fish.traits.dist

fish.traits.40NA.dist <- gowdis(fish.traits.40NA, asym.bin = NULL)
fish.traits.40NA.dist

fish.traits.60NA.dist <- gowdis(fish.traits.60NA, asym.bin = NULL)
fish.traits.60NA.dist


##ONEHOT

fish.traits.40NA.hot <- one_hot(as.data.table(fish.traits.40NA))
fish.traits.40NA.hot.dist <-gowdis(fish.traits.40NA.hot, asym.bin = NULL)

fish.traits.60NA.hot <- one_hot(as.data.table(fish.traits.60NA))
fish.traits.60NA.hot.dist <-gowdis(fish.traits.60NA.hot, asym.bin = NULL)

# PCOA two function produce similar plots, cmdscale has open circles and pcoa uses the species names 

# CMDSCALE pcoa method with ggplot (RUN SHEET 5)---------------------------------------


fish.cmd <- cmdscale(fish.traits.40NA.dist, eig = TRUE, x.ret = TRUE, k =2)

fish.cdm.eig <- round (fish.cmd$eig/sum(fish.cmd$eig)*100,1)

fish.cmd.values <- fish.cmd$points

fish.cmd.data <- data.frame(Sample = rownames(fish.cmd.values), 
                            X= fish.cmd.values [,1], 
                            Y= fish.cmd.values [,2])
fish.cmd.data$clusterID <- clust.fish$clusterID

plot2 <- ggplot(data =fish.cmd.data, aes (x=X, y=Y, color = clusterID)
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

plot2

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


#old_var_name <- U$var_name
#new_var_name <- c("Eel-like", "Elongated", "Normal", "Other","Short","Fresh", "Brack","Bathydemersal", "Bathypelagic","Benthopelagic", "Demersal","Neritic","Oceanic", "Shallow", "Vertical N", "Vertical Y", "Length", "Larval demersal", "Larval direct", "Larval pelagic")

#U$var_name <- new_var_name

plot2 + geom_segment(data = U,
                     aes(x = 0, xend = x, 
                         y = 0, yend = y),
                     arrow = arrow(length = unit(0.25, "cm")), colour = "red") +
  geom_text(data = U,
             aes(x = x*1.1, y = y*1.1, 
                 label = stringr::str_wrap(var_name,20)), colour = "red")

#60na------
fish.cmd <- cmdscale(fish.traits.60NA.dist, eig = TRUE, x.ret = TRUE, k =2)

fish.cdm.eig <- round (fish.cmd$eig/sum(fish.cmd$eig)*100,1)

fish.cmd.values <- fish.cmd$points

fish.cmd.data <- data.frame(Sample = rownames(fish.cmd.values), 
                            X= fish.cmd.values [,1], 
                            Y= fish.cmd.values [,2])
fish.cmd.data$clusterID <- clust.fish$clusterID

plot2 <- ggplot(data =fish.cmd.data, aes (x=X, y=Y)
) +
  geom_point(size =3)+
  #geom_text()+
  theme_bw()+
  xlab(paste("PCoA 1 -",fish.cdm.eig[1], "%", sep ="" ))+
  ylab(paste("PCoA 2 -",fish.cdm.eig[2], "%", sep ="" ))+
  ggtitle("PCOA 60") +
  #geom_path() + 
  theme(legend.position="none", asp=1) +
  scale_color_grey()

plot2

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
  get_U_matrix(scale(fish.traits.60NA.hot), fish.cdm.eig, fish.cmd.values)))/5 
colnames(U) <- c("x", "y")
U$var_name <- rownames(U)

plot2 + geom_segment(data = U,
                     aes(x = 0, xend = x, 
                         y = 0, yend = y),
                     arrow = arrow(length = unit(0.25, "cm")), colour = "red") +
  geom_text(data = U,
            aes(x = x*1.1, y = y*1.1, 
                label = stringr::str_wrap(var_name,20)), colour = "red")

# Biplot with pcoa method using lingoes correction ----------------------------------------


fish.pcoa <- pcoa(fish.traits.40NA.dist, correction = "lingoes")

numeric_columns_id <- attributes(fish.traits.40NA.dist)$Types == "C"
fish.traits.60percent.numeric <- fish.traits.60percent[,numeric_columns_id]

# biplot.pcoa(fish.pcoa, scale(fish.traits.60percent))
# abline(h=0, lty=3)
# abline(v=0, lty=3)

hot.dist.pcoa <- pcoa(fish.traits.40NA.hot.dist, correction = "lingoes")

biplot.pcoa(hot.dist.pcoa,scale(fish.traits.40NA.hot))
abline(h=0, lty=3)
abline(v=0, lty=3)



# PCOA two function produce similar plots, cmdscale has open circles and pcoa uses the species names 

# CMDSCALE pcoa method with ggplot####


fish.cmd <- cmdscale(fish.traits.40NA.dist, eig = TRUE, x.ret = TRUE, k =2)

fish.cdm.eig <- round (fish.cmd$eig/sum(fish.cmd$eig)*100,1)

fish.cmd.values <- fish.cmd$points

fish.cmd.data <- data.frame(Sample = rownames(fish.cmd.values), 
                            X= fish.cmd.values [,1], 
                            Y= fish.cmd.values [,2])
fish.cmd.data$clusterID <- fish.traits.40NA.clust$clusterID

#PCoA plot 1 ----
plot1 <- ggplot(data =fish.cmd.data, aes (x=X, y=Y, color = clusterID)
) +
  geom_point(size =3)+
  stat_ellipse()+
  #geom_text()+
  theme_bw()+
  xlab(paste("PCoA 1: ",fish.cdm.eig[1], "%", sep ="" ))+
  ylab(paste("PCoA 2: ",fish.cdm.eig[2], "%", sep ="" ))+
  ggtitle("PCoA Marine fish dispersal") +
  theme(legend.position="none", asp=1)


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
  get_U_matrix(scale(fish.traits.40NA.hot), fish.cdm.eig, fish.cmd.values)))/7 
colnames(U) <- c("x", "y")
U$var_name <- rownames(U)


old_var_name <- U$var_name
new_var_name <- c("Eel-like", "Elongated", "Normal", "Other","Short","Fresh", "Brackish","Bathydemersal", "Bathypelagic","Benthopelagic", "Demersal","Neritic","Oceanic","CShallow","Vertical N", "Vertical Y", "Length", "Larval demersal", "Larval direct", "Larval pelagic")

U$var_name <- new_var_name

plot1 + geom_segment(data = U,
                     aes(x = 0, xend = x, 
                         y = 0, yend = y),
                     arrow = arrow(length = unit(0.25, "cm")), colour = "black") +
  geom_text(data = U,
            aes(x = x*1.1, y = y*1.1, 
                label = stringr::str_wrap(var_name,20)), colour = "black")

#PCoA plot2----
fish.cmd <- cmdscale(fish.traits.60NA.dist, eig = TRUE, x.ret = TRUE, k =2)

fish.cdm.eig <- round (fish.cmd$eig/sum(fish.cmd$eig)*100,1)

fish.cmd.values <- fish.cmd$points

fish.cmd.data <- data.frame(Sample = rownames(fish.cmd.values), 
                            X= fish.cmd.values [,1], 
                            Y= fish.cmd.values [,2])
fish.cmd.data$clusterID <- fish.traits.60NA.clust$clusterID

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



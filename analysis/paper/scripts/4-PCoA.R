# DATA REQUIRED FROM PAGE 2
# PRODUCES THE PCOA PLOTS FOR TRAIT SPACE VISUALIZATION
# NOT USED FOR ANOTHER PAGE
# NOTE PCOA two function produce similar plots, cmdscale has open circles and
# pcoa uses the species names

# Install and load packages ----
# install.packages("ggplot2")

library(ggplot2)
library(here)

# CMDSCALE - PCoA ----
# multidimensional scaling, to return a set of points that reflect the
# distances/dissimilarities between the trait values#
# returns a matrix of coordinates-points for dissimilarities and
# eigenvalues computed from scaling
fish.cmd <-
  cmdscale(fish.traits.40NA.dist,
           eig = TRUE,
           x.ret = TRUE,
           k = 2)


fish.cdm.eig <- round(fish.cmd$eig / sum(fish.cmd$eig) * 100, 1)

fish.cmd.values <- fish.cmd$points

fish.cmd.data <- data.frame(Sample = rownames(fish.cmd.values),
                            X = fish.cmd.values[, 1],
                            Y = fish.cmd.values[, 2])
fish.cmd.data$clusterID <- fish.traits.40NA.clust$clusterID

fish.cmd.data <-
  fish.cmd.data |>
  group_by(clusterID) |>
  mutate(cluster_size = n()) |>
  ungroup()

#PLOT - PCoA W/ Ellipses ----
# missing vectors #
PCOA.plot <-
  ggplot(data = fish.cmd.data, aes(x = X, y = Y)) +
  geom_point(size = 2,
             data = filter(fish.cmd.data, cluster_size < 4),
             color = "grey") +
  geom_point(size = 2,
             data = filter(fish.cmd.data, cluster_size >= 4),
             aes(color = clusterID)) +
  stat_ellipse(data = filter(fish.cmd.data, cluster_size >= 4),
              aes(color = clusterID)) +
  theme_bw() +
  xlab(paste("PCoA 1: ", fish.cdm.eig[1], "%", sep = "")) +
  ylab(paste("PCoA 2: ", fish.cdm.eig[2], "%", sep = "")) +
  theme(legend.position = "none", asp = 1)

# Code taken from ape::biplot.pcoa
# getting data for vector names and arrows
get_U_matrix <- function(Y,
                         Eig,
                         vectors,
                         plot.axes = c(1, 2),
                         dir.axis1 = 1,
                         dir.axis2 = 1) {
  # browser()
  diag.dir <- diag(c(dir.axis1, dir.axis2))
  vectors[, plot.axes] <- vectors[, plot.axes] %*% diag.dir
  n <- nrow(Y)
  points.stand <- scale(vectors[, plot.axes])
  S <- cov(Y, points.stand)
  U <- S %*% diag((Eig[plot.axes] / (n - 1)) ^ (-0.5))
  colnames(U) <- colnames(vectors[, plot.axes])
  return(U)
}

U <-  tidyr::drop_na(as.data.frame(get_U_matrix(
  scale(fish.traits.40NA.hot), fish.cdm.eig, fish.cmd.values
))) / 7
colnames(U) <- c("x", "y")
U$var_name <- rownames(U)


old_var_name <- U$var_name
new_var_name <-
  c(
    "Eel-like",
    "Elongated",
    "Normal",
    "Other",
    "Short",
    "Fresh",
    "Brackish",
    "Bathydemersal",
    "Bathypelagic",
    "Benthopelagic",
    "Demersal",
    "Neritic",
    "Oceanic",
    "Shallow",
    "Deep",
    "Length",
    "Larval demersal",
    "Larval direct",
    "Larval pelagic"
  )

U$var_name <- new_var_name

#add vectors to the base plot
PCOAplot <- PCOA.plot + geom_segment(
  data = U,
  aes(x = 0,
      xend = x,
      y = 0,
      yend = y),
  arrow = arrow(length = unit(0.25, "cm")),
  colour = "red") +
  geom_text(
    data = U,
    aes(x = x * 1.1,
        y = y * 1.1,
        label = stringr::str_wrap(var_name, 20)),
    colour = "black"
  )


# uncomment to save PCOAplot
# ggsave(path = here("analysis","figures"), "PCOAplot.png",
#        PCOAplot, width =  10, height = 10)


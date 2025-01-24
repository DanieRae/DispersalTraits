##Might delete this page. seems irrelevant.
library(ggplot2)

ggplot(data = fish.traits.40NA, aes(x = BodyShape)) +
  geom_bar()

ggplot(data = fish.traits.40NA, aes(x = DemersPelag)) +
  geom_bar()

ggplot(data = fish.traits.40NA, aes(x = LarvalStrategy)) +
  geom_bar()

ggplot(data = fish.traits.40NA, aes(x = DepthRangeDeep)) +
  geom_histogram()

ggplot(data = fish.traits.40NA, aes(x = DepthRangeShallow)) +
  geom_histogram()

ggplot(data = fish.traits.40NA, aes(x = Length)) +
  geom_histogram()

#Trying to make a stacked bar plot#

fish.traits.40NA.rownames <- fish.traits.40NA
taxa_name <- rownames(fish.traits.40NA)

#fish.traits.40NA.rownames <- cbind(taxa_name,fish.traits.40NA.rownames)
num <- unlist(lapply(fish.traits.40NA.rownames, is.numeric))
cat <- unlist(lapply(fish.traits.40NA.rownames, is.factor))

fish.traits.40NA.num <- fish.traits.40NA.rownames[,num]
fish.traits.40NA.num <- cbind(taxa_name,fish.traits.40NA.num)

fish.traits.40NA.cat <- fish.traits.40NA.rownames[,cat]
fish.traits.40NA.cat <- cbind(taxa_name,fish.traits.40NA.cat)

fish.traits.40.num.long <- fish.traits.40NA.num %>%
  pivot_longer(!taxa_name, names_to = "Trait", values_to = "count")

fish.traits.40.cat.long <- fish.traits.40NA.cat %>%
  pivot_longer(!taxa_name, names_to = "Trait", values_to = "count")

# Stacked + percent
ggplot(fish.traits.40.cat.long, aes(fill = taxa_name, y = count)) +
  geom_bar(position = "fill", stat = "count") +
  theme(legend.position = "none")



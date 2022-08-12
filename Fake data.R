library(ggplot2)
library(ggformula)
library(mosaicData)


#creating dataframe
x <- seq(0, 1, length=100)
m = 0.3
b = 0.483

df1 <- data.frame(x = x, y = dnorm(x, mean = 0.5, sd =0.5))
df2 <- data.frame (x = x, y=b+m*x)


p1 <- ggplot() +
  geom_path(data = df1, aes(x = x, y = y, color = "Prediction 1", size = 0.2 )) +
  geom_path(data = df2, aes(x = x, y = y, color = "Prediction 2", size = 0.2)) +
  ylab("Stability") +
  xlab("Dispersal Diversity ")+ 
  scale_y_continuous(breaks = NULL)+ 
  scale_x_continuous(breaks = NULL)+
  labs(color = "Predictions") 
  
p1


#Create a sequence of 100 equally spaced numbers between -4 and 4

#create a vector of values that shows the height of the probability distribution
#for each value in x



p1 <- gf_line(dnorm(x, mean = 0.5) ~ x, color = "orange", size = 1.5, xlab = "Dispersal diversity", ylab = "Stabilty") 

p2 <- gf_line((m*x+b) ~ x, color = "green", size = 1.5, xlab = "Dispersal diversity", ylab = "Stabilty") 

p3 <- p1+p2



gf_line()
gf_point(age ~ sex, alpha = 0.25, data = mosaicData::HELPrct)
gf_point(births ~ date, color = ~wday, data = mosaicData::Births78)
# lines make the exceptions stand out more prominently
gf_line(births ~ date, color = ~wday, data = mosaicData::Births78)
gf_path()
if (require(dplyr)) {
  data.frame(t = seq(1, 10 * pi, length.out = 400)) %>%
    mutate(x = t * cos(t), y = t * sin(t)) %>%
    gf_path(y ~ x, color = ~t)
}

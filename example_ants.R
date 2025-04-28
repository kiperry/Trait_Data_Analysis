###################################################################################
#
# Functional trait data analysis workshop
#
# Community Trait Means
# Functional Diversity
# PLS
# RLQ/Fourth Corner
#
# KI Perry; 28 April 2025
#
###################################################################################

# For this workshop, we will work with ant species collected in urban greenspaces
# Ants were collected in 2015 and 2016. We will work with both years pooled together
# Because there are biases when collecting ants via pitfall traps, we will use presence/absence data
# However, abundance data can also be used for these analyses

## Species data
a <- read.csv("Data/ant_assemblages_pooled.csv", row.names=1)
str(a)

## Trait data
t <- read.csv("Data/ant_traits_pooled.csv", row.names=1)
str(t)

# are any traits highly correlated?
library(ggplot2)
library(GGally)
cp <- ggpairs(t, upper = list(continuous = wrap("cor", size = 3, color = "black")))
cp + theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 10))

plot(t)
cor(t, method = c("pearson"), use = "complete.obs")

# will need to make a decision about removing any traits that are highly correlated or
# accounting for them in the indices

## Environmental data
env <- read.csv("Data/env_local_landscape.csv", row.names=1)
str(env)

# change categorical predictors to factors
env$Neighborhood <- as.factor(env$Neighborhood)
env$trmt <- as.factor(env$trmt)
str(env)







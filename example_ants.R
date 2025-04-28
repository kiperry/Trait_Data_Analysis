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

# check traits for normality
hist(t$wl)
hist(t$rhw)
hist(t$rml)
hist(log(t$rml + 1))
hist(t$rpw)
hist(t$rew)
hist(t$rsl)
hist(t$rfl)
hist(log(t$rfl + 1 ))
hist(t$rcl)

## Environmental data
env <- read.csv("Data/env_local_landscape.csv", row.names=1)
str(env)

# change categorical predictors to factors
env$Neighborhood <- as.factor(env$Neighborhood)
env$trmt <- as.factor(env$trmt)
str(env)

# Double check all species names are in the same order in both data sets
rownames(t) == colnames(a) 

library(FD)
library(picante)
library(gawdis)
library(betapart)
library(ade4)

#####################################################################################
# Community (weighted) means

#observed functional CWM
cwm.obs <- functcomp(t, as.matrix(a), CWM.type = "all")
cwm.obs

# CWM values for each site can be used in other analyses to test for treatment effects
boxplot(cwm.obs$wl ~ env$trmt)
boxplot(cwm.obs$rml ~ env$trmt)

#####################################################################################
# Functional diversity indices
# If dataset has a categorical variable, then must calculate an appropriate dissimilarity matrix

## Gower dissimilarity matrix
td.go <- as.matrix(gowdis(t))

# Indices will only calculate value if community has at least 3 species
# Communities that have fewer than 3 species, output will be NA

#Functional richness
fric <- dbFD(td.go[colnames(a), ], a, w.abun = F, stand.x = T)$FRic
fric

#Functional evenness
feve <- dbFD(td.go[colnames(a), ], a, w.abun = F, stand.x = T)$FEve
feve

#Functional divergence
fdiv <- dbFD(td.go[colnames(a), ], a, w.abun = F, stand.x = T)$FDiv
fdiv

#Functional dispersion
fdis <- dbFD(td.go[colnames(a), ], a, w.abun = F, stand.x = T)$FDis
fdis

fun.div <- as.data.frame(cbind(fric, feve, fdiv, fdis))
plot(fun.div)

# Functional diversity indices can be used in other analyses to test for treatment effects
boxplot(fun.div$fric ~ env$trmt)
boxplot(fun.div$feve ~ env$trmt)
boxplot(fun.div$fdiv ~ env$trmt)
boxplot(fun.div$fdis ~ env$trmt)

#####################################################################################
# Functional beta diversity

# Can use gower dissimilarity or calculate a gawdis dissimilarity matrix

## Gower dissimilarity matrix
td.go.2 <- gowdis(t)

## Gawdis dissimilarity matrix
td.ga <- gawdis(t, w.type = "optimized", opti.maxiter = 300,
                groups.weight = T, groups = c(1, 2, 2, 2, 2, 2, 3, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
attr(td.ga, "weights")

pco.ga <- dudi.pco(sqrt(td.ga), scannf = FALSE, nf = 5) # nf is number of axes
scatter(pco.ga)

pco.ga$li
sum(pco.ga$eig[1:5]) / sum(pco.ga$eig) 
sum(pco.ga$eig[1:3]) / sum(pco.ga$eig)
sum(pco.ga$eig[1:3]) / sum(pco.ga$eig)
sum(pco.ga$eig[1:2]) / sum(pco.ga$eig)

# check correlations among axes and traits
cor(pco.ga$li, t, use = "complete.obs")

# the number of axes is dependent on the lowest number of species collected across sites
# if we wanted to select 5 axes, the site with the lowest ant richness must equal 6
# our lowest richness site equals 5, so we have to select 4 axes
t.ax <- as.matrix(pco.ga$li[1:4])

# returns pairwise between-site values of each functional beta-diversity component
fun.beta <- functional.beta.pair(a, t.ax, index.family = "jaccard") # family = sorensen for abundance data
str(fun.beta)
# jtu = functional turnover
# jne = functional nestedness
# jac = total functional beta diversity

# If using sorensen dissimilarity:
# sim = functional turnover
# sne = functional nestedness
# sor = total functional beta diversity

# Functional beta diversity can be used to assess differences in functional composition among treatments

# total beta-diversity
fun.nmds <- metaMDS(fun.beta$funct.beta.jac, trymax = 500, autotransform = TRUE, k = 2)
fun.nmds
plot(fun.nmds)

adonis2(fun.beta$funct.beta.jac ~ env$trmt, permutations = 999)

ordiplot(fun.nmds, disp = "sites", type = "n", xlim = c(-0.5, 1), ylim = c(-0.5, 0.5))
points(fun.nmds, disp = "sites", select = which(env$trmt=="T1"), pch = 17, cex = 2, col = "#27AD81FF")
points(fun.nmds, disp = "sites", select = which(env$trmt=="T2"), pch = 16, cex = 2, col = "#481567FF")
ordiellipse(fun.nmds, env$trmt, draw = "lines", col = c("#27AD81FF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)
legend("bottomright", legend = c("Vacant Lot", "Urban Meadow"),
       pch = c(17, 16), cex = 1.5, bty = "n", col = c("#27AD81FF", "#481567FF"))

#####################################################################################
# Partial Least Squares (PLS) Analysis

# Check out this paper for descriptions of methods
# https://doi.org/10.1007/s11252-020-01069-0

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mixOmics")
library(mixOmics)

## community trait means
local.env <- pls(env[, 9:13], cwm.obs[,1:27], mode = c("canonical"), 
              ncomp = 2, scale = TRUE, max.iter = 100)
local.env

local.env$prop_expl_var
local.env$loadings # identify a cutoff level

plotVar(local.env)
plotLoadings(local.env)

## I will usually run a reduced pls with variables that made the determined cutoff value
## I have used 0.3, 0.4, and 0.5 as cutoff values, depending on the data set

# for this example, let's use 0.3
local.env.red <- pls(env[, 9:13], cwm.obs[, c(2,5,9,13:16,19:22)], mode = c("canonical"), 
                  ncomp = 2, scale = TRUE, max.iter = 100)
local.env.red

plotVar(local.env.red)
local.env.red$loadings
local.env.red$prop_expl_var

# then run a network analysis
# this analysis provides similarity values which are similar to pearson correlation coefficients
# keep the cutoff threshold at 0.5 for this analysis
nw.local.env <- network(local.env.red, cutoff = 0.5, color.edge = color.spectral(2), 
                         lty.edge = c("solid", "dashed"), lwd.edge = 2)
nw.local.env


# save the loadings as a new data frame to make a better figure
local.env.red.pred <- as.data.frame(local.env.red$loadings$X)
local.env.red.resp <- as.data.frame(local.env.red$loadings$Y)

plot(local.env.red.pred$comp2 ~ local.env.red.pred$comp1, pch = 1, ylim = c(-1, 1), xlim = c(-1, 1), col = "black",
     xlab = "PLS Axis 1", ylab = "PLS Axis 2", main = "Local", cex = 2)
points(local.env.red.resp$comp2 ~ local.env.red.resp$comp1, pch = 15, cex = 2)
abline(h = 0.0, v = 0.0, col = "black", lwd = 1, lty=1)

# add text for predictor variables
text(0.27, 0.69, "PLI", pos = 4, font = 1, cex = 1)

# add text for response variables
text(0.07, 0.44, "Aboveground Nests", pos = 4, font = 2, cex = 1)
text(-0.07, -0.44, "Belowground Nests", pos = 2, font = 2, cex = 1)

#####################################################################################
# RLQ and Fourth Corner Analysis

# check out this paper for description of methods
# https://doi.org/10.1002/eap.2191

# we are missing trait data for some species, so have to trim the data set
t.na <- na.omit(t)
trait.Q <- dudi.hillsmith(t.na, row.w = abund.L$cw, scannf = FALSE)

select.species <- rownames(t) %in% rownames (t.na)
a.na <- a[,select.species]
abund.L <- dudi.coa(a.na, scannf = FALSE)

env.R <- dudi.pca(env[,2:13], row.w = abund.L$lw, scannf = FALSE) # use with all continuous variables
env.R <- dudi.hillsmith(env[,2:13], row.w = abund.L$lw, scannf = FALSE) # use with categorical variables

rlq.ants.env <- rlq(env.R, abund.L, trait.Q, scannf = FALSE)

summary(rlq.ants.env)
print(rlq.ants.env)

## weighted correlations axes / env.
t(env.R$tab)%*%(diag(env.R$lw))%*%as.matrix(rlq.ants.env$mR)

## weighted correlations axes / traits.
t(trait.Q$tab)%*%(diag(trait.Q$lw))%*%as.matrix(rlq.ants.env$mQ)

## biplot representing traits and environmental variables
s.arrow(rlq.ants.env$c1, xlim=c(-0.8,0.8), boxes = FALSE)
s.label(rlq.ants.env$l1, add.plot=T, clab=1.0)

## correlations traits / env.
rlq.ants.env$tab
rlq.ants.env$aR
rlq.ants.env$aQ
rlq.ants.env$l1
rlq.ants.env$c1
rlq.ants.env$li
rlq.ants.env$co

plot(rlq.ants.env)
s.arrow(rlq.ants.env$l1)
s.arrow(rlq.ants.env$c1)
s.label(rlq.ants.env$lQ, boxes = FALSE)

## Fourth-corner analysis
nrepet <- 999

##P-values not adjusted for multiple comparisons
four.comb.ant.env <- fourthcorner(env, a.na, t.na, modeltype = 6, p.adjust.method.G = "none",
                                     p.adjust.method.D = "none", nrepet = nrepet)
plot(four.comb.ant.env, alpha = 0.05, stat = "D2")
four.comb.ant.env

##P-values adjusted for multiple comparisons
four.comb.ant.env.adj <- fourthcorner(env, a.na, t.na, modeltype = 6, p.adjust.method.G = "fdr",
                                         p.adjust.method.D = "fdr", nrepet = nrepet)
plot(four.comb.ant.env.adj, alpha = 0.1, stat = "D2")
four.comb.ant.env.adj

## Combining both approaches
## Evaluates the global significance of the traits-environment relationships, based on total
## inertia of the RLQ analysis
testrlq.ant.env <- randtest(rlq.ants.env, modeltype = 6, nrepet = nrepet)
testrlq.ant.env
plot(testrlq.ant.env)

## Calculate the total inertia of RLQ analysis (Srlq)
Srlq.ant.env <- fourthcorner2(env, a.na, t.na,
                                 modeltype = 6, p.adjust.method.G = "fdr", nrepet = nrepet)
Srlq.ant.env$trRLQ

plot(four.comb.ant.env, x.rlq = rlq.ants.env, alpha = 0.1,
     stat = "D2", type = "biplot")
plot(four.comb.ant.env.adj, x.rlq = rlq.ants.env, alpha = 0.05,
     stat = "D2", type = "biplot")

##tests the links between RLQ axes and traits (Qaxes) or environmental variables (Raxes)
testQaxes.comb.ant.env <- fourthcorner.rlq(rlq.ants.env, modeltype = 6,
                                              typetest = "Q.axes", nrepet = nrepet, p.adjust.method.G = "fdr",
                                              p.adjust.method.D = "fdr")
print(testQaxes.comb.ant.env, stat = "D2")

testRaxes.comb.ant.env <- fourthcorner.rlq(rlq.ants.env, modeltype = 6,
                                              typetest = "R.axes", nrepet = nrepet, p.adjust.method.G = "fdr",
                                              p.adjust.method.D = "fdr")
print(testRaxes.comb.ant.env, stat = "D2")

par(mfrow = c(1, 2))
plot(testQaxes.comb.ant.env, alpha = 0.1, type = "table",
     stat = "D2")
plot(testRaxes.comb.ant.env, alpha = 0.1, type = "table",
     stat = "D2")

par(mfrow = c(1, 2))
plot(testQaxes.comb.ant.env, alpha = 0.1, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))
plot(testRaxes.comb.ant.env, alpha = 0.1, type = "biplot",
     stat = "D2", col = c("black", "blue", Holly M. Martinson, 
Michael J. Raupp
Ecospher"orange", "green"))

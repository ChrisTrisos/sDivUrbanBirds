###########################################################################
#######  Testing changes in biodiversity along urbanisation gradients #####
###########################################################################


### GOALS: ###

# Model how the biodiversity metrics vary with urbanisation
    # Quadratic Entropy (QE) for all morphological and behavioral traits + phylogeny
    # Uniqueness* for all morphological and behavioral traits + phylogeny
    # Redundancies (CR = 1-Uniqueness) for all morphological and behavioral traits + phylogeny
    # Simpson index (QE taxonomy)
    # Species richness


### INPUTS: ### 

# Morphologically-based diversity metrics for communities ("Morphological diversity metrics for communities.txt")
# Diet-based diversity metrics for communities ("Morphological diversity metrics for communities.txt")
# Foraging behaviour data-based diversity metrics for communities ("Morphological diversity metrics for communities.txt")


### OUTPUTS: ###

# Tables showing the results of the models
# Figures describing the FD changes along urbanisation gradients

# We will repeat all analyses 1) for native species, 2) exluding exotics and 3) exluding exotics and strict exploiters (natives not present or rare in the surroundings).
# Each of the above analyses will be repeated with 1) all the habitats as well as coding urban habitats as 2) urban garden, suburbs and urban centre and 3) pooling all the categories together  


### ANALYSES START HERE ### 

## libraries

library(nlme)
library(phia)
library(effects)
library(ggplot2)



########### Analyses for morphological diversity ##############
###############################################################

## Download biodiversity metrics for communities

setwd("~/ownCloud2/Science/Research/Urbanisation/Functional diversity and urbanization/Data and analyses")

x <- read.table("Morphological diversity metrics for communities.txt") # metrics estimated including all species
# x <- read.table("Morphological diversity metrics for communities natives.txt")  # metrics estimated excluding exotics


# Three ways to code habitats

# 1. All habitats separated
# levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")

# 2. All urban habitats separated, all non-urban habitat together
levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")


# 3. All urban and non-urban habitats pooled together
# levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")



## Tests of the effect of urbanization on functional diversity*

# * Note that the best random structure have been previously evaluated with the method=REML
# * We need to include confounds (city age, human density, coordinates)


# Simpson's index

spp.richness = lme(Species.richness ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(spp.richness)
plot(spp.richness)
qqnorm(spp.richness)
spp.richness.I <- interactionMeans(spp.richness) # effect plots
plot(spp.richness.I, errorbar="ci95")

QE.taxonomy = lme(QE.taxonomy ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(QE.taxonomy)
plot(QE.taxonomy)
qqnorm(QE.taxonomy)
QE.taxonomy.I <- interactionMeans(QE.taxonomy) # effect plots
plot(QE.taxonomy.I, errorbar="ci95")


# 8 functional traits 

QE.all.morph = lme(QE.all.morph ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(QE.all.morph)
plot(QE.all.morph)
qqnorm(QE.all.morph)
QE.all.morph.I <- interactionMeans(QE.all.morph) # effect plots
plot(QE.all.morph.I, errorbar="ci95")

CR.all.morph = lme(CR.all.morph ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(CR.all.morph)  # urban habitats have lower redundancy
plot(CR.all.morph)
qqnorm(CR.all.morph)
CR.all.morph.I <- interactionMeans(CR.all.morph) # effect plots
plot(CR.all.morph.I, errorbar="ci95")


# All three PCAs (body size, beak shape and locomotory shape)

QE.PCA3 = lme(QE.PCA3 ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(QE.PCA3)
plot(QE.PCA3)
qqnorm(QE.PCA3)
QE.PCA3.I <- interactionMeans(QE.PCA3) # effect plots
plot(QE.PCA3.I, errorbar="ci95")

CR.PCA3 = lme(CR.PCA3 ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(CR.PCA3)  # urban habitats have lower redundancy
plot(CR.PCA3)
qqnorm(CR.PCA3)
CR.PCA3.I <- interactionMeans(CR.PCA3) # effect plots
plot(CR.PCA3.I, errorbar="ci95")


# beak shape

QE.beak = lme(QE.beak ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(QE.beak)
plot(QE.beak)
qqnorm(QE.beak)
QE.beak.I <- interactionMeans(QE.beak) # effect plots
plot(QE.beak.I, errorbar="ci95")

CR.beak = lme(CR.beak ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(CR.beak)  # urban habitats have lower redundancy
plot(CR.beak)
qqnorm(CR.beak)
CR.beak.I <- interactionMeans(CR.beak) # effect plots
plot(CR.beak.I, errorbar="ci95")


# locomotory shape (tarsus vs tail shape)

QE.locom = lme(QE.locom ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(QE.locom)
plot(QE.locom)
qqnorm(QE.locom)
QE.locom.I <- interactionMeans(QE.locom) # effect plots
plot(QE.locom.I, errorbar="ci95")

CR.locom = lme(CR.locom ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(CR.locom)  # urban habitats have lower redundancy
plot(CR.locom)
qqnorm(CR.locom)
CR.locom.I <- interactionMeans(CR.locom) # effect plots
plot(CR.locom.I, errorbar="ci95")


# size

QE.size = lme(QE.size ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(QE.size)
plot(QE.size)
qqnorm(QE.size)
QE.size.I <- interactionMeans(QE.size) # effect plots
plot(QE.size.I, errorbar="ci95")

CR.size = lme(CR.size ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(CR.size)  # urban habitats have lower redundancy
plot(CR.size)
qqnorm(CR.size)
CR.size.I <- interactionMeans(CR.size) # effect plots
plot(CR.size.I, errorbar="ci95")


# winghand

QE.winghand = lme(QE.winghand ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(QE.winghand)
plot(QE.winghand)
qqnorm(QE.winghand)
QE.winghand.I <- interactionMeans(QE.winghand) # effect plots
plot(QE.winghand.I, errorbar="ci95")

CR.winghand = lme(CR.winghand ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(CR.winghand)  # urban habitats have lower redundancy
plot(CR.winghand)
qqnorm(CR.winghand)
CR.winghand.I <- interactionMeans(CR.winghand) # effect plots
plot(CR.winghand.I, errorbar="ci95")


# Phylogeny Ericksson

QE.phyE = lme(QE.phyE ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(QE.phyE)
plot(QE.phyE)
qqnorm(QE.phyE)
QE.phyE.I <- interactionMeans(QE.phyE) # effect plots
plot(QE.phyE.I, errorbar="ci95")

CR.phyE = lme(CR.phyE ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(CR.phyE)  # urban habitats have lower redundancy
plot(CR.phyE)
qqnorm(CR.phyE)
CR.phyE.I <- interactionMeans(CR.phyE) # effect plots
plot(CR.phyE.I, errorbar="ci95")



# Phylogeny Hackett

QE.phyH = lme(QE.phyH ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(QE.phyH)
plot(QE.phyH)
qqnorm(QE.phyH)
QE.phyH.I <- interactionMeans(QE.phyH) # effect plots
plot(QE.phyH.I, errorbar="ci95")

CR.phyH = lme(CR.phyH ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(CR.phyH)  # urban habitats have lower redundancy
plot(CR.phyH)
qqnorm(CR.phyH)
CR.phyH.I <- interactionMeans(CR.phyH) # effect plots
plot(CR.phyH.I, errorbar="ci95")


###########################################################################
#######  Testing changes in biodiversity along urbanisation gradients #####
###########################################################################


### GOALS: ###

# Model how the biodiversity metrics of section 6 vary with urbanisation
  
  # Quadratic Entropy (QE) for all morphological and behavioral traits + phylogeny
  # Uniqueness* for all morphological and behavioral traits + phylogeny
  # Redundancies (CR = 1-Uniqueness) for all morphological and behavioral traits + phylogeny
  # Simpson index (QE taxonomy or HGS)
  # Species richness
  # The meanD
  # The balance factor

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


########### Analyses for morphological diversity ##############
###############################################################

## Download biodiversity metrics for communities

x <- read.table("/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis/Morphological diversity metrics for communities.txt") # metrics estimated including all species
# x <- read.table("Morphological diversity metrics for communities natives.txt")  # metrics estimated excluding exotics


# Three ways to code habitats

# 1. All habitats separated
# levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")

# 2. All urban habitats separated, all non-urban habitat together
levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
# habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))


# 3. All urban and non-urban habitats pooled together
# levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")



## Tests of the effect of urbanization on biodiversity*
################################################################

# * Note that the best random structure has been previously evaluated with the method=REML
# * We need to include confounds (city age, human density, coordinates)



# Species richness

spp.richness = lme(Species.richness ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(spp.richness)
plot(spp.richness)
qqnorm(spp.richness)
spp.richness.I <- interactionMeans(spp.richness) # effect plots
plot(spp.richness.I, errorbar="ci95")


# Simpson's index

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

meanD.all.morph = lme(all.morph.meanD ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.all.morph)  # urban habitats have lower redundancy
plot(meanD.all.morph)
qqnorm(meanD.all.morph)
meanD.all.morph.I <- interactionMeans(meanD.all.morph) # effect plots
plot(meanD.all.morph.I, errorbar="ci95")

Balance.all.morph = lme(all.morph.Balance ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.all.morph)  # urban habitats have lower redundancy
plot(Balance.all.morph)
qqnorm(Balance.all.morph)
Balance.all.morph.I <- interactionMeans(Balance.all.morph) # effect plots
plot(Balance.all.morph.I, errorbar="ci95")


# All three PCAs (body size, size shape and locomotory shape)

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

meanD.PCA3 = lme(PCA3.meanD ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.PCA3)  # urban habitats have lower redundancy
plot(meanD.PCA3)
qqnorm(meanD.PCA3)
meanD.PCA3.I <- interactionMeans(meanD.PCA3) # effect plots
plot(meanD.PCA3.I, errorbar="ci95")

Balance.PCA3 = lme(PCA3.Balance ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.PCA3)  # urban habitats have lower redundancy
plot(Balance.PCA3)
qqnorm(Balance.PCA3)
Balance.PCA3.I <- interactionMeans(Balance.PCA3) # effect plots
plot(Balance.PCA3.I, errorbar="ci95")


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

meanD.beak = lme(beak.meanD ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.beak)  # urban habitats have lower redundancy
plot(meanD.beak)
qqnorm(meanD.beak)
meanD.beak.I <- interactionMeans(meanD.beak) # effect plots
plot(meanD.beak.I, errorbar="ci95")

Balance.beak = lme(beak.Balance ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.beak)  # urban habitats have lower redundancy
plot(Balance.beak)
qqnorm(Balance.beak)
Balance.beak.I <- interactionMeans(Balance.beak) # effect plots
plot(Balance.beak.I, errorbar="ci95")


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

meanD.locom = lme(locom.meanD ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.locom)  # urban habitats have lower redundancy
plot(meanD.locom)
qqnorm(meanD.locom)
meanD.locom.I <- interactionMeans(meanD.locom) # effect plots
plot(meanD.locom.I, errorbar="ci95")

Balance.locom = lme(locom.Balance ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.locom)  # urban habitats have lower redundancy
plot(Balance.locom)
qqnorm(Balance.locom)
Balance.locom.I <- interactionMeans(Balance.locom) # effect plots
plot(Balance.locom.I, errorbar="ci95")


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

meanD.size = lme(size.meanD ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.size)  # urban habitats have lower redundancy
plot(meanD.size)
qqnorm(meanD.size)
meanD.size.I <- interactionMeans(meanD.size) # effect plots
plot(meanD.size.I, errorbar="ci95")

Balance.size = lme(size.Balance ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.size)  # urban habitats have lower redundancy
plot(Balance.size)
qqnorm(Balance.size)
Balance.size.I <- interactionMeans(Balance.size) # effect plots
plot(Balance.size.I, errorbar="ci95")



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

meanD.winghand = lme(winghand.meanD ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.winghand)  # urban habitats have lower redundancy
plot(meanD.winghand)
qqnorm(meanD.winghand)
meanD.winghand.I <- interactionMeans(meanD.winghand) # effect plots
plot(meanD.winghand.I, errorbar="ci95")

Balance.winghand = lme(winghand.Balance ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.winghand)  # urban habitats have lower redundancy
plot(Balance.winghand)
qqnorm(Balance.winghand)
Balance.winghand.I <- interactionMeans(Balance.winghand) # effect plots
plot(Balance.winghand.I, errorbar="ci95")


# Phylogenetic diversity Ericksson

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

meanD.phyE = lme(phyE.meanD ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.phyE)  # urban habitats have lower redundancy
plot(meanD.phyE)
qqnorm(meanD.phyE)
meanD.phyE.I <- interactionMeans(meanD.phyE) # effect plots
plot(meanD.phyE.I, errorbar="ci95")

Balance.phyE = lme(phyE.Balance ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.phyE)  # urban habitats have lower redundancy
plot(Balance.phyE)
qqnorm(Balance.phyE)
Balance.phyE.I <- interactionMeans(Balance.phyE) # effect plots
plot(Balance.phyE.I, errorbar="ci95")



# Phylogenetic diversity Hackett

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

meanD.phyH = lme(phyH.meanD ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.phyH)  # urban habitats have lower redundancy
plot(meanD.phyH)
qqnorm(meanD.phyH)
meanD.phyH.I <- interactionMeans(meanD.phyH) # effect plots
plot(meanD.phyH.I, errorbar="ci95")

Balance.phyH = lme(phyH.Balance ~ habitat, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.phyH)  # urban habitats have lower redundancy
plot(Balance.phyH)
qqnorm(Balance.phyH)
Balance.phyH.I <- interactionMeans(Balance.phyH) # effect plots
plot(Balance.phyH.I, errorbar="ci95")
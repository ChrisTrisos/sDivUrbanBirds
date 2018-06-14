
####################################################
#######  ANALYSES FINALLY USED IN THE PAPER  #######
####################################################


### LOSS OF DIVERSITY ###
#########################


## All analyses are for native species


########### 1. Analyses of species loss  ##############
#######################################################

{

# data
  
x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities ocurrences natives.txt"))
x <- subset(x, used.urban.nonurban=="yes")
x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
x <- cbind(x,habitat.ordered)


# Differences in species richness across habitats

spp.richness = lme(Species.richness ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(spp.richness)
testInteractions(spp.richness, pairwise="habitat.ordered", adjustment="BH")
plot(spp.richness)
qqnorm(spp.richness)
spp.richness.I <- interactionMeans(spp.richness) # effect plots
plot(spp.richness.I, errorbar="ci95")

}


###########   2. Analyses of Simpson index and abundance eveness  ##############
################################################################################

{
x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities natives.txt"))

# We restrict the analyses to studies where there is information inside and outside the city
x <- subset(x, used.urban.nonurban=="yes")
x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
x <- cbind(x,habitat.ordered)


# Simpson's index

QE.taxonomy = lme(QE.taxonomy ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.taxonomy)
testInteractions(QE.taxonomy, pairwise="habitat.ordered", adjustment="BH")
plot(QE.taxonomy)
qqnorm(QE.taxonomy)
QE.taxonomy.I <- interactionMeans(QE.taxonomy) # effect plots
plot(QE.taxonomy.I, errorbar="ci95")


# Simpson*R/(R-1), an index of abundance evenness, independent of species richness. 

Abundance.Eveness <- x$QE.taxonomy*x$Species.richness/(x$Species.richness-1)
x <- cbind(x,Abundance.Eveness)

AEveness = lme(Abundance.Eveness ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(AEveness)
testInteractions(AEveness, pairwise="habitat.ordered", adjustment="BH")
plot(AEveness)
qqnorm(AEveness)
AEveness.I <- interactionMeans(AEveness) # effect plots
plot(AEveness.I, errorbar="ci95")

}


########### 3. Analyses for morphological diversity ##############
##################################################################

{

  x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities natives.txt"))
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  
  
  # 8 functional traits 
  
  QE.all.morph = lme(QE.all.morph ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.all.morph)
  testInteractions(QE.all.morph, pairwise="habitat.ordered", adjustment="BH")
  plot(QE.all.morph)
  qqnorm(QE.all.morph)
  QE.all.morph.I <- interactionMeans(QE.all.morph) # effect plots
  plot(QE.all.morph.I, errorbar="ci95")
  
  CR.all.morph = lme(CR.all.morph ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.all.morph)  
  testInteractions(CR.all.morph, pairwise="habitat.ordered", adjustment="BH")
  plot(CR.all.morph)
  qqnorm(CR.all.morph)
  CR.all.morph.I <- interactionMeans(CR.all.morph) # effect plots
  plot(CR.all.morph.I, errorbar="ci95")
  
  meanD.all.morph = lme(all.morph.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.all.morph)  
  testInteractions(meanD.all.morph, pairwise="habitat.ordered", adjustment="BH")
  plot(meanD.all.morph)
  qqnorm(meanD.all.morph)
  meanD.all.morph.I <- interactionMeans(meanD.all.morph) # effect plots
  plot(meanD.all.morph.I, errorbar="ci95")
  
  Balance.all.morph = lme(all.morph.Balance.cor ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(Balance.all.morph)  
  testInteractions(Balance.all.morph, pairwise="habitat.ordered", adjustment="BH")
  plot(Balance.all.morph)
  qqnorm(Balance.all.morph)
  Balance.all.morph.I <- interactionMeans(Balance.all.morph) # effect plots
  plot(Balance.all.morph.I, errorbar="ci95")
  
}



########### 4. Analyses for diet diversity ##########
#####################################################

{
  x<-read.table(paste0(workingData,"/Diet diversity metrics for communities natives.txt"))
   x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)

  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
 
  QE.diet = lme(QE.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.diet)
  testInteractions(QE.diet, pairwise="habitat.ordered", adjustment="BH")
  plot(QE.diet)
  qqnorm(QE.diet)
  QE.diet.I <- interactionMeans(QE.diet) # effect plots
  plot(QE.diet.I, errorbar="ci95")
  
  CR.diet = lme(CR.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.diet)  
  testInteractions(CR.diet, pairwise="habitat.ordered", adjustment="BH")
  plot(CR.diet)
  qqnorm(CR.diet)
  CR.diet.I <- interactionMeans(CR.diet) # effect plots
  plot(CR.diet.I, errorbar="ci95")
  
  meanD.diet = lme(diet.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.diet)  
  testInteractions(meanD.diet, pairwise="habitat.ordered", adjustment="BH")
  plot(meanD.diet)
  qqnorm(meanD.diet)
  meanD.diet.I <- interactionMeans(meanD.diet) # effect plots
  plot(meanD.diet.I, errorbar="ci95")
  
  Balance.diet = lme(diet.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
  summary(Balance.diet)  
  testInteractions(Balance.diet, pairwise="habitat.ordered", adjustment="BH")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
  plot(Balance.diet)
  qqnorm(Balance.diet)
  Balance.diet.I <- interactionMeans(Balance.diet) # effect plots
  plot(Balance.diet.I, errorbar="ci95")
  
  
}


########### 5. Analyses for foraging diversity #########
########################################################

{
  x<-read.table(paste0(workingData,"/Forag diversity metrics for communities natives.txt"))
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  
  QE.forag = lme(QE.forag ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.forag)
  testInteractions(QE.forag, pairwise="habitat.ordered", adjustment="BH")
  plot(QE.forag)
  qqnorm(QE.forag)
  QE.forag.I <- interactionMeans(QE.forag) # effect plots
  plot(QE.forag.I, errorbar="ci95")
  
  CR.forag = lme(CR.forag ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.forag)  
  testInteractions(CR.forag, pairwise="habitat.ordered", adjustment="BH")
  plot(CR.forag)
  qqnorm(CR.forag)
  CR.forag.I <- interactionMeans(CR.forag) # effect plots
  plot(CR.forag.I, errorbar="ci95")
  
  meanD.forag = lme(forag.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.forag)  
  testInteractions(meanD.forag, pairwise="habitat.ordered", adjustment="BH")
  plot(meanD.forag)
  qqnorm(meanD.forag)
  meanD.forag.I <- interactionMeans(meanD.forag) # effect plots
  plot(meanD.forag.I, errorbar="ci95")
  
  Balance.forag = lme(forag.Balance.cor ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
  summary(Balance.forag)  
  testInteractions(Balance.forag, pairwise="habitat.ordered", adjustment="BH")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
  plot(Balance.forag)
  qqnorm(Balance.forag)
  Balance.forag.I <- interactionMeans(Balance.forag) # effect plots
  plot(Balance.forag.I, errorbar="ci95")
  
  
}



########### 6. Analyses for species-level functional uniqueness ################
################################################################################

{

# data
  
dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)
dat$community <- factor(dat$community)
diet.Kbar <- read.table(paste0(workingData,"/vulnerability_species_diet_natives.txt"))
morphology.Kbar <- read.table(paste0(workingData,"/vulnerability_species_morphology_natives.txt"))
forag.Kbar <- read.table(paste0(workingData,"/vulnerability_species_foraging_natives.txt"))


library(moments)

tmp.hab <- ddply(dat, c("country", "city", "community", "habitat", "used.urban.nonurban"), summarise,
             Regional.richness = length(relative.abundance))

tmp.diet <- ddply(diet.Kbar, c("community"), summarise,
                  Kbar.diet.skew = skewness(Kbar.diet))

tmp2 <- merge(tmp.diet,tmp.hab, by="community")

tmp.morph <- ddply(morphology.Kbar, c("community"), summarise,
                   Kbar.morph.skew = skewness(Kbar.morphology))

x <- merge(tmp.morph,tmp2, by="community")

tmp.forag <- ddply(forag.Kbar, c("community"), summarise,
                   Kbar.forag.skew = skewness(Kbar.foraging))

x <- merge(tmp.forag,x, by="community")


levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
x <- cbind(x,habitat.ordered)


# diet

Kbar.diet.skew = lme(Kbar.diet.skew ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(Kbar.diet.skew)
r2<- rsquared(Kbar.diet.skew)  #  marginal R^2 (based on fixed effects only) and conditional R^2 (based on fixed and random effects, if present)
testInteractions(Kbar.diet.skew, pairwise="habitat.ordered", adjustment="BH")
plot(Kbar.diet.skew)
qqnorm(Kbar.diet.skew)
Kbar.diet.skew.I <- interactionMeans(Kbar.diet.skew) # effect plots
plot(Kbar.diet.skew, errorbar="ci95")
hist(residuals(Kbar.diet.skew))


# Morphology

Kbar.morph.skew = lme(Kbar.morph.skew ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(Kbar.morph.skew)
r2<- rsquared(Kbar.morph.skew)  #  marginal R^2 (based on fixed effects only) and conditional R^2 (based on fixed and random effects, if present)
testInteractions(Kbar.morph.skew, pairwise="habitat.ordered", adjustment="BH")
plot(Kbar.morph.skew)
qqnorm(Kbar.morph.skew)
Kbar.morph.skew.I <- interactionMeans(Kbar.morph.skew) # effect plots
plot(Kbar.morph.skew, errorbar="ci95")
hist(residuals(Kbar.morph.skew))


# Foraging

Kbar.forag.skew = lme(Kbar.forag.skew ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(Kbar.forag.skew)
r2<- rsquared(Kbar.forag.skew)  #  marginal R^2 (based on fixed effects only) and conditional R^2 (based on fixed and random effects, if present)
testInteractions(Kbar.forag.skew, pairwise="habitat.ordered", adjustment="BH")
plot(Kbar.forag.skew)
qqnorm(Kbar.forag.skew)
Kbar.forag.skew.I <- interactionMeans(Kbar.forag.skew) # effect plots
plot(Kbar.forag.skew, errorbar="ci95")
hist(residuals(Kbar.forag.skew))

}


########### 7. Analyses rarefied species ################
#########################################################


{
  dat<-read.table(paste0(workingData,"/Rarefied.communities.txt"), header=TRUE)
  levels(dat$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(dat$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  dat <- cbind(dat,habitat.ordered)
  
  
  # Species richness rarefied
  
  spp.richnessR = lme(SRraref ~ habitat.ordered, random = ~ 1|country/city, data = dat, method="ML")
  summary(spp.richnessR)
  testInteractions(spp.richnessR, pairwise="habitat.ordered")
  plot(spp.richnessR)
  qqnorm(spp.richnessR)
  spp.richnessR.I <- interactionMeans(spp.richnessR) # effect plots
  plot(spp.richnessR.I, errorbar="ci95")
}



########### 8. Analyses of community-weighted mean morphological traits ################
#######################################################################################


{
  
  x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities natives.txt"))
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  
  
  # traits
  
  MCWM.beak.shape = lme(CWM.beak.shape ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWM.beak.shape)
  testInteractions(MCWM.beak.shape, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWM.beak.shape)
  qqnorm(MCWM.beak.shape)
  MCWM.beak.shape.I <- interactionMeans(MCWM.beak.shape) # effect plots
  plot(MCWM.beak.shape.I, errorbar="ci95")
  
  MCWM.locom.shape = lme(CWM.locom.shape ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWM.locom.shape)
  testInteractions(MCWM.locom.shape, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWM.locom.shape)
  qqnorm(MCWM.locom.shape)
  MCWM.locom.shape.I <- interactionMeans(MCWM.locom.shape) # effect plots
  plot(MCWM.locom.shape.I, errorbar="ci95")
  
  MCWM.body.size = lme(CWM.body.size ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWM.body.size)
  testInteractions(MCWM.body.size, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWM.body.size)
  qqnorm(MCWM.body.size)
  MCWM.body.size.I <- interactionMeans(MCWM.body.size) # effect plots
  plot(MCWM.body.size.I, errorbar="ci95")
  
  MCWM.hand.wing = lme(CWM.hand.wing ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWM.hand.wing)
  testInteractions(MCWM.hand.wing, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWM.hand.wing)
  qqnorm(MCWM.hand.wing)
  MCWM.hand.wing.I <- interactionMeans(MCWM.hand.wing) # effect plots
  plot(MCWM.hand.wing.I, errorbar="ci95")
  
  
  
  
}


########### 9. Analyses of community-weighted mean diet traits ################
##############################################################################

{
  x<-read.table(paste0(workingData,"/Diet diversity metrics for communities natives.txt"))
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  
  MCWMInv_SeedFruit= lme(CWMInv_SeedFruit~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWMInv_SeedFruit)
  testInteractions(MCWMInv_SeedFruit, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWMInv_SeedFruit)
  qqnorm(MCWMInv_SeedFruit)
  MCWMInv_SeedFruit.I <- interactionMeans(MCWMInv_SeedFruit) # effect plots
  plot(MCWMInv_SeedFruit.I, errorbar="ci95")
  
  MCWMInv_Seed= lme(CWMInv_Seed~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWMInv_Seed)
  testInteractions(MCWMInv_Seed, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWMInv_Seed)
  qqnorm(MCWMInv_Seed)
  MCWMInv_Seed.I <- interactionMeans(MCWMInv_Seed) # effect plots
  plot(MCWMInv_Seed.I, errorbar="ci95")
  
  MCWMSeed= lme(CWMSeed~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWMSeed)
  testInteractions(MCWMSeed, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWMSeed)
  qqnorm(MCWMSeed)
  MCWMSeed.I <- interactionMeans(MCWMSeed) # effect plots
  plot(MCWMSeed.I, errorbar="ci95")

  
  }


########### 10. Analyses of community-weighted mean foraging traits ################
###################################################################################

{
  x<-read.table(paste0(workingData,"/Forag diversity metrics for communities natives.txt"))
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  
  MCWM_PCo1= lme(CWM_PCo1~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWM_PCo1)
  testInteractions(MCWM_PCo1, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWM_PCo1)
  qqnorm(MCWM_PCo1)
  MCWM_PCo1.I <- interactionMeans(MCWM_PCo1) # effect plots
  plot(MCWM_PCo1.I, errorbar="ci95")
  
  MCWM_PCo2= lme(CWM_PCo2~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWM_PCo2)
  testInteractions(MCWM_PCo2, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWM_PCo2)
  qqnorm(MCWM_PCo2)
  MCWM_PCo2.I <- interactionMeans(MCWM_PCo2) # effect plots
  plot(MCWM_PCo2.I, errorbar="ci95")
  
  MCWM_PCo3= lme(CWM_PCo3~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWM_PCo3)
  testInteractions(MCWM_PCo3, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWM_PCo3)
  qqnorm(MCWM_PCo3)
  MCWM_PCo3.I <- interactionMeans(MCWM_PCo3) # effect plots
  plot(MCWM_PCo3.I, errorbar="ci95")
  
  MCWM_PCo4= lme(CWM_PCo4~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(MCWM_PCo4)
  testInteractions(MCWM_PCo4, pairwise="habitat.ordered", adjustment="BH")
  plot(MCWM_PCo4)
  qqnorm(MCWM_PCo4)
  MCWM_PCo4.I <- interactionMeans(MCWM_PCo4) # effect plots
  plot(MCWM_PCo4.I, errorbar="ci95")
  
}


########## 11. Morphological PCAs #########
###########################################

{
  
  morph0<- read.table(paste0(workingData,"/Morphological traits urban birds 24 Feb 2018 for R.txt"), h=TRUE)
  morph <- morph0[,c(8,10,12:18)] # we exclude redundant traits
  names(morph) <- c("mass", "bill_length", "bill_width", "bill_depth", "tarsus", "second","wing","hand_wing","tail")
  
  # Log10-transformed and standardized to zero mean and unit variance (z-transformation).
  
  morph$mass <- log(morph$mass)
  morph$bill_length <- log(morph$bill_length)
  morph$bill_width <- log(morph$bill_width)
  morph$bill_depth <- log(morph$bill_depth)
  morph$tarsus <- log(morph$tarsus)
  morph$second <- log(morph$second)
  morph$wing <- log(morph$wing)
  morph$hand_wing <- log(morph$hand_wing)
  morph$tail <- log(morph$tail)
  
  morph <- as.data.frame(scale(morph))   # we scale the data
  # colnames(morph) <- c("Body mass", "Bill length", "Bill width", "Bill depth", "Tarsus length","Second wing length","Wing length","Hand wing index","Tail length")       
  colnames(morph) <- c("BM", "BL", "BW", "BD", "TS","SW","WL","HWI","TL")       

  
  pca.9 <- dudi.pca(morph[,-8], scannf = F, nf = 3)  # We exclude hand_wing index
  Ppca.9 <- scatter(pca.9, clab.row=1, clab.col = 1.5, xax = 1, yax = 2, posieig = "none", ylab="PCA 2", xlab="PCA 1")
  
  beak <- morph[,c(2:4)]
  pca.beak <- dudi.pca(beak, scannf = F, nf = 3)
  Pbeak <- scatter(pca.beak, clab.row=0, clab.col = 1.5, xax = 1, yax = 2, posieig = "none", ylab="PCA 2", xlab="PCA 1")
  
  locom <- morph[,c(5:7,9)]
  pca.locom <- dudi.pca(locom, scannf = F, nf = 3)
  Plocom <- scatter(pca.locom, clab.row=0, clab.col = 1.5, xax = 1, yax = 2, posieig = "none", ylab="PCA 2", xlab="PCA 1")
  
  
  
}



#######################
#### Other analyses ###
#######################

### Alpha Q diversity with occurrence data ###
##############################################

{

########### Analyses for morphological diversity natives, ocurrence data ######
##################################################################################


{
  x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities ocurrences natives.txt"))
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  # 8 functional traits 
  
  QE.all.morphO = lme(QE.all.morph ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.all.morphO)
  testInteractions(QE.all.morphO, pairwise="habitat.ordered", adjustment="BH")
  plot(QE.all.morphO)
  qqnorm(QE.all.morphO)
  QE.all.morphO.I <- interactionMeans(QE.all.morphO) # effect plots
  plot(QE.all.morphO.I, errorbar="ci95")
  
}



########### Analyses for diet diversity ocurrence natives ###############
#############################################################################

{
  x<-read.table(paste0(workingData,"/Diet diversity metrics for communities ocurrences natives.txt"))
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  
  QE.dietO = lme(QE.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.dietO)
  testInteractions(QE.dietO, pairwise="habitat.ordered", adjustment="BH")
  plot(QE.dietO)
  qqnorm(QE.dietO)
  QE.dietO.I <- interactionMeans(QE.dietO) # effect plots
  plot(QE.dietO.I, errorbar="ci95")
  
  
}

  
  ########### Analyses for foraging diversity ocurrence natives ###############
  #############################################################################
  
{
    x<-read.table(paste0(workingData,"/Forag diversity metrics for communities ocurrences natives.txt"))
    x <- subset(x, used.urban.nonurban=="yes")
    x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
    
    levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
    habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
    x <- cbind(x,habitat.ordered)
    
    
    QE.foragO = lme(QE.forag ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
    summary(QE.foragO)
    testInteractions(QE.foragO, pairwise="habitat.ordered", adjustment="BH")
    plot(QE.foragO)
    qqnorm(QE.foragO)
    QE.foragO.I <- interactionMeans(QE.foragO) # effect plots
    plot(QE.foragO.I, errorbar="ci95")
    
    
  }
  

}


### Beta Q diversity ###
########################

{
### Morphological Beta Q test

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Rao beta diversity morphology natives.txt"))
  
  y <- subset(x, hab.comp!="Urban_Urban_Park" & hab.comp!="Urban_Park_Urban_Park" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Urban_Park_Wildland" & hab.comp!="Urban_Park_Rural" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Wildland_Wildland" & hab.comp!="Suburban_Suburban" & hab.comp!="Urban_Urban" & hab.comp!="Rural_Rural")
  y[] <- lapply(y, function(x) if(is.factor(x)) factor(x) else x)
  
  habitat.ordered  = factor(y$hab.comp, levels=c("Rural_Wildland","Suburban_Rural","Suburban_Wildland","Urban_Wildland","Urban_Rural","Urban_Suburban"))
  y <- cbind(y,habitat.ordered)
  
  habitat.ordered  = factor(y$hab.comp, levels=c("Rural_Wildland","Suburban_Rural","Suburban_Wildland","Urban_Wildland","Urban_Rural","Urban_Suburban"))
  y <- cbind(y,habitat.ordered)
  
  
  
  ## Tests of the effect of urbanization on beta biodiversity*
  ################################################################
  
  
  # Changes in functional composition
  
  Qbetast.diff.morph = lme(Qbetast ~ habitat.ordered, random = ~ 1|city, data = y, method="ML")
  summary(Qbetast.diff.morph)
  testInteractions(Qbetast.diff.morph, pairwise="habitat.ordered")
  plot(Qbetast.diff.morph)
  qqnorm(Qbetast.diff.morph)
  Qbetast.diff.morph.I <- interactionMeans(Qbetast.diff.morph) # effect plots
  plot(Qbetast.diff.morph.I, errorbar="ci95")
  
  
  
  # Final figure (effect +/- standard error)
  
  ab <- ggplot(Qbetast.diff.morph.I, aes(x= habitat.ordered, y=Qbetast.diff.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Qbetast.diff.morph.I[,2]-Qbetast.diff.morph.I[,3], ymax=Qbetast.diff.morph.I[,2]+Qbetast.diff.morph.I[,3]), width=.2) +
    geom_point(data=Qbetast.diff.morph.I, mapping=aes(x=habitat.ordered, y=Qbetast.diff.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Beta Q morphology", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b","a")))
  
} 


### Foraging Beta Q test

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Rao beta diversity foraging natives.txt"))
  
  y <- subset(x, hab.comp!="Urban_Urban_Park" & hab.comp!="Urban_Park_Urban_Park" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Urban_Park_Wildland" & hab.comp!="Urban_Park_Rural" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Wildland_Wildland" & hab.comp!="Suburban_Suburban" & hab.comp!="Urban_Urban" & hab.comp!="Rural_Rural")
  y[] <- lapply(y, function(x) if(is.factor(x)) factor(x) else x)
  
  habitat.ordered  = factor(y$hab.comp, levels=c("Rural_Wildland","Suburban_Rural","Suburban_Wildland","Urban_Wildland","Urban_Rural","Urban_Suburban"))
  y <- cbind(y,habitat.ordered)
  
  ## Tests of the effect of urbanization on beta biodiversity*
  ################################################################
  
  
  # Changes in functional composition
  
  Qbetast.diff.forag = lme(Qbetast ~ habitat.ordered, random = ~ 1|city, data = y, method="ML")
  summary(Qbetast.diff.forag)
  testInteractions(Qbetast.diff.forag, pairwise="habitat.ordered")
  plot(Qbetast.diff.forag)
  qqnorm(Qbetast.diff.forag)
  Qbetast.diff.forag.I <- interactionMeans(Qbetast.diff.forag) # effect plots
  plot(Qbetast.diff.forag.I, errorbar="ci95")
  
  
  
  # Final figure (effect +/- standard error)
  
  bc <- ggplot(Qbetast.diff.forag.I, aes(x= habitat.ordered, y=Qbetast.diff.forag.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Qbetast.diff.forag.I[,2]-Qbetast.diff.forag.I[,3], ymax=Qbetast.diff.forag.I[,2]+Qbetast.diff.forag.I[,3]), width=.2) +
    geom_point(data=Qbetast.diff.forag.I, mapping=aes(x=habitat.ordered, y=Qbetast.diff.forag.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Beta Q foraging", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b","a")))
  
}


### Diet Beta Q test

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Rao beta diversity diet natives.txt"))
  
  y <- subset(x, hab.comp!="Urban_Urban_Park" & hab.comp!="Urban_Park_Urban_Park" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Urban_Park_Wildland" & hab.comp!="Urban_Park_Rural" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Wildland_Wildland" & hab.comp!="Suburban_Suburban" & hab.comp!="Urban_Urban" & hab.comp!="Rural_Rural")
  y[] <- lapply(y, function(x) if(is.factor(x)) factor(x) else x)
  
  habitat.ordered  = factor(y$hab.comp, levels=c("Rural_Wildland","Suburban_Rural","Suburban_Wildland","Urban_Wildland","Urban_Rural","Urban_Suburban"))
  y <- cbind(y,habitat.ordered)
  
  ## Tests of the effect of urbanization on beta biodiversity*
  ################################################################
  
  
  # Changes in functional composition
  
  Qbetast.diff.diet = lme(Qbetast ~ habitat.ordered, random = ~ 1|city, data = y, method="ML")
  summary(Qbetast.diff.diet)
  testInteractions(Qbetast.diff.diet, pairwise="habitat.ordered")
  plot(Qbetast.diff.diet)
  qqnorm(Qbetast.diff.diet)
  Qbetast.diff.diet.I <- interactionMeans(Qbetast.diff.diet) # effect plots
  plot(Qbetast.diff.diet.I, errorbar="ci95")
  
  
  
  # Final figure (effect +/- standard error)
  
  cd <- ggplot(Qbetast.diff.diet.I, aes(x= habitat.ordered, y=Qbetast.diff.diet.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Qbetast.diff.diet.I[,2]-Qbetast.diff.diet.I[,3], ymax=Qbetast.diff.diet.I[,2]+Qbetast.diff.diet.I[,3]), width=.2) +
    geom_point(data=Qbetast.diff.diet.I, mapping=aes(x=habitat.ordered, y=Qbetast.diff.diet.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Beta Q diet", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b","a")))
  
}

}


###########  All Figures  ################
##########################################

{
# Species richness

a <- ggplot(spp.richness.I, aes(x= habitat.ordered, y=spp.richness.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=spp.richness.I[,2]-spp.richness.I[,3], ymax=spp.richness.I[,2]+spp.richness.I[,3]), width=.2) +
  geom_point(data=spp.richness.I, mapping=aes(x=habitat.ordered, y=spp.richness.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Species richness", cex=16) +
  geom_text(aes(label= c("a","a","a","a","b")))


# Simpson

b <- ggplot(QE.taxonomy.I, aes(x= habitat.ordered, y=QE.taxonomy.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.taxonomy.I[,2]-QE.taxonomy.I[,3], ymax=QE.taxonomy.I[,2]+QE.taxonomy.I[,3]), width=.2) +
  geom_point(data=QE.taxonomy.I, mapping=aes(x=habitat.ordered, y=QE.taxonomy.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Gini-Simpson's index", cex=16) +
  geom_text(aes(label= c("a","a","a","b","c")))


# Eveness

c <- ggplot(AEveness.I, aes(x= habitat.ordered, y=AEveness.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=AEveness.I[,2]-AEveness.I[,3], ymax=AEveness.I[,2]+AEveness.I[,3]), width=.2) +
  geom_point(data=AEveness.I, mapping=aes(x=habitat.ordered, y=AEveness.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Abundance evenness index", cex=16) +
  geom_text(aes(label= c("a","a","a","b","b")))


# Morphology, abundance weighted

d <- ggplot(QE.all.morph.I, aes(x= habitat.ordered, y=QE.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.all.morph.I[,2]-QE.all.morph.I[,3], ymax=QE.all.morph.I[,2]+QE.all.morph.I[,3]), width=.2) +
  geom_point(data=QE.all.morph.I, mapping=aes(x=habitat.ordered, y=QE.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Q morphology", cex=16) +
  geom_text(aes(label= c("a","ab","ab","b","b")))

e <- ggplot(CR.all.morph.I, aes(x= habitat.ordered, y=CR.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.all.morph.I[,2]-CR.all.morph.I[,3], ymax=CR.all.morph.I[,2]+CR.all.morph.I[,3]), width=.2) +
  geom_point(data=CR.all.morph.I, mapping=aes(x=habitat.ordered, y=CR.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Morphological redundancies", cex=16) +
  geom_text(aes(label= c("ab","a","ab","ab","b")))

f <- ggplot(meanD.all.morph.I, aes(x= habitat.ordered, y=meanD.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.all.morph.I[,2]-meanD.all.morph.I[,3], ymax=meanD.all.morph.I[,2]+meanD.all.morph.I[,3]), width=.2) +
  geom_point(data=meanD.all.morph.I, mapping=aes(x=habitat.ordered, y=meanD.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Morphological meanD", cex=16) +
  geom_text(aes(label= c("a","b","ab","ab","a")))

g <- ggplot(Balance.all.morph.I, aes(x= habitat.ordered, y=Balance.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.all.morph.I[,2]-Balance.all.morph.I[,3], ymax=Balance.all.morph.I[,2]+Balance.all.morph.I[,3]), width=.2) +
  geom_point(data=Balance.all.morph.I, mapping=aes(x=habitat.ordered, y=Balance.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Morphological BC", cex=16) +
  geom_text(aes(label= c("a","a","a","a","b")))


# Figure decomposition Q morphology

QE.taxonomy.I[,4] <- "2:HGS"
QE.all.morph.I[,4] <- "1:Q"
meanD.all.morph.I[,4] <- "3:meanD"
Balance.all.morph.I[,4] <- "4:BC"
tgc2 <- rbind(QE.taxonomy.I,meanD.all.morph.I,Balance.all.morph.I)

colnames(tgc2) <- c("habitat.ordered","adjusted.mean","std.error","Metric")


decomp.morph <- ggplot(tgc2, aes(x=habitat.ordered, y=adjusted.mean, fill=Metric)) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  labs(x = "", y = "Effect size, morphology", cex=16) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=adjusted.mean-std.error, ymax=adjusted.mean+std.error),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))



# diet, abundance weighted

h <- ggplot(QE.diet.I , aes(x= habitat.ordered, y=QE.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.diet.I [,2]-QE.diet.I [,3], ymax=QE.diet.I [,2]+QE.diet.I [,3]), width=.2) +
  geom_point(data=QE.diet.I , mapping=aes(x=habitat.ordered, y=QE.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Q diet", cex=16) +
  geom_text(aes(label= c("a","a","a","a","b")))

i <- ggplot(CR.diet.I , aes(x= habitat.ordered, y=CR.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.diet.I [,2]-CR.diet.I [,3], ymax=CR.diet.I [,2]+CR.diet.I [,3]), width=.2) +
  geom_point(data=CR.diet.I , mapping=aes(x=habitat.ordered, y=CR.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Diet redundancies", cex=16) +
  geom_text(aes(label= c("a","a","ab","bc","c")))

j <- ggplot(meanD.diet.I, aes(x= habitat.ordered, y=meanD.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.diet.I [,2]-meanD.diet.I [,3], ymax=meanD.diet.I [,2]+meanD.diet.I [,3]), width=.2) +
  geom_point(data=meanD.diet.I , mapping=aes(x=habitat.ordered, y=meanD.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Diet meanD ", cex=16) +
  geom_text(aes(label= c("a","ab","ab","ab","b")))

k <- ggplot(Balance.diet.I , aes(x= habitat.ordered, y=Balance.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.diet.I [,2]-Balance.diet.I [,3], ymax=Balance.diet.I [,2]+Balance.diet.I [,3]), width=.2) +
  geom_point(data=Balance.diet.I , mapping=aes(x=habitat.ordered, y=Balance.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Diet BC", cex=16) +
  geom_text(aes(label= c("ab","ab","ab","a","b")))



# Figure decomposition Q diet

QE.taxonomy.I[,4] <- "2:HGS"
QE.diet.I[,4] <- "1:Q"
meanD.diet.I[,4] <- "3:meanD"
Balance.diet.I[,4] <- "4:BC"
tgc3 <- rbind(QE.taxonomy.I,meanD.diet.I,Balance.diet.I)

colnames(tgc3) <- c("habitat.ordered","adjusted.mean","std.error","Metric")


decomp.diet <- ggplot(tgc3, aes(x=habitat.ordered, y=adjusted.mean, fill=Metric)) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  labs(x = "", y = "Effect size, diet", cex=16) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=adjusted.mean-std.error, ymax=adjusted.mean+std.error),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))





# foraging, abundance weighted

h1 <- ggplot(QE.forag.I , aes(x= habitat.ordered, y=QE.forag.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.forag.I [,2]-QE.forag.I [,3], ymax=QE.forag.I [,2]+QE.forag.I [,3]), width=.2) +
  geom_point(data=QE.forag.I , mapping=aes(x=habitat.ordered, y=QE.forag.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Q foraging", cex=16) +
  geom_text(aes(label= c("a","a","b","a","c")))

i1 <- ggplot(CR.forag.I , aes(x= habitat.ordered, y=CR.forag.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.forag.I [,2]-CR.forag.I [,3], ymax=CR.forag.I [,2]+CR.forag.I [,3]), width=.2) +
  geom_point(data=CR.forag.I , mapping=aes(x=habitat.ordered, y=CR.forag.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Foraging redundancies", cex=16) +
  geom_text(aes(label= c("a","a","ab","b","b")))

j1 <- ggplot(meanD.forag.I, aes(x= habitat.ordered, y=meanD.forag.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.forag.I [,2]-meanD.forag.I [,3], ymax=meanD.forag.I [,2]+meanD.forag.I [,3]), width=.2) +
  geom_point(data=meanD.forag.I , mapping=aes(x=habitat.ordered, y=meanD.forag.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Foraging meanD ", cex=16) +
  geom_text(aes(label= c("ab","a","ab","ab","b")))

k1 <- ggplot(Balance.forag.I , aes(x= habitat.ordered, y=Balance.forag.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.forag.I [,2]-Balance.forag.I [,3], ymax=Balance.forag.I [,2]+Balance.forag.I [,3]), width=.2) +
  geom_point(data=Balance.forag.I , mapping=aes(x=habitat.ordered, y=Balance.forag.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Foraging BC", cex=16) +
  geom_text(aes(label= c("a","a","ab","b","ab")))



# Figure decomposition Q forag

QE.taxonomy.I[,4] <- "2:HGS"
QE.forag.I[,4] <- "1:Q"
meanD.forag.I[,4] <- "3:meanD"
Balance.forag.I[,4] <- "4:BC"
tgc3 <- rbind(QE.taxonomy.I,meanD.forag.I,Balance.forag.I)

colnames(tgc3) <- c("habitat.ordered","adjusted.mean","std.error","Metric")


decomp.forag <- ggplot(tgc3, aes(x=habitat.ordered, y=adjusted.mean, fill=Metric)) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  labs(x = "", y = "Effect size, foraging", cex=16) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=adjusted.mean-std.error, ymax=adjusted.mean+std.error),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))




# Morphology, occurrences

l <- ggplot(QE.all.morphO.I, aes(x= habitat.ordered, y=QE.all.morphO.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.all.morphO.I[,2]-QE.all.morphO.I[,3], ymax=QE.all.morphO.I[,2]+QE.all.morphO.I[,3]), width=.2) +
  geom_point(data=QE.all.morphO.I, mapping=aes(x=habitat.ordered, y=QE.all.morphO.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Q morphology, unweigthed", cex=16) +
  geom_text(aes(label= c("a","b","ab","ab","ab")))


# Diet, occurrences

m <- ggplot(QE.dietO.I, aes(x= habitat.ordered, y=QE.dietO.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.dietO.I[,2]-QE.dietO.I[,3], ymax=QE.dietO.I[,2]+QE.dietO.I[,3]), width=.2) +
  geom_point(data=QE.dietO.I, mapping=aes(x=habitat.ordered, y=QE.dietO.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Q diet, unweigthed", cex=16) +
  geom_text(aes(label= c("a","ab","ab","ab","b")))



# Foraging, occurrences

m1 <- ggplot(QE.foragO.I, aes(x= habitat.ordered, y=QE.foragO.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.foragO.I[,2]-QE.foragO.I[,3], ymax=QE.foragO.I[,2]+QE.foragO.I[,3]), width=.2) +
  geom_point(data=QE.foragO.I, mapping=aes(x=habitat.ordered, y=QE.foragO.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Q foraging, unweigthed", cex=16) +
  geom_text(aes(label= c("ab","ab","a","b","c")))



# Kbar morphology

n <- ggplot(Kbar.morph.skew.I, aes(x= habitat.ordered, y=Kbar.morph.skew.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Kbar.morph.skew.I[,2]-Kbar.morph.skew.I[,3], ymax=Kbar.morph.skew.I[,2]+Kbar.morph.skew.I[,3]), width=.2) +
  geom_point(data=Kbar.morph.skew.I, mapping=aes(x=habitat.ordered, y=Kbar.morph.skew.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Kbar morphology", cex=16) +
  geom_text(aes(label= c("a","a","ab","a","b")))

# Kbar diet

o <- ggplot(Kbar.diet.skew.I, aes(x= habitat.ordered, y=Kbar.diet.skew.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Kbar.diet.skew.I[,2]-Kbar.diet.skew.I[,3], ymax=Kbar.diet.skew.I[,2]+Kbar.diet.skew.I[,3]), width=.2) +
  geom_point(data=Kbar.diet.skew.I, mapping=aes(x=habitat.ordered, y=Kbar.diet.skew.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Kbar diet", cex=16) +
  geom_text(aes(label= c("a","a","a","b","b")))


# Kbar foraging

o1 <- ggplot(Kbar.forag.skew.I, aes(x= habitat.ordered, y=Kbar.forag.skew.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Kbar.forag.skew.I[,2]-Kbar.forag.skew.I[,3], ymax=Kbar.forag.skew.I[,2]+Kbar.forag.skew.I[,3]), width=.2) +
  geom_point(data=Kbar.forag.skew.I, mapping=aes(x=habitat.ordered, y=Kbar.forag.skew.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Kbar foraging", cex=16) +
  geom_text(aes(label= c("a","a","ab","ab","b")))




# SR rarefied

p <- ggplot(spp.richnessR.I, aes(x= habitat.ordered, y=spp.richnessR.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=spp.richnessR.I[,2]-spp.richnessR.I[,3], ymax=spp.richnessR.I[,2]+spp.richnessR.I[,3]), width=.2) +
  geom_point(data=spp.richnessR.I, mapping=aes(x=habitat.ordered, y=spp.richnessR.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Species richness rarefied", cex=16) +
  geom_text(aes(label= c("a","ab","ab","bc","c")))




# MCWM.beak.shape

q <- ggplot(MCWM.beak.shape.I, aes(x= habitat.ordered, y=MCWM.beak.shape.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWM.beak.shape.I[,2]-MCWM.beak.shape.I[,3], ymax=MCWM.beak.shape.I[,2]+MCWM.beak.shape.I[,3]), width=.2) +
  geom_point(data=MCWM.beak.shape.I, mapping=aes(x=habitat.ordered, y=MCWM.beak.shape.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM beak shape", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))



# MCWM.locom.shape

r <- ggplot(MCWM.locom.shape.I, aes(x= habitat.ordered, y=MCWM.locom.shape.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWM.locom.shape.I[,2]-MCWM.locom.shape.I[,3], ymax=MCWM.locom.shape.I[,2]+MCWM.locom.shape.I[,3]), width=.2) +
  geom_point(data=MCWM.locom.shape.I, mapping=aes(x=habitat.ordered, y=MCWM.locom.shape.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM tarsus-to-tail ratio", cex=16) +
  geom_text(aes(label= c("a","a","b","b","c")))


# MCWM.body.size

s <- ggplot(MCWM.body.size.I, aes(x= habitat.ordered, y=MCWM.body.size.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWM.body.size.I[,2]-MCWM.body.size.I[,3], ymax=MCWM.body.size.I[,2]+MCWM.body.size.I[,3]), width=.2) +
  geom_point(data=MCWM.body.size.I, mapping=aes(x=habitat.ordered, y=MCWM.body.size.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM body size", cex=16) +
  geom_text(aes(label= c("a","a","bc","ab","c")))


# MCWM.hand.wing

u <- ggplot(MCWM.hand.wing.I, aes(x= habitat.ordered, y=MCWM.hand.wing.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWM.hand.wing.I[,2]-MCWM.hand.wing.I[,3], ymax=MCWM.hand.wing.I[,2]+MCWM.hand.wing.I[,3]), width=.2) +
  geom_point(data=MCWM.hand.wing.I, mapping=aes(x=habitat.ordered, y=MCWM.hand.wing.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM hand-wing", cex=16) +
  geom_text(aes(label= c("a","a","b","b","c")))


# MCWMInv_SeedFruit

v <- ggplot(MCWMInv_SeedFruit.I, aes(x= habitat.ordered, y=MCWMInv_SeedFruit.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWMInv_SeedFruit.I[,2]-MCWMInv_SeedFruit.I[,3], ymax=MCWMInv_SeedFruit.I[,2]+MCWMInv_SeedFruit.I[,3]), width=.2) +
  geom_point(data=MCWMInv_SeedFruit.I, mapping=aes(x=habitat.ordered, y=MCWMInv_SeedFruit.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM diet PCo1", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))


# MCWMInv_Seed

w <- ggplot(MCWMInv_Seed.I, aes(x= habitat.ordered, y=MCWMInv_Seed.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWMInv_Seed.I[,2]-MCWMInv_Seed.I[,3], ymax=MCWMInv_Seed.I[,2]+MCWMInv_Seed.I[,3]), width=.2) +
  geom_point(data=MCWMInv_Seed.I, mapping=aes(x=habitat.ordered, y=MCWMInv_Seed.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM diet PCo2", cex=16) +
  geom_text(aes(label= c("a","ab","abc","bc","c")))


# MCWMSeed

x <- ggplot(MCWMSeed.I, aes(x= habitat.ordered, y=MCWMSeed.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWMSeed.I[,2]-MCWMSeed.I[,3], ymax=MCWMSeed.I[,2]+MCWMSeed.I[,3]), width=.2) +
  geom_point(data=MCWMSeed.I, mapping=aes(x=habitat.ordered, y=MCWMSeed.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM diet PCo3", cex=16) +
  geom_text(aes(label= c("a","b","c","b","b")))


# Foraging PCo1

x0 <- ggplot(MCWM_PCo1.I, aes(x= habitat.ordered, y=MCWM_PCo1.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWM_PCo1.I[,2]-MCWM_PCo1.I[,3], ymax=MCWM_PCo1.I[,2]+MCWM_PCo1.I[,3]), width=.2) +
  geom_point(data=MCWM_PCo1.I, mapping=aes(x=habitat.ordered, y=MCWM_PCo1.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM foraging PCo1", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))


# Foraging PCo2

x1 <- ggplot(MCWM_PCo2.I, aes(x= habitat.ordered, y=MCWM_PCo2.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWM_PCo2.I[,2]-MCWM_PCo2.I[,3], ymax=MCWM_PCo2.I[,2]+MCWM_PCo2.I[,3]), width=.2) +
  geom_point(data=MCWM_PCo2.I, mapping=aes(x=habitat.ordered, y=MCWM_PCo2.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM foraging PCo2", cex=16) +
  geom_text(aes(label= c("a","a","ab","ab","b")))


# Foraging PCo3

x2 <- ggplot(MCWM_PCo3.I, aes(x= habitat.ordered, y=MCWM_PCo3.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWM_PCo3.I[,2]-MCWM_PCo3.I[,3], ymax=MCWM_PCo3.I[,2]+MCWM_PCo3.I[,3]), width=.2) +
  geom_point(data=MCWM_PCo3.I, mapping=aes(x=habitat.ordered, y=MCWM_PCo3.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM foraging PCo3", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))


# Foraging PCo4

x3 <- ggplot(MCWM_PCo4.I, aes(x= habitat.ordered, y=MCWM_PCo4.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=MCWM_PCo4.I[,2]-MCWM_PCo4.I[,3], ymax=MCWM_PCo4.I[,2]+MCWM_PCo4.I[,3]), width=.2) +
  geom_point(data=MCWM_PCo4.I, mapping=aes(x=habitat.ordered, y=MCWM_PCo4.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CWM foraging PCo4", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))


}



###########  Main Figures  ################
###########################################



# Q loss with urbanization

tiff(paste0(GoogleFigs,"/Fig. 1.tiff"), width = 11, height = 8, units = 'in', res = 50)
ggplot2.multiplot(d,h,h1, cols=1)
dev.off()




# Q decomposition

tiff(paste0(GoogleFigs,"/Fig. 2.tiff"), width = 16, height = 8, units = 'in', res = 50)
ggplot2.multiplot(decomp.morph, decomp.diet, decomp.forag, cols=1) 
dev.off()



# Beta diversity

tiff(paste0(GoogleFigs,"/Fig. 3.tiff"), width = 11, height = 8, units = 'in', res = 100)
ggplot2.multiplot(ab,  cd, bc, cols=1)
dev.off()








###########  Extended data Figures  #################
#####################################################

# Species richness

tiff(paste0(GoogleFigs,"/Extended data Fig. 2.tiff"), width = 11, height = 8, units = 'in', res = 100)
ggplot2.multiplot(a,p,c,b)
dev.off()


# Q for occurrence data

tiff(paste0(GoogleFigs,"/Extended data Fig. 3.tiff"), width = 11, height = 8, units = 'in', res = 50)
ggplot2.multiplot(l,m,m1, cols=1)
dev.off()


# Redundancies and Kbar

tiff(paste0(GoogleFigs,"/Extended data Fig. 4.tiff"), width = 16, height = 8, units = 'in', res = 50)
ggplot2.multiplot(e,i,i1,n,o,o1, cols=3) 
dev.off()



# PCA
Ppca.9
Pbeak
Plocom


# Community weighted means

tiff(paste0(GoogleFigs,"/Extended data Fig. 5.tiff"), width = 11, height = 8, units = 'in', res = 100)
ggplot2.multiplot(s, q,r,u,v,w,x0,x1,x2, cols=3)
dev.off()




# details on Q decomposition

tiff(paste0(GoogleFigs,"/Extended data Fig. 5.tiff"), width = 11, height = 8, units = 'in', res = 100)
ggplot2.multiplot(f,j,j1,g,k,k1, cols=3)
dev.off()




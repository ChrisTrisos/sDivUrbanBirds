
####################################################
#######  ANALYSES FINALLY USED IN THE PAPER  #######
####################################################


### LOSS OF DIVERSITY ###
#########################


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


########### 2. Analyses for morphological diversity, only natives ##############
################################################################################

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
  
  Balance.all.morph = lme(all.morph.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(Balance.all.morph)  
  testInteractions(Balance.all.morph, pairwise="habitat.ordered", adjustment="BH")
  plot(Balance.all.morph)
  qqnorm(Balance.all.morph)
  Balance.all.morph.I <- interactionMeans(Balance.all.morph) # effect plots
  plot(Balance.all.morph.I, errorbar="ci95")
  
}



########### 3. Analyses for diet diversity natives ###############
##################################################################

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



########### 4. Analyses for morphological diversity natives, ocurrence data ######
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



########### 5. Analyses for diet diversity ocurrence natives ###############
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



########### 6. Analyses for species-level functional uniqueness natives ################
########################################################################################

{

# data
  
dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)
dat$community <- factor(dat$community)
diet.Kbar <- read.table(paste0(workingData,"/vulnerability_species_diet_natives.txt"))
morphology.Kbar <- read.table(paste0(workingData,"/vulnerability_species_morphology_natives.txt"))


library(moments)

tmp <- ddply(dat, c("country", "city", "community", "habitat", "used.urban.nonurban"), summarise,
             Regional.Kbar.diet.skew = length(relative.abundance))

tmp.diet <- ddply(diet.Kbar, c("community"), summarise,
                  Kbar.diet.skew = skewness(Kbar.diet))

tmp2 <- merge(tmp.diet,tmp, by="community")

tmp.morph <- ddply(morphology.Kbar, c("community"), summarise,
                   Kbar.morph.skew = skewness(Kbar.morphology))

x <- merge(tmp.morph,tmp2, by="community")

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



########### 8. Analyses of commuity-weighted mean morphological traits ################
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


########### 9. Analyses of commuity-weighted mean diet traits ################
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
  labs(x = "", y = "Simpson's index", cex=16) +
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
  labs(x = "", y = "Morphological redundancy", cex=16) +
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
  geom_text(aes(label= c("a","a","a","a","a")))


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
  geom_text(aes(label= c("a","a","a","a","a")))

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
  geom_text(aes(label= c("a","a","a","a","b")))


# Kbar morphology

n <- ggplot(Kbar.morph.skew.I, aes(x= habitat.ordered, y=Kbar.morph.skew.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Kbar.morph.skew.I[,2]-Kbar.morph.skew.I[,3], ymax=Kbar.morph.skew.I[,2]+Kbar.morph.skew.I[,3]), width=.2) +
  geom_point(data=Kbar.morph.skew.I, mapping=aes(x=habitat.ordered, y=Kbar.morph.skew.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Kbar skeweness, morphology", cex=16) +
  geom_text(aes(label= c("a","a","ab","a","b")))

# Kbar diet

o <- ggplot(Kbar.diet.skew.I, aes(x= habitat.ordered, y=Kbar.diet.skew.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Kbar.diet.skew.I[,2]-Kbar.diet.skew.I[,3], ymax=Kbar.diet.skew.I[,2]+Kbar.diet.skew.I[,3]), width=.2) +
  geom_point(data=Kbar.diet.skew.I, mapping=aes(x=habitat.ordered, y=Kbar.diet.skew.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Kbar skeweness, diet", cex=16) +
  geom_text(aes(label= c("a","a","a","b","b")))


# SR rarefied

p <- ggplot(spp.richnessR.I, aes(x= habitat.ordered, y=spp.richnessR.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=spp.richnessR.I[,2]-spp.richnessR.I[,3], ymax=spp.richnessR.I[,2]+spp.richnessR.I[,3]), width=.2) +
  geom_point(data=spp.richnessR.I, mapping=aes(x=habitat.ordered, y=spp.richnessR.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Species richness rarefied", cex=16) +
  geom_text(aes(label= c("a","ab","ab","bc","c")))


}



###########  Figure 1  ################
#######################################


tiff(paste0(GoogleFigs,"/Fig. 1.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(d,h,l,m,cols=2)
dev.off()


tiff(paste0(GoogleFigs,"/Fig. 2.tiff"), width = 16, height = 12, units = 'in', res = 200)
ggplot2.multiplot(e,i,f,j,g,k,n,o, cols=3) 
dev.off()


###########  Extended data Figure 1  ################
#####################################################

# Species richness

tiff(paste0(GoogleFigs,"/Extended data Fig. 1.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,p,c,b)
dev.off()

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

# Morphologically-based diversity metrics for communities, based on abundances ("Morphological diversity metrics for communities.txt")
# Morphologically-based diversity metrics for communities including only natives, based on abundances ("Morphological diversity metrics for communities.txt")
# Morphologically-based diversity metrics for communities, based on occurrences (Morphological diversity metrics for communities ocurrences.txt)


### OUTPUTS: ###

# Tables showing the results of the models
# Figures describing the FD changes along urbanisation gradients

# We will repeat all analyses 1) for native species, 2) exluding exotics and 3) exluding exotics and strict exploiters (natives not present or rare in the surroundings).
# Each of the above analyses will be repeated with 1) all the habitats as well as coding urban habitats as 2) urban garden, suburbs and urban centre and 3) pooling all the categories together  


### ANALYSES: ###

## 1. Analyses for morphological diversity
## 2. Analyses for morphological diversity, only natives
## 3. Analyses for morphological diversity, with occurrence data
## 4. Analyses for diet diversity
## 5. Analyses for diet diversity, only natives
## 6. Analyses for diet diversity, with occurrence data
## 7. Analyses for diet diversity, with occurrence data natives
## 8. Analyses for morphological diversity natives, with occurrence data
## 9. Analyses for morphology-diet diversity, only natives



########### 1. Analyses for morphological diversity ##############
##################################################################

{## Import biodiversity metrics for communities

x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities.txt"))

# We restrict the analyses to studies where there is information inside and outside the city
x <- subset(x, used.urban.nonurban=="yes")
x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)


# Three ways to code habitats

# 1. All habitats separated
# levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")

# 2. All urban habitats separated, all non-urban habitat together
levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
x <- cbind(x,habitat.ordered)

# 3. All urban and non-urban habitats pooled together
# levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")



## Tests of the effect of urbanization on biodiversity*
################################################################

# * Note that the best random structure has been previously evaluated with the method=REML
# * We need to include confounds (city age, human density, coordinates)



# Species richness

spp.richness = lme(Species.richness ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(spp.richness)
r2<- rsquared(spp.richness)  #  marginal R^2 (based on fixed effects only) and conditional R^2 (based on fixed and random effects, if present)
testInteractions(spp.richness, pairwise="habitat.ordered")
plot(spp.richness)
qqnorm(spp.richness)
spp.richness.I <- interactionMeans(spp.richness) # effect plots
plot(spp.richness.I, errorbar="ci95")
hist(residuals(spp.richness))

# the same with lme4

library(lme4)
library(lmerTest)
library(piecewiseSEM)

m <- lmer(Species.richness ~ habitat.ordered + (1 | country/city), data = x)
summary(m)
anova(m) # with p-values from F-tests using Satterthwaite's denominator df
testInteractions(m, pairwise="habitat.ordered")
(lsm <- ls_means(m))
library(effects)
ef <- effect("habitat.ordered", m)
plot(ef)
y <- as.data.frame(ef)
coef(m) #intercept for each level in Batch 

# including spatial analyses
library(geoR)
spp.richness.ratio.geo<-update(spp.richness, correlation = corRatio(0.01, form= ~ east + south), method="ML") #east and south are variables in the dataset for latitude and longitude 
anova(spp.richness, spp.richness.ratio.geo) #Spatial autocorrelation was included in the final models if the model including spatial #autocorrelation had a significantly better fit to the data than the model with no spatial effect.



# Simpson's index

QE.taxonomy = lme(QE.taxonomy ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.taxonomy)
testInteractions(QE.taxonomy, pairwise="habitat.ordered")
plot(QE.taxonomy)
qqnorm(QE.taxonomy)
QE.taxonomy.I <- interactionMeans(QE.taxonomy) # effect plots
plot(QE.taxonomy.I, errorbar="ci95")


# Simpson*R/(R-1), an index of abundance evenness, independent of species richness. 

Abundance.Eveness <- x$QE.taxonomy*x$Species.richness/(x$Species.richness-1)
x <- cbind(x,Abundance.Eveness)
  
AEveness = lme(Abundance.Eveness ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(AEveness)
testInteractions(AEveness, pairwise="habitat.ordered")
plot(AEveness)
qqnorm(AEveness)
AEveness.I <- interactionMeans(AEveness) # effect plots
plot(AEveness.I, errorbar="ci95")




# Final figure (effect +/- standard error)

a <- ggplot(spp.richness.I, aes(x= habitat.ordered, y=spp.richness.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=spp.richness.I[,2]-spp.richness.I[,3], ymax=spp.richness.I[,2]+spp.richness.I[,3]), width=.2) +
  geom_point(data=spp.richness.I, mapping=aes(x=habitat.ordered, y=spp.richness.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Species richness", cex=16) +
  geom_text(aes(label= c("a","a","a","a","b")))

b <- ggplot(QE.taxonomy.I, aes(x= habitat.ordered, y=QE.taxonomy.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.taxonomy.I[,2]-QE.taxonomy.I[,3], ymax=QE.taxonomy.I[,2]+QE.taxonomy.I[,3]), width=.2) +
  geom_point(data=QE.taxonomy.I, mapping=aes(x=habitat.ordered, y=QE.taxonomy.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Simpson's index", cex=16) +
  geom_text(aes(label= c("a","a","ab","b","c")))

c <- ggplot(AEveness.I, aes(x= habitat.ordered, y=AEveness.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=AEveness.I[,2]-AEveness.I[,3], ymax=AEveness.I[,2]+AEveness.I[,3]), width=.2) +
  geom_point(data=AEveness.I, mapping=aes(x=habitat.ordered, y=AEveness.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Abundance evenness index", cex=16) +
  geom_text(aes(label= c("a","a","a","b","c")))

tiff(paste0(GoogleFigs,"/plot_taxonomic_diversity.tiff"), width = 9, height = 11, units = 'in', res = 200)
ggplot2.multiplot(a,b,c, cols=1)
dev.off()





# 8 functional traits 

QE.all.morph = lme(QE.all.morph ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.all.morph)
testInteractions(QE.all.morph, pairwise="habitat.ordered")
plot(QE.all.morph)
qqnorm(QE.all.morph)
QE.all.morph.I <- interactionMeans(QE.all.morph) # effect plots
plot(QE.all.morph.I, errorbar="ci95")

CR.all.morph = lme(CR.all.morph ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.all.morph)  
testInteractions(CR.all.morph, pairwise="habitat.ordered")
plot(CR.all.morph)
qqnorm(CR.all.morph)
CR.all.morph.I <- interactionMeans(CR.all.morph) # effect plots
plot(CR.all.morph.I, errorbar="ci95")

meanD.all.morph = lme(all.morph.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.all.morph)  
testInteractions(meanD.all.morph, pairwise="habitat.ordered")
plot(meanD.all.morph)
qqnorm(meanD.all.morph)
meanD.all.morph.I <- interactionMeans(meanD.all.morph) # effect plots
plot(meanD.all.morph.I, errorbar="ci95")

Balance.all.morph = lme(all.morph.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.all.morph)  
testInteractions(Balance.all.morph, pairwise="habitat.ordered")
plot(Balance.all.morph)
qqnorm(Balance.all.morph)
Balance.all.morph.I <- interactionMeans(Balance.all.morph) # effect plots
plot(Balance.all.morph.I, errorbar="ci95")


# Final figure (effect +/- standard error)

a <- ggplot(QE.all.morph.I, aes(x= habitat.ordered, y=QE.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.all.morph.I[,2]-QE.all.morph.I[,3], ymax=QE.all.morph.I[,2]+QE.all.morph.I[,3]), width=.2) +
  geom_point(data=QE.all.morph.I, mapping=aes(x=habitat.ordered, y=QE.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional diversity (QE), all traits", cex=16) +
  geom_text(aes(label= c("a","b","ab","b","b")))

b <- ggplot(CR.all.morph.I, aes(x= habitat.ordered, y=CR.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.all.morph.I[,2]-CR.all.morph.I[,3], ymax=CR.all.morph.I[,2]+CR.all.morph.I[,3]), width=.2) +
  geom_point(data=CR.all.morph.I, mapping=aes(x=habitat.ordered, y=CR.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, all traits", cex=16) +
  geom_text(aes(label= c("ab","c","abc","ac","b")))

c <- ggplot(meanD.all.morph.I, aes(x= habitat.ordered, y=meanD.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.all.morph.I[,2]-meanD.all.morph.I[,3], ymax=meanD.all.morph.I[,2]+meanD.all.morph.I[,3]), width=.2) +
  geom_point(data=meanD.all.morph.I, mapping=aes(x=habitat.ordered, y=meanD.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, all traits", cex=16) +
  geom_text(aes(label= c("a","b","ab","ab","a")))
  
d <- ggplot(Balance.all.morph.I, aes(x= habitat.ordered, y=Balance.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.all.morph.I[,2]-Balance.all.morph.I[,3], ymax=Balance.all.morph.I[,2]+Balance.all.morph.I[,3]), width=.2) +
  geom_point(data=Balance.all.morph.I, mapping=aes(x=habitat.ordered, y=Balance.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, all traits", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

tiff(paste0(GoogleFigs,"/plot_FD_all_traits.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()



# beak shape

QE.beak = lme(QE.beak ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.beak)
testInteractions(QE.beak, pairwise="habitat.ordered")
plot(QE.beak)
qqnorm(QE.beak)
QE.beak.I <- interactionMeans(QE.beak) # effect plots
plot(QE.beak.I, errorbar="ci95")

CR.beak = lme(CR.beak ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.beak)  
testInteractions(CR.beak, pairwise="habitat.ordered")
plot(CR.beak)
qqnorm(CR.beak)
CR.beak.I <- interactionMeans(CR.beak) # effect plots
plot(CR.beak.I, errorbar="ci95")

meanD.beak = lme(beak.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.beak)  
testInteractions(meanD.beak, pairwise="habitat.ordered")
plot(meanD.beak)
qqnorm(meanD.beak)
meanD.beak.I <- interactionMeans(meanD.beak) # effect plots
plot(meanD.beak.I, errorbar="ci95")

Balance.beak = lme(beak.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.beak)  
testInteractions(Balance.beak, pairwise="habitat.ordered")
plot(Balance.beak)
qqnorm(Balance.beak)
Balance.beak.I <- interactionMeans(Balance.beak) # effect plots
plot(Balance.beak.I, errorbar="ci95")

# Final figure (effect +/- standard error)

a <- ggplot(QE.beak.I, aes(x= habitat.ordered, y=QE.beak.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.beak.I[,2]-QE.beak.I[,3], ymax=QE.beak.I[,2]+QE.beak.I[,3]), width=.2) +
  geom_point(data=QE.beak.I, mapping=aes(x=habitat.ordered, y=QE.beak.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional diversity (QE), beak", cex=16) +
  geom_text(aes(label= c("a","a","ab","b","b")))

b <- ggplot(CR.beak.I, aes(x= habitat.ordered, y=CR.beak.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.beak.I[,2]-CR.beak.I[,3], ymax=CR.beak.I[,2]+CR.beak.I[,3]), width=.2) +
  geom_point(data=CR.beak.I, mapping=aes(x=habitat.ordered, y=CR.beak.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, beak", cex=16) +
  geom_text(aes(label= c("a","a","a","b","b")))

c <- ggplot(meanD.beak.I, aes(x= habitat.ordered, y=meanD.beak.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.beak.I[,2]-meanD.beak.I[,3], ymax=meanD.beak.I[,2]+meanD.beak.I[,3]), width=.2) +
  geom_point(data=meanD.beak.I, mapping=aes(x=habitat.ordered, y=meanD.beak.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, beak", cex=16) +
  geom_text(aes(label= c("a","a","ab","b","c")))

d <- ggplot(Balance.beak.I, aes(x= habitat.ordered, y=Balance.beak.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.beak.I[,2]-Balance.beak.I[,3], ymax=Balance.beak.I[,2]+Balance.beak.I[,3]), width=.2) +
  geom_point(data=Balance.beak.I, mapping=aes(x=habitat.ordered, y=Balance.beak.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, beak", cex=16) +
  geom_text(aes(label= c("a","a","a","b","b")))

tiff(paste0(GoogleFigs,"/plot_FD_beak.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()


# locomotory shape (tarsus vs tail shape)

QE.locom = lme(QE.locom ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.locom)
testInteractions(QE.locom, pairwise="habitat.ordered")
plot(QE.locom)
qqnorm(QE.locom)
QE.locom.I <- interactionMeans(QE.locom) # effect plots
plot(QE.locom.I, errorbar="ci95")

CR.locom = lme(CR.locom ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.locom)  
testInteractions(CR.locom, pairwise="habitat.ordered")
plot(CR.locom)
qqnorm(CR.locom)
CR.locom.I <- interactionMeans(CR.locom) # effect plots
plot(CR.locom.I, errorbar="ci95")

meanD.locom = lme(locom.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.locom)  
testInteractions(meanD.locom, pairwise="habitat.ordered")
plot(meanD.locom)
qqnorm(meanD.locom)
meanD.locom.I <- interactionMeans(meanD.locom) # effect plots
plot(meanD.locom.I, errorbar="ci95")

Balance.locom = lme(locom.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.locom)  
testInteractions(Balance.locom, pairwise="habitat.ordered")
plot(Balance.locom)
qqnorm(Balance.locom)
Balance.locom.I <- interactionMeans(Balance.locom) # effect plots
plot(Balance.locom.I, errorbar="ci95")


# Final figure (effect +/- standard error)

a <- ggplot(QE.locom.I , aes(x= habitat.ordered, y=QE.locom.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.locom.I [,2]-QE.locom.I [,3], ymax=QE.locom.I [,2]+QE.locom.I [,3]), width=.2) +
  geom_point(data=QE.locom.I , mapping=aes(x=habitat.ordered, y=QE.locom.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional diversity (QE), locomotory system", cex=16) +
  geom_text(aes(label= c("a","ab","ab","ab","b")))

b <- ggplot(CR.locom.I , aes(x= habitat.ordered, y=CR.locom.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.locom.I [,2]-CR.locom.I [,3], ymax=CR.locom.I [,2]+CR.locom.I [,3]), width=.2) +
  geom_point(data=CR.locom.I , mapping=aes(x=habitat.ordered, y=CR.locom.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, locomotory system", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

c <- ggplot(meanD.locom.I , aes(x= habitat.ordered, y=meanD.locom.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.locom.I [,2]-meanD.locom.I [,3], ymax=meanD.locom.I [,2]+meanD.locom.I [,3]), width=.2) +
  geom_point(data=meanD.locom.I , mapping=aes(x=habitat.ordered, y=meanD.locom.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, locomotory system", cex=16) +
  geom_text(aes(label= c("ab","ab","a","b","b")))

d <- ggplot(Balance.locom.I , aes(x= habitat.ordered, y=Balance.locom.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.locom.I [,2]-Balance.locom.I [,3], ymax=Balance.locom.I [,2]+Balance.locom.I [,3]), width=.2) +
  geom_point(data=Balance.locom.I , mapping=aes(x=habitat.ordered, y=Balance.locom.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, locomotory system", cex=16) +
  geom_text(aes(label= c("ab","ab","a","ab","b")))

tiff(paste0(GoogleFigs,"/plot_FD_Locomotory.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()



# size

QE.size = lme(QE.size ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.size)
testInteractions(QE.size, pairwise="habitat.ordered")
plot(QE.size)
qqnorm(QE.size)
QE.size.I <- interactionMeans(QE.size) # effect plots
plot(QE.size.I, errorbar="ci95")

CR.size = lme(CR.size ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.size)  
testInteractions(CR.size, pairwise="habitat.ordered")
plot(CR.size)
qqnorm(CR.size)
CR.size.I <- interactionMeans(CR.size) # effect plots
plot(CR.size.I, errorbar="ci95")

meanD.size = lme(size.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.size)  
testInteractions(meanD.size, pairwise="habitat.ordered")
plot(meanD.size)
qqnorm(meanD.size)
meanD.size.I <- interactionMeans(meanD.size) # effect plots
plot(meanD.size.I, errorbar="ci95")

Balance.size = lme(size.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.size)  
testInteractions(Balance.size, pairwise="habitat.ordered")
plot(Balance.size)
qqnorm(Balance.size)
Balance.size.I <- interactionMeans(Balance.size) # effect plots
plot(Balance.size.I, errorbar="ci95")

# Final figure (effect +/- standard error)

a <- ggplot(QE.size.I , aes(x= habitat.ordered, y=QE.size.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.size.I [,2]-QE.size.I [,3], ymax=QE.size.I [,2]+QE.size.I [,3]), width=.2) +
  geom_point(data=QE.size.I , mapping=aes(x=habitat.ordered, y=QE.size.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional diversity (QE), body size", cex=16) +
  geom_text(aes(label= c("a","b","a","b","b")))

b <- ggplot(CR.size.I , aes(x= habitat.ordered, y=CR.size.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.size.I [,2]-CR.size.I [,3], ymax=CR.size.I [,2]+CR.size.I [,3]), width=.2) +
  geom_point(data=CR.size.I , mapping=aes(x=habitat.ordered, y=CR.size.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, body size", cex=16) +
  geom_text(aes(label= c("a","b","ab","b","a")))

c <- ggplot(meanD.size.I , aes(x= habitat.ordered, y=meanD.size.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.size.I [,2]-meanD.size.I [,3], ymax=meanD.size.I [,2]+meanD.size.I [,3]), width=.2) +
  geom_point(data=meanD.size.I , mapping=aes(x=habitat.ordered, y=meanD.size.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, body size", cex=16) +
  geom_text(aes(label= c("a","b","abc","abc","ac")))

d <- ggplot(Balance.size.I , aes(x= habitat.ordered, y=Balance.size.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.size.I [,2]-Balance.size.I [,3], ymax=Balance.size.I [,2]+Balance.size.I [,3]), width=.2) +
  geom_point(data=Balance.size.I , mapping=aes(x=habitat.ordered, y=Balance.size.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, body size", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

tiff(paste0(GoogleFigs,"/plot_FD_Body_size.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()



# winghand

QE.winghand = lme(QE.winghand ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.winghand)
testInteractions(QE.winghand, pairwise="habitat.ordered")
plot(QE.winghand)
qqnorm(QE.winghand)
QE.winghand.I <- interactionMeans(QE.winghand) # effect plots
plot(QE.winghand.I, errorbar="ci95")

CR.winghand = lme(CR.winghand ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.winghand)  
testInteractions(CR.winghand, pairwise="habitat.ordered")
plot(CR.winghand)
qqnorm(CR.winghand)
CR.winghand.I <- interactionMeans(CR.winghand) # effect plots
plot(CR.winghand.I, errorbar="ci95")

meanD.winghand = lme(winghand.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.winghand)  
testInteractions(meanD.winghand, pairwise="habitat.ordered")
plot(meanD.winghand)
qqnorm(meanD.winghand)
meanD.winghand.I <- interactionMeans(meanD.winghand) # effect plots
plot(meanD.winghand.I, errorbar="ci95")

Balance.winghand = lme(winghand.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.winghand)  
testInteractions(Balance.winghand, pairwise="habitat.ordered")
plot(Balance.winghand)
qqnorm(Balance.winghand)
Balance.winghand.I <- interactionMeans(Balance.winghand) # effect plots
plot(Balance.winghand.I, errorbar="ci95")


# Final figure (effect +/- standard error)

a <- ggplot(QE.winghand.I , aes(x= habitat.ordered, y=QE.winghand.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.winghand.I [,2]-QE.winghand.I [,3], ymax=QE.winghand.I [,2]+QE.winghand.I [,3]), width=.2) +
  geom_point(data=QE.winghand.I , mapping=aes(x=habitat.ordered, y=QE.winghand.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional diversity (QE), wing hand", cex=16) +
  geom_text(aes(label= c("a","ab","ab","ab","b")))

b <- ggplot(CR.winghand.I , aes(x= habitat.ordered, y=CR.winghand.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.winghand.I [,2]-CR.winghand.I [,3], ymax=CR.winghand.I [,2]+CR.winghand.I [,3]), width=.2) +
  geom_point(data=CR.winghand.I , mapping=aes(x=habitat.ordered, y=CR.winghand.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, wing hand ", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

c <- ggplot(meanD.winghand.I , aes(x= habitat.ordered, y=meanD.winghand.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.winghand.I [,2]-meanD.winghand.I [,3], ymax=meanD.winghand.I [,2]+meanD.winghand.I [,3]), width=.2) +
  geom_point(data=meanD.winghand.I , mapping=aes(x=habitat.ordered, y=meanD.winghand.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, wing hand ", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

d <- ggplot(Balance.winghand.I , aes(x= habitat.ordered, y=Balance.winghand.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.winghand.I [,2]-Balance.winghand.I [,3], ymax=Balance.winghand.I [,2]+Balance.winghand.I [,3]), width=.2) +
  geom_point(data=Balance.winghand.I , mapping=aes(x=habitat.ordered, y=Balance.winghand.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, wing hand ", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

tiff(paste0(GoogleFigs,"/plot_FD_wing_hand.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()



# Phylogenetic diversity Ericksson

QE.phyE = lme(QE.phyE ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.phyE)
testInteractions(QE.phyE, pairwise="habitat.ordered")
plot(QE.phyE)
qqnorm(QE.phyE)
QE.phyE.I <- interactionMeans(QE.phyE) # effect plots
plot(QE.phyE.I, errorbar="ci95")

CR.phyE = lme(CR.phyE ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.phyE)  
testInteractions(CR.phyE, pairwise="habitat.ordered")
plot(CR.phyE)
qqnorm(CR.phyE)
CR.phyE.I <- interactionMeans(CR.phyE) # effect plots
plot(CR.phyE.I, errorbar="ci95")

meanD.phyE = lme(phyE.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.phyE)  
testInteractions(meanD.phyE, pairwise="habitat.ordered")
plot(meanD.phyE)
qqnorm(meanD.phyE)
meanD.phyE.I <- interactionMeans(meanD.phyE) # effect plots
plot(meanD.phyE.I, errorbar="ci95")

Balance.phyE = lme(phyE.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.phyE)  
testInteractions(Balance.phyE, pairwise="habitat.ordered")
plot(Balance.phyE)
qqnorm(Balance.phyE)
Balance.phyE.I <- interactionMeans(Balance.phyE) # effect plots
plot(Balance.phyE.I, errorbar="ci95")

# Final figure (effect +/- standard error)

a <- ggplot(QE.phyE.I , aes(x= habitat.ordered, y=QE.phyE.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.phyE.I [,2]-QE.phyE.I [,3], ymax=QE.phyE.I [,2]+QE.phyE.I [,3]), width=.2) +
  geom_point(data=QE.phyE.I , mapping=aes(x=habitat.ordered, y=QE.phyE.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Phylogenetic diversity (QE), Erickson’s phylogeny", cex=16) +
  geom_text(aes(label= c("ac","b","ac","abc","bc")))

b <- ggplot(CR.phyE.I , aes(x= habitat.ordered, y=CR.phyE.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.phyE.I [,2]-CR.phyE.I [,3], ymax=CR.phyE.I [,2]+CR.phyE.I [,3]), width=.2) +
  geom_point(data=CR.phyE.I , mapping=aes(x=habitat.ordered, y=CR.phyE.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, Erickson’s phylogeny ", cex=16) +
  geom_text(aes(label= c("ab","a","bc","bc","c")))

c <- ggplot(meanD.phyE.I , aes(x= habitat.ordered, y=meanD.phyE.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.phyE.I [,2]-meanD.phyE.I [,3], ymax=meanD.phyE.I [,2]+meanD.phyE.I [,3]), width=.2) +
  geom_point(data=meanD.phyE.I , mapping=aes(x=habitat.ordered, y=meanD.phyE.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, Erickson’s phylogeny ", cex=16) +
  geom_text(aes(label= c("ab","a","b","ab","ab")))

d <- ggplot(Balance.phyE.I , aes(x= habitat.ordered, y=Balance.phyE.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.phyE.I [,2]-Balance.phyE.I [,3], ymax=Balance.phyE.I [,2]+Balance.phyE.I [,3]), width=.2) +
  geom_point(data=Balance.phyE.I , mapping=aes(x=habitat.ordered, y=Balance.phyE.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, Erickson’s phylogeny ", cex=16) +
  geom_text(aes(label= c("a","a","a","b","b")))

tiff(paste0(GoogleFigs,"/plot_FD_Erickson_phylogeny.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()




# Phylogenetic diversity Hackett

QE.phyH = lme(QE.phyH ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.phyH)
testInteractions(QE.phyH, pairwise="habitat.ordered")
plot(QE.phyH)
qqnorm(QE.phyH)
QE.phyH.I <- interactionMeans(QE.phyH) # effect plots
plot(QE.phyH.I, errorbar="ci95")

CR.phyH = lme(CR.phyH ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.phyH)  
testInteractions(CR.phyH, pairwise="habitat.ordered")
plot(CR.phyH)
qqnorm(CR.phyH)
CR.phyH.I <- interactionMeans(CR.phyH) # effect plots
plot(CR.phyH.I, errorbar="ci95")

meanD.phyH = lme(phyH.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.phyH)  
testInteractions(meanD.phyH, pairwise="habitat.ordered")
plot(meanD.phyH)
qqnorm(meanD.phyH)
meanD.phyH.I <- interactionMeans(meanD.phyH) # effect plots
plot(meanD.phyH.I, errorbar="ci95")

Balance.phyH = lme(phyH.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(Balance.phyH)  
testInteractions(Balance.phyH, pairwise="habitat.ordered")
plot(Balance.phyH)
qqnorm(Balance.phyH)
Balance.phyH.I <- interactionMeans(Balance.phyH) # effect plots
plot(Balance.phyH.I, errorbar="ci95")

# Final figure (effect +/- standard error)

a <- ggplot(QE.phyH.I , aes(x= habitat.ordered, y=QE.phyH.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.phyH.I [,2]-QE.phyH.I [,3], ymax=QE.phyH.I [,2]+QE.phyH.I [,3]), width=.2) +
  geom_point(data=QE.phyH.I , mapping=aes(x=habitat.ordered, y=QE.phyH.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Phylogenetic diversity (QE), Hackett’s phylogeny", cex=16) +
  geom_text(aes(label= c("a","b","ac","abc","b")))

b <- ggplot(CR.phyH.I , aes(x= habitat.ordered, y=CR.phyH.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.phyH.I [,2]-CR.phyH.I [,3], ymax=CR.phyH.I [,2]+CR.phyH.I [,3]), width=.2) +
  geom_point(data=CR.phyH.I , mapping=aes(x=habitat.ordered, y=CR.phyH.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, Hackett’s phylogeny ", cex=16) +
  geom_text(aes(label= c("a","b","ac","abc","c")))

c <- ggplot(meanD.phyH.I, aes(x= habitat.ordered, y=meanD.phyH.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.phyH.I [,2]-meanD.phyH.I [,3], ymax=meanD.phyH.I [,2]+meanD.phyH.I [,3]), width=.2) +
  geom_point(data=meanD.phyH.I , mapping=aes(x=habitat.ordered, y=meanD.phyH.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, Hackett’s phylogeny ", cex=16) +
  geom_text(aes(label= c("abc","b","c","abc","abc")))

d <- ggplot(Balance.phyH.I , aes(x= habitat.ordered, y=Balance.phyH.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.phyH.I [,2]-Balance.phyH.I [,3], ymax=Balance.phyH.I [,2]+Balance.phyH.I [,3]), width=.2) +
  geom_point(data=Balance.phyH.I , mapping=aes(x=habitat.ordered, y=Balance.phyH.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, Hackett’s phylogeny ", cex=16) +
  geom_text(aes(label= c("a","a","ab","bc","c")))

tiff(paste0(GoogleFigs,"/plot_FD_Hackett_phylogeny.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()
}





########### 2. Analyses for morphological diversity, only natives ##############
################################################################################

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities natives.txt"))
  
  # We restrict the analyses to studies where there is information inside and outside the city
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  
  
  # Three ways to code habitats
  
  # 1. All habitats separated
  # levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")
  
  # 2. All urban habitats separated, all non-urban habitat together
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  # 3. All urban and non-urban habitats pooled together
  # levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")
  
  
  
  ## Tests of the effect of urbanization on biodiversity*
  ################################################################
  
  # * Note that the best random structure has been previously evaluated with the method=REML
  # * We need to include confounds (city age, human density, coordinates)
  
  
  # Species richness
  
  spp.richness = lme(Species.richness ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(spp.richness)
  testInteractions(spp.richness, pairwise="habitat.ordered")
  plot(spp.richness)
  qqnorm(spp.richness)
  spp.richness.I <- interactionMeans(spp.richness) # effect plots
  plot(spp.richness.I, errorbar="ci95")
  
  
  # Simpson's index
  
  QE.taxonomy = lme(QE.taxonomy ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.taxonomy)
  testInteractions(QE.taxonomy, pairwise="habitat.ordered")
  plot(QE.taxonomy)
  qqnorm(QE.taxonomy)
  QE.taxonomy.I <- interactionMeans(QE.taxonomy) # effect plots
  plot(QE.taxonomy.I, errorbar="ci95")
  
  
  # Simpson*R/(R-1), an index of abundance evenness, independent of species richness. 
  
  Abundance.Eveness <- x$QE.taxonomy*x$Species.richness/(x$Species.richness-1)
  x <- cbind(x,Abundance.Eveness)
  
  AEveness = lme(Abundance.Eveness ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(AEveness)
  testInteractions(AEveness, pairwise="habitat.ordered")
  plot(AEveness)
  qqnorm(AEveness)
  AEveness.I <- interactionMeans(AEveness) # effect plots
  plot(AEveness.I, errorbar="ci95")
  
  
  
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(spp.richness.I, aes(x= habitat.ordered, y=spp.richness.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=spp.richness.I[,2]-spp.richness.I[,3], ymax=spp.richness.I[,2]+spp.richness.I[,3]), width=.2) +
    geom_point(data=spp.richness.I, mapping=aes(x=habitat.ordered, y=spp.richness.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Species richness", cex=16) +
    geom_text(aes(label= c("a","a","a","a","b")))
  
  b <- ggplot(QE.taxonomy.I, aes(x= habitat.ordered, y=QE.taxonomy.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.taxonomy.I[,2]-QE.taxonomy.I[,3], ymax=QE.taxonomy.I[,2]+QE.taxonomy.I[,3]), width=.2) +
    geom_point(data=QE.taxonomy.I, mapping=aes(x=habitat.ordered, y=QE.taxonomy.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Simpson's index", cex=16) +
    geom_text(aes(label= c("a","a","a","b","c")))
  
  c <- ggplot(AEveness.I, aes(x= habitat.ordered, y=AEveness.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=AEveness.I[,2]-AEveness.I[,3], ymax=AEveness.I[,2]+AEveness.I[,3]), width=.2) +
    geom_point(data=AEveness.I, mapping=aes(x=habitat.ordered, y=AEveness.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Abundance evenness index", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b")))
  
  tiff(paste0(GoogleFigs,"/plot_taxonomic_diversity_natives.tiff"), width = 9, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c, cols=1)
  dev.off()
  
  
  
  
  
  
  
  
  # 8 functional traits 
  
  QE.all.morph = lme(QE.all.morph ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.all.morph)
  testInteractions(QE.all.morph, pairwise="habitat.ordered")
  plot(QE.all.morph)
  qqnorm(QE.all.morph)
  QE.all.morph.I <- interactionMeans(QE.all.morph) # effect plots
  plot(QE.all.morph.I, errorbar="ci95")
  
  CR.all.morph = lme(CR.all.morph ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.all.morph)  
  testInteractions(CR.all.morph, pairwise="habitat.ordered")
  plot(CR.all.morph)
  qqnorm(CR.all.morph)
  CR.all.morph.I <- interactionMeans(CR.all.morph) # effect plots
  plot(CR.all.morph.I, errorbar="ci95")
  
  meanD.all.morph = lme(all.morph.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.all.morph)  
  testInteractions(meanD.all.morph, pairwise="habitat.ordered")
  plot(meanD.all.morph)
  qqnorm(meanD.all.morph)
  meanD.all.morph.I <- interactionMeans(meanD.all.morph) # effect plots
  plot(meanD.all.morph.I, errorbar="ci95")
  
  Balance.all.morph = lme(all.morph.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(Balance.all.morph)  
  testInteractions(Balance.all.morph, pairwise="habitat.ordered")
  plot(Balance.all.morph)
  qqnorm(Balance.all.morph)
  Balance.all.morph.I <- interactionMeans(Balance.all.morph) # effect plots
  plot(Balance.all.morph.I, errorbar="ci95")
  
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(QE.all.morph.I, aes(x= habitat.ordered, y=QE.all.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.all.morph.I[,2]-QE.all.morph.I[,3], ymax=QE.all.morph.I[,2]+QE.all.morph.I[,3]), width=.2) +
    geom_point(data=QE.all.morph.I, mapping=aes(x=habitat.ordered, y=QE.all.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional diversity (QE), all traits", cex=16) +
    geom_text(aes(label= c("a","ab","ab","b","b")))
  
  b <- ggplot(CR.all.morph.I, aes(x= habitat.ordered, y=CR.all.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.all.morph.I[,2]-CR.all.morph.I[,3], ymax=CR.all.morph.I[,2]+CR.all.morph.I[,3]), width=.2) +
    geom_point(data=CR.all.morph.I, mapping=aes(x=habitat.ordered, y=CR.all.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, all traits", cex=16) +
    geom_text(aes(label= c("ab","a","ab","ab","b")))
  
  c <- ggplot(meanD.all.morph.I, aes(x= habitat.ordered, y=meanD.all.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.all.morph.I[,2]-meanD.all.morph.I[,3], ymax=meanD.all.morph.I[,2]+meanD.all.morph.I[,3]), width=.2) +
    geom_point(data=meanD.all.morph.I, mapping=aes(x=habitat.ordered, y=meanD.all.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, all traits", cex=16) +
    geom_text(aes(label= c("a","b","ab","ab","a")))
  
  d <- ggplot(Balance.all.morph.I, aes(x= habitat.ordered, y=Balance.all.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.all.morph.I[,2]-Balance.all.morph.I[,3], ymax=Balance.all.morph.I[,2]+Balance.all.morph.I[,3]), width=.2) +
    geom_point(data=Balance.all.morph.I, mapping=aes(x=habitat.ordered, y=Balance.all.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, all traits", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_all_traits_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
  
  
  # All three PCAs (body size, size shape and locomotory shape)
  
  QE.PCA3 = lme(QE.PCA3 ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.PCA3)
  testInteractions(QE.PCA3, pairwise="habitat.ordered")
  plot(QE.PCA3)
  qqnorm(QE.PCA3)
  QE.PCA3.I <- interactionMeans(QE.PCA3) # effect plots
  plot(QE.PCA3.I, errorbar="ci95")
  
  CR.PCA3 = lme(CR.PCA3 ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.PCA3)  
  testInteractions(CR.PCA3, pairwise="habitat.ordered")
  plot(CR.PCA3)
  qqnorm(CR.PCA3)
  CR.PCA3.I <- interactionMeans(CR.PCA3) # effect plots
  plot(CR.PCA3.I, errorbar="ci95")
  
  meanD.PCA3 = lme(PCA3.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.PCA3)  
  testInteractions(meanD.PCA3, pairwise="habitat.ordered")
  plot(meanD.PCA3)
  qqnorm(meanD.PCA3)
  meanD.PCA3.I <- interactionMeans(meanD.PCA3) # effect plots
  plot(meanD.PCA3.I, errorbar="ci95")
  
  Balance.PCA3 = lme(PCA3.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(Balance.PCA3)  
  testInteractions(Balance.PCA3, pairwise="habitat.ordered")
  plot(Balance.PCA3)
  qqnorm(Balance.PCA3)
  Balance.PCA3.I <- interactionMeans(Balance.PCA3) # effect plots
  plot(Balance.PCA3.I, errorbar="ci95")
  
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(QE.PCA3.I, aes(x= habitat.ordered, y=QE.PCA3.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.PCA3.I[,2]-QE.PCA3.I[,3], ymax=QE.PCA3.I[,2]+QE.PCA3.I[,3]), width=.2) +
    geom_point(data=QE.PCA3.I, mapping=aes(x=habitat.ordered, y=QE.PCA3.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional diversity (QE), PCA3", cex=16) +
    geom_text(aes(label= c("a","b","ab","b","b")))
  
  b <- ggplot(CR.PCA3.I, aes(x= habitat.ordered, y=CR.PCA3.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.PCA3.I[,2]-CR.PCA3.I[,3], ymax=CR.PCA3.I[,2]+CR.PCA3.I[,3]), width=.2) +
    geom_point(data=CR.PCA3.I, mapping=aes(x=habitat.ordered, y=CR.PCA3.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, PCA3", cex=16) +
    geom_text(aes(label= c("ab","a","ab","ab","b")))
  
  c <- ggplot(meanD.PCA3.I, aes(x= habitat.ordered, y=meanD.PCA3.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.PCA3.I[,2]-meanD.PCA3.I[,3], ymax=meanD.PCA3.I[,2]+meanD.PCA3.I[,3]), width=.2) +
    geom_point(data=meanD.PCA3.I, mapping=aes(x=habitat.ordered, y=meanD.PCA3.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, PCA3", cex=16) +
    geom_text(aes(label= c("a","b","ab","ab","a")))
  
  d <- ggplot(Balance.PCA3.I, aes(x= habitat.ordered, y=Balance.PCA3.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.PCA3.I[,2]-Balance.PCA3.I[,3], ymax=Balance.PCA3.I[,2]+Balance.PCA3.I[,3]), width=.2) +
    geom_point(data=Balance.PCA3.I, mapping=aes(x=habitat.ordered, y=Balance.PCA3.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, PCA3", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_PCA3_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
  
  
  
  
  # beak shape
  
  QE.beak = lme(QE.beak ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.beak)
  testInteractions(QE.beak, pairwise="habitat.ordered")
  plot(QE.beak)
  qqnorm(QE.beak)
  QE.beak.I <- interactionMeans(QE.beak) # effect plots
  plot(QE.beak.I, errorbar="ci95")
  
  CR.beak = lme(CR.beak ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.beak)  
  testInteractions(CR.beak, pairwise="habitat.ordered")
  plot(CR.beak)
  qqnorm(CR.beak)
  CR.beak.I <- interactionMeans(CR.beak) # effect plots
  plot(CR.beak.I, errorbar="ci95")
  
  meanD.beak = lme(beak.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.beak)  
  testInteractions(meanD.beak, pairwise="habitat.ordered")
  plot(meanD.beak)
  qqnorm(meanD.beak)
  meanD.beak.I <- interactionMeans(meanD.beak) # effect plots
  plot(meanD.beak.I, errorbar="ci95")
  
  Balance.beak = lme(beak.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(Balance.beak)  
  testInteractions(Balance.beak, pairwise="habitat.ordered")
  plot(Balance.beak)
  qqnorm(Balance.beak)
  Balance.beak.I <- interactionMeans(Balance.beak) # effect plots
  plot(Balance.beak.I, errorbar="ci95")
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(QE.beak.I, aes(x= habitat.ordered, y=QE.beak.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.beak.I[,2]-QE.beak.I[,3], ymax=QE.beak.I[,2]+QE.beak.I[,3]), width=.2) +
    geom_point(data=QE.beak.I, mapping=aes(x=habitat.ordered, y=QE.beak.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional diversity (QE), beak", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  b <- ggplot(CR.beak.I, aes(x= habitat.ordered, y=CR.beak.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.beak.I[,2]-CR.beak.I[,3], ymax=CR.beak.I[,2]+CR.beak.I[,3]), width=.2) +
    geom_point(data=CR.beak.I, mapping=aes(x=habitat.ordered, y=CR.beak.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, beak", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b")))
  
  c <- ggplot(meanD.beak.I, aes(x= habitat.ordered, y=meanD.beak.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.beak.I[,2]-meanD.beak.I[,3], ymax=meanD.beak.I[,2]+meanD.beak.I[,3]), width=.2) +
    geom_point(data=meanD.beak.I, mapping=aes(x=habitat.ordered, y=meanD.beak.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, beak", cex=16) +
    geom_text(aes(label= c("a","ab","ab","ab","bc")))
  
  d <- ggplot(Balance.beak.I, aes(x= habitat.ordered, y=Balance.beak.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.beak.I[,2]-Balance.beak.I[,3], ymax=Balance.beak.I[,2]+Balance.beak.I[,3]), width=.2) +
    geom_point(data=Balance.beak.I, mapping=aes(x=habitat.ordered, y=Balance.beak.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, beak", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b")))  # check why the summary of the model() and the testInteractions() give different results
  
  tiff(paste0(GoogleFigs,"/plot_FD_beak_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
  
  
  # locomotory shape (tarsus vs tail shape)
  
  QE.locom = lme(QE.locom ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.locom)
  testInteractions(QE.locom, pairwise="habitat.ordered")
  plot(QE.locom)
  qqnorm(QE.locom)
  QE.locom.I <- interactionMeans(QE.locom) # effect plots
  plot(QE.locom.I, errorbar="ci95")
  
  CR.locom = lme(CR.locom ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.locom)  
  testInteractions(CR.locom, pairwise="habitat.ordered")
  plot(CR.locom)
  qqnorm(CR.locom)
  CR.locom.I <- interactionMeans(CR.locom) # effect plots
  plot(CR.locom.I, errorbar="ci95")
  
  meanD.locom = lme(locom.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.locom)  
  testInteractions(meanD.locom, pairwise="habitat.ordered")
  plot(meanD.locom)
  qqnorm(meanD.locom)
  meanD.locom.I <- interactionMeans(meanD.locom) # effect plots
  plot(meanD.locom.I, errorbar="ci95")
  
  Balance.locom = lme(locom.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(Balance.locom)  
  testInteractions(Balance.locom, pairwise="habitat.ordered")
  plot(Balance.locom)
  qqnorm(Balance.locom)
  Balance.locom.I <- interactionMeans(Balance.locom) # effect plots
  plot(Balance.locom.I, errorbar="ci95")
  
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(QE.locom.I , aes(x= habitat.ordered, y=QE.locom.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.locom.I [,2]-QE.locom.I [,3], ymax=QE.locom.I [,2]+QE.locom.I [,3]), width=.2) +
    geom_point(data=QE.locom.I , mapping=aes(x=habitat.ordered, y=QE.locom.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional diversity (QE), locomotory system", cex=16) +
    geom_text(aes(label= c("a","ab","ab","ab","b")))
  
  b <- ggplot(CR.locom.I , aes(x= habitat.ordered, y=CR.locom.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.locom.I [,2]-CR.locom.I [,3], ymax=CR.locom.I [,2]+CR.locom.I [,3]), width=.2) +
    geom_point(data=CR.locom.I , mapping=aes(x=habitat.ordered, y=CR.locom.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, locomotory system", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  c <- ggplot(meanD.locom.I , aes(x= habitat.ordered, y=meanD.locom.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.locom.I [,2]-meanD.locom.I [,3], ymax=meanD.locom.I [,2]+meanD.locom.I [,3]), width=.2) +
    geom_point(data=meanD.locom.I , mapping=aes(x=habitat.ordered, y=meanD.locom.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, locomotory system", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  d <- ggplot(Balance.locom.I , aes(x= habitat.ordered, y=Balance.locom.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.locom.I [,2]-Balance.locom.I [,3], ymax=Balance.locom.I [,2]+Balance.locom.I [,3]), width=.2) +
    geom_point(data=Balance.locom.I , mapping=aes(x=habitat.ordered, y=Balance.locom.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, locomotory system", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_Locomotory_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
  
  
  
  # size
  
  QE.size = lme(QE.size ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.size)
  testInteractions(QE.size, pairwise="habitat.ordered")
  plot(QE.size)
  qqnorm(QE.size)
  QE.size.I <- interactionMeans(QE.size) # effect plots
  plot(QE.size.I, errorbar="ci95")
  
  CR.size = lme(CR.size ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.size)  
  testInteractions(CR.size, pairwise="habitat.ordered")
  plot(CR.size)
  qqnorm(CR.size)
  CR.size.I <- interactionMeans(CR.size) # effect plots
  plot(CR.size.I, errorbar="ci95")
  
  meanD.size = lme(size.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.size)  
  testInteractions(meanD.size, pairwise="habitat.ordered")
  plot(meanD.size)
  qqnorm(meanD.size)
  meanD.size.I <- interactionMeans(meanD.size) # effect plots
  plot(meanD.size.I, errorbar="ci95")
  
  Balance.size = lme(size.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(Balance.size)  
  testInteractions(Balance.size, pairwise="habitat.ordered")
  plot(Balance.size)
  qqnorm(Balance.size)
  Balance.size.I <- interactionMeans(Balance.size) # effect plots
  plot(Balance.size.I, errorbar="ci95")
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(QE.size.I , aes(x= habitat.ordered, y=QE.size.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.size.I [,2]-QE.size.I [,3], ymax=QE.size.I [,2]+QE.size.I [,3]), width=.2) +
    geom_point(data=QE.size.I , mapping=aes(x=habitat.ordered, y=QE.size.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional diversity (QE), body size", cex=16) +
    geom_text(aes(label= c("a","ab","ab","b","b")))
  
  b <- ggplot(CR.size.I , aes(x= habitat.ordered, y=CR.size.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.size.I [,2]-CR.size.I [,3], ymax=CR.size.I [,2]+CR.size.I [,3]), width=.2) +
    geom_point(data=CR.size.I , mapping=aes(x=habitat.ordered, y=CR.size.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, body size", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  c <- ggplot(meanD.size.I , aes(x= habitat.ordered, y=meanD.size.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.size.I [,2]-meanD.size.I [,3], ymax=meanD.size.I [,2]+meanD.size.I [,3]), width=.2) +
    geom_point(data=meanD.size.I , mapping=aes(x=habitat.ordered, y=meanD.size.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, body size", cex=16) +
    geom_text(aes(label= c("a","b","ab","ab","a")))
  
  d <- ggplot(Balance.size.I , aes(x= habitat.ordered, y=Balance.size.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.size.I [,2]-Balance.size.I [,3], ymax=Balance.size.I [,2]+Balance.size.I [,3]), width=.2) +
    geom_point(data=Balance.size.I , mapping=aes(x=habitat.ordered, y=Balance.size.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, body size", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_Body_size_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
  
  
  
  # winghand
  
  QE.winghand = lme(QE.winghand ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.winghand)
  testInteractions(QE.winghand, pairwise="habitat.ordered")
  plot(QE.winghand)
  qqnorm(QE.winghand)
  QE.winghand.I <- interactionMeans(QE.winghand) # effect plots
  plot(QE.winghand.I, errorbar="ci95")
  
  CR.winghand = lme(CR.winghand ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.winghand)  
  testInteractions(CR.winghand, pairwise="habitat.ordered")
  plot(CR.winghand)
  qqnorm(CR.winghand)
  CR.winghand.I <- interactionMeans(CR.winghand) # effect plots
  plot(CR.winghand.I, errorbar="ci95")
  
  meanD.winghand = lme(winghand.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.winghand)  
  testInteractions(meanD.winghand, pairwise="habitat.ordered")
  plot(meanD.winghand)
  qqnorm(meanD.winghand)
  meanD.winghand.I <- interactionMeans(meanD.winghand) # effect plots
  plot(meanD.winghand.I, errorbar="ci95")
  
  Balance.winghand = lme(winghand.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(Balance.winghand)  
  testInteractions(Balance.winghand, pairwise="habitat.ordered")
  plot(Balance.winghand)
  qqnorm(Balance.winghand)
  Balance.winghand.I <- interactionMeans(Balance.winghand) # effect plots
  plot(Balance.winghand.I, errorbar="ci95")
  
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(QE.winghand.I , aes(x= habitat.ordered, y=QE.winghand.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.winghand.I [,2]-QE.winghand.I [,3], ymax=QE.winghand.I [,2]+QE.winghand.I [,3]), width=.2) +
    geom_point(data=QE.winghand.I , mapping=aes(x=habitat.ordered, y=QE.winghand.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional diversity (QE), wing hand", cex=16) +
    geom_text(aes(label= c("a","ab","ab","ab","b")))
  
  b <- ggplot(CR.winghand.I , aes(x= habitat.ordered, y=CR.winghand.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.winghand.I [,2]-CR.winghand.I [,3], ymax=CR.winghand.I [,2]+CR.winghand.I [,3]), width=.2) +
    geom_point(data=CR.winghand.I , mapping=aes(x=habitat.ordered, y=CR.winghand.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, wing hand ", cex=16) +
    geom_text(aes(label= c("ab","ab","a","ab","b")))
  
  c <- ggplot(meanD.winghand.I , aes(x= habitat.ordered, y=meanD.winghand.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.winghand.I [,2]-meanD.winghand.I [,3], ymax=meanD.winghand.I [,2]+meanD.winghand.I [,3]), width=.2) +
    geom_point(data=meanD.winghand.I , mapping=aes(x=habitat.ordered, y=meanD.winghand.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, wing hand ", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  d <- ggplot(Balance.winghand.I , aes(x= habitat.ordered, y=Balance.winghand.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.winghand.I [,2]-Balance.winghand.I [,3], ymax=Balance.winghand.I [,2]+Balance.winghand.I [,3]), width=.2) +
    geom_point(data=Balance.winghand.I , mapping=aes(x=habitat.ordered, y=Balance.winghand.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, wing hand ", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_wing_hand_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
  
  
  
  # Phylogenetic diversity Ericksson
  
  QE.phyE = lme(QE.phyE ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.phyE)
  testInteractions(QE.phyE, pairwise="habitat.ordered")
  plot(QE.phyE)
  qqnorm(QE.phyE)
  QE.phyE.I <- interactionMeans(QE.phyE) # effect plots
  plot(QE.phyE.I, errorbar="ci95")
  
  CR.phyE = lme(CR.phyE ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.phyE)  
  testInteractions(CR.phyE, pairwise="habitat.ordered")
  plot(CR.phyE)
  qqnorm(CR.phyE)
  CR.phyE.I <- interactionMeans(CR.phyE) # effect plots
  plot(CR.phyE.I, errorbar="ci95")
  
  meanD.phyE = lme(phyE.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.phyE)  
  testInteractions(meanD.phyE, pairwise="habitat.ordered")
  plot(meanD.phyE)
  qqnorm(meanD.phyE)
  meanD.phyE.I <- interactionMeans(meanD.phyE) # effect plots
  plot(meanD.phyE.I, errorbar="ci95")
  
  Balance.phyE = lme(phyE.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(Balance.phyE)  
  testInteractions(Balance.phyE, pairwise="habitat.ordered")
  plot(Balance.phyE)
  qqnorm(Balance.phyE)
  Balance.phyE.I <- interactionMeans(Balance.phyE) # effect plots
  plot(Balance.phyE.I, errorbar="ci95")
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(QE.phyE.I , aes(x= habitat.ordered, y=QE.phyE.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.phyE.I [,2]-QE.phyE.I [,3], ymax=QE.phyE.I [,2]+QE.phyE.I [,3]), width=.2) +
    geom_point(data=QE.phyE.I , mapping=aes(x=habitat.ordered, y=QE.phyE.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Phylogenetic diversity (QE), Erickson’s phylogeny", cex=16) +
    geom_text(aes(label= c("a","b","a","ab","b")))
  
  b <- ggplot(CR.phyE.I , aes(x= habitat.ordered, y=CR.phyE.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.phyE.I [,2]-CR.phyE.I [,3], ymax=CR.phyE.I [,2]+CR.phyE.I [,3]), width=.2) +
    geom_point(data=CR.phyE.I , mapping=aes(x=habitat.ordered, y=CR.phyE.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, Erickson’s phylogeny ", cex=16) +
    geom_text(aes(label= c("ab","b","a","a","a")))
  
  c <- ggplot(meanD.phyE.I , aes(x= habitat.ordered, y=meanD.phyE.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.phyE.I [,2]-meanD.phyE.I [,3], ymax=meanD.phyE.I [,2]+meanD.phyE.I [,3]), width=.2) +
    geom_point(data=meanD.phyE.I , mapping=aes(x=habitat.ordered, y=meanD.phyE.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, Erickson’s phylogeny ", cex=16) +
    geom_text(aes(label= c("ab","a","b","ab","ab")))
  
  d <- ggplot(Balance.phyE.I , aes(x= habitat.ordered, y=Balance.phyE.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.phyE.I [,2]-Balance.phyE.I [,3], ymax=Balance.phyE.I [,2]+Balance.phyE.I [,3]), width=.2) +
    geom_point(data=Balance.phyE.I , mapping=aes(x=habitat.ordered, y=Balance.phyE.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, Erickson’s phylogeny ", cex=16) +
    geom_text(aes(label= c("a","a","ab","b","b")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_Erickson_phylogeny_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
  
  
  
  
  # Phylogenetic diversity Hackett
  
  QE.phyH = lme(QE.phyH ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.phyH)
  testInteractions(QE.phyH, pairwise="habitat.ordered")
  plot(QE.phyH)
  qqnorm(QE.phyH)
  QE.phyH.I <- interactionMeans(QE.phyH) # effect plots
  plot(QE.phyH.I, errorbar="ci95")
  
  CR.phyH = lme(CR.phyH ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(CR.phyH)  
  testInteractions(CR.phyH, pairwise="habitat.ordered")
  plot(CR.phyH)
  qqnorm(CR.phyH)
  CR.phyH.I <- interactionMeans(CR.phyH) # effect plots
  plot(CR.phyH.I, errorbar="ci95")
  
  meanD.phyH = lme(phyH.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(meanD.phyH)  
  testInteractions(meanD.phyH, pairwise="habitat.ordered")
  plot(meanD.phyH)
  qqnorm(meanD.phyH)
  meanD.phyH.I <- interactionMeans(meanD.phyH) # effect plots
  plot(meanD.phyH.I, errorbar="ci95")
  
  Balance.phyH = lme(phyH.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(Balance.phyH)  
  testInteractions(Balance.phyH, pairwise="habitat.ordered")
  plot(Balance.phyH)
  qqnorm(Balance.phyH)
  Balance.phyH.I <- interactionMeans(Balance.phyH) # effect plots
  plot(Balance.phyH.I, errorbar="ci95")
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(QE.phyH.I , aes(x= habitat.ordered, y=QE.phyH.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.phyH.I [,2]-QE.phyH.I [,3], ymax=QE.phyH.I [,2]+QE.phyH.I [,3]), width=.2) +
    geom_point(data=QE.phyH.I , mapping=aes(x=habitat.ordered, y=QE.phyH.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Phylogenetic diversity (QE), Hackett’s phylogeny", cex=16) +
    geom_text(aes(label= c("a","b","ac","abc","b")))
  
  b <- ggplot(CR.phyH.I , aes(x= habitat.ordered, y=CR.phyH.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.phyH.I [,2]-CR.phyH.I [,3], ymax=CR.phyH.I [,2]+CR.phyH.I [,3]), width=.2) +
    geom_point(data=CR.phyH.I , mapping=aes(x=habitat.ordered, y=CR.phyH.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, Hackett’s phylogeny ", cex=16) +
    geom_text(aes(label= c("a","a","b","b","b")))  # check significances
  
  c <- ggplot(meanD.phyH.I, aes(x= habitat.ordered, y=meanD.phyH.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.phyH.I [,2]-meanD.phyH.I [,3], ymax=meanD.phyH.I [,2]+meanD.phyH.I [,3]), width=.2) +
    geom_point(data=meanD.phyH.I , mapping=aes(x=habitat.ordered, y=meanD.phyH.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, Hackett’s phylogeny ", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  d <- ggplot(Balance.phyH.I , aes(x= habitat.ordered, y=Balance.phyH.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.phyH.I [,2]-Balance.phyH.I [,3], ymax=Balance.phyH.I [,2]+Balance.phyH.I [,3]), width=.2) +
    geom_point(data=Balance.phyH.I , mapping=aes(x=habitat.ordered, y=Balance.phyH.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, Hackett’s phylogeny ", cex=16) +
    geom_text(aes(label= c("a","a","ab","b","b")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_Hackett_phylogeny_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
}






########### 3. Analyses for morphological diversity, ocurrence data ##############
##################################################################################


{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities ocurrences.txt"))
  
# We restrict the analyses to studies where there is information inside and outside the city
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  
  # Three ways to code habitats
  
  # 1. All habitats separated
  # levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")
  
  # 2. All urban habitats separated, all non-urban habitat together
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  # 3. All urban and non-urban habitats pooled together
  # levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")
  
  
  
  ## Tests of the effect of urbanization on biodiversity*
  ################################################################
  
  # * Note that the best random structure has been previously evaluated with the method=REML
  # * We need to include confounds (city age, human density, coordinates)
  
  
  
  # Species richness
  
  spp.richness = lme(Species.richness ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(spp.richness)
  testInteractions(spp.richness, pairwise="habitat.ordered")
  plot(spp.richness)
  qqnorm(spp.richness)
  spp.richness.I <- interactionMeans(spp.richness) # effect plots
  plot(spp.richness.I, errorbar="ci95")
  
  
  
  # 8 functional traits 
  
  QE.all.morph = lme(QE.all.morph ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.all.morph)
  testInteractions(QE.all.morph, pairwise="habitat.ordered")
  plot(QE.all.morph)
  qqnorm(QE.all.morph)
  QE.all.morph.I <- interactionMeans(QE.all.morph) # effect plots
  plot(QE.all.morph.I, errorbar="ci95")
  
  
  
  a <- ggplot(spp.richness.I, aes(x= habitat.ordered, y=spp.richness.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=spp.richness.I[,2]-spp.richness.I[,3], ymax=spp.richness.I[,2]+spp.richness.I[,3]), width=.2) +
    geom_point(data=spp.richness.I, mapping=aes(x=habitat.ordered, y=spp.richness.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Species richness", cex=16) +
    geom_text(aes(label= c("a","a","a","a","b")))
  
  b <- ggplot(QE.all.morph.I, aes(x= habitat.ordered, y=QE.all.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.all.morph.I[,2]-QE.all.morph.I[,3], ymax=QE.all.morph.I[,2]+QE.all.morph.I[,3]), width=.2) +
    geom_point(data=QE.all.morph.I, mapping=aes(x=habitat.ordered, y=QE.all.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional diversity (QE), all traits", cex=16) +
    geom_text(aes(label= c("a","b","ab","ab","ab")))
  
  tiff(paste0(GoogleFigs,"/plot_Biodiversity_loss_Occurrences.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b)
  dev.off()
  }
  


########### 4. Analyses for diet diversity ##############
##################################################################

{
## Import biodiversity metrics for communities

x<-read.table(paste0(workingData,"/Diet diversity metrics for communities.txt"))

# We restrict the analyses to studies where there is information inside and outside the city
x <- subset(x, used.urban.nonurban=="yes")
x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)


# Three ways to code habitats

# 1. All habitats separated
# levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")

# 2. All urban habitats separated, all non-urban habitat together
levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
x <- cbind(x,habitat.ordered)

# 3. All urban and non-urban habitats pooled together
# levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")



## Tests of the effect of urbanization on biodiversity*
################################################################

# * Note that the best random structure has been previously evaluated with the method=REML
# * We need to include confounds (city age, human density, coordinates)

QE.diet = lme(QE.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.diet)
testInteractions(QE.diet, pairwise="habitat.ordered")
plot(QE.diet)
qqnorm(QE.diet)
QE.diet.I <- interactionMeans(QE.diet) # effect plots
plot(QE.diet.I, errorbar="ci95")

CR.diet = lme(CR.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.diet)  
testInteractions(CR.diet, pairwise="habitat.ordered")
plot(CR.diet)
qqnorm(CR.diet)
CR.diet.I <- interactionMeans(CR.diet) # effect plots
plot(CR.diet.I, errorbar="ci95")

meanD.diet = lme(diet.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.diet)  
testInteractions(meanD.diet, pairwise="habitat.ordered")
plot(meanD.diet)
qqnorm(meanD.diet)
meanD.diet.I <- interactionMeans(meanD.diet) # effect plots
plot(meanD.diet.I, errorbar="ci95")

Balance.diet = lme(diet.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
summary(Balance.diet)  
testInteractions(Balance.diet, pairwise="habitat.ordered")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
plot(Balance.diet)
qqnorm(Balance.diet)
Balance.diet.I <- interactionMeans(Balance.diet) # effect plots
plot(Balance.diet.I, errorbar="ci95")

# Final figure (effect +/- standard error)

a <- ggplot(QE.diet.I , aes(x= habitat.ordered, y=QE.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.diet.I [,2]-QE.diet.I [,3], ymax=QE.diet.I [,2]+QE.diet.I [,3]), width=.2) +
  geom_point(data=QE.diet.I , mapping=aes(x=habitat.ordered, y=QE.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Phylogenetic diversity (QE), diet", cex=16) +
  geom_text(aes(label= c("a","a","a","a","b")))

b <- ggplot(CR.diet.I , aes(x= habitat.ordered, y=CR.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.diet.I [,2]-CR.diet.I [,3], ymax=CR.diet.I [,2]+CR.diet.I [,3]), width=.2) +
  geom_point(data=CR.diet.I , mapping=aes(x=habitat.ordered, y=CR.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, diet ", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

c <- ggplot(meanD.diet.I, aes(x= habitat.ordered, y=meanD.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.diet.I [,2]-meanD.diet.I [,3], ymax=meanD.diet.I [,2]+meanD.diet.I [,3]), width=.2) +
  geom_point(data=meanD.diet.I , mapping=aes(x=habitat.ordered, y=meanD.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, diet ", cex=16) +
  geom_text(aes(label= c("a","ab","ab","ab","b")))

d <- ggplot(Balance.diet.I , aes(x= habitat.ordered, y=Balance.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.diet.I [,2]-Balance.diet.I [,3], ymax=Balance.diet.I [,2]+Balance.diet.I [,3]), width=.2) +
  geom_point(data=Balance.diet.I , mapping=aes(x=habitat.ordered, y=Balance.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, diet ", cex=16) +
  geom_text(aes(label= c("ab","ab","ab","a","b")))

tiff(paste0(GoogleFigs,"/plot_FD_diet.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()

}


########### 5. Analyses for diet diversity natives ###############
##################################################################

{
## Import biodiversity metrics for communities

x<-read.table(paste0(workingData,"/Diet diversity metrics for communities natives.txt"))

# We restrict the analyses to studies where there is information inside and outside the city
x <- subset(x, used.urban.nonurban=="yes")
x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)


# Three ways to code habitats

# 1. All habitats separated
# levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")

# 2. All urban habitats separated, all non-urban habitat together
levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
x <- cbind(x,habitat.ordered)

# 3. All urban and non-urban habitats pooled together
# levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")



## Tests of the effect of urbanization on biodiversity*
################################################################

# * Note that the best random structure has been previously evaluated with the method=REML
# * We need to include confounds (city age, human density, coordinates)

QE.diet = lme(QE.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.diet)
testInteractions(QE.diet, pairwise="habitat.ordered")
plot(QE.diet)
qqnorm(QE.diet)
QE.diet.I <- interactionMeans(QE.diet) # effect plots
plot(QE.diet.I, errorbar="ci95")

CR.diet = lme(CR.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.diet)  
testInteractions(CR.diet, pairwise="habitat.ordered")
plot(CR.diet)
qqnorm(CR.diet)
CR.diet.I <- interactionMeans(CR.diet) # effect plots
plot(CR.diet.I, errorbar="ci95")

meanD.diet = lme(diet.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.diet)  
testInteractions(meanD.diet, pairwise="habitat.ordered")
plot(meanD.diet)
qqnorm(meanD.diet)
meanD.diet.I <- interactionMeans(meanD.diet) # effect plots
plot(meanD.diet.I, errorbar="ci95")

Balance.diet = lme(diet.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
summary(Balance.diet)  
testInteractions(Balance.diet, pairwise="habitat.ordered")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
plot(Balance.diet)
qqnorm(Balance.diet)
Balance.diet.I <- interactionMeans(Balance.diet) # effect plots
plot(Balance.diet.I, errorbar="ci95")

# Final figure (effect +/- standard error)

a <- ggplot(QE.diet.I , aes(x= habitat.ordered, y=QE.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.diet.I [,2]-QE.diet.I [,3], ymax=QE.diet.I [,2]+QE.diet.I [,3]), width=.2) +
  geom_point(data=QE.diet.I , mapping=aes(x=habitat.ordered, y=QE.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Phylogenetic diversity (QE), diet", cex=16) +
  geom_text(aes(label= c("a","a","a","a","b")))

b <- ggplot(CR.diet.I , aes(x= habitat.ordered, y=CR.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.diet.I [,2]-CR.diet.I [,3], ymax=CR.diet.I [,2]+CR.diet.I [,3]), width=.2) +
  geom_point(data=CR.diet.I , mapping=aes(x=habitat.ordered, y=CR.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, diet ", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

c <- ggplot(meanD.diet.I, aes(x= habitat.ordered, y=meanD.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.diet.I [,2]-meanD.diet.I [,3], ymax=meanD.diet.I [,2]+meanD.diet.I [,3]), width=.2) +
  geom_point(data=meanD.diet.I , mapping=aes(x=habitat.ordered, y=meanD.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, diet ", cex=16) +
  geom_text(aes(label= c("a","ab","ab","ab","b")))

d <- ggplot(Balance.diet.I , aes(x= habitat.ordered, y=Balance.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.diet.I [,2]-Balance.diet.I [,3], ymax=Balance.diet.I [,2]+Balance.diet.I [,3]), width=.2) +
  geom_point(data=Balance.diet.I , mapping=aes(x=habitat.ordered, y=Balance.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, diet ", cex=16) +
  geom_text(aes(label= c("ab","ab","ab","a","b")))

tiff(paste0(GoogleFigs,"/plot_FD_diet_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()
}



########### 6. Analyses for diet diversity ocurrences ###############
######################################################################

{
## Import biodiversity metrics for communities

x<-read.table(paste0(workingData,"/Diet diversity metrics for communities ocurrences.txt"))

# We restrict the analyses to studies where there is information inside and outside the city
x <- subset(x, used.urban.nonurban=="yes")
x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)


# Three ways to code habitats

# 1. All habitats separated
# levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")

# 2. All urban habitats separated, all non-urban habitat together
levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
x <- cbind(x,habitat.ordered)

# 3. All urban and non-urban habitats pooled together
# levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")



## Tests of the effect of urbanization on biodiversity*
################################################################

# * Note that the best random structure has been previously evaluated with the method=REML
# * We need to include confounds (city age, human density, coordinates)

QE.diet = lme(QE.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.diet)
testInteractions(QE.diet, pairwise="habitat.ordered")
plot(QE.diet)
qqnorm(QE.diet)
QE.diet.I <- interactionMeans(QE.diet) # effect plots
plot(QE.diet.I, errorbar="ci95")

CR.diet = lme(CR.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.diet)  
testInteractions(CR.diet, pairwise="habitat.ordered")
plot(CR.diet)
qqnorm(CR.diet)
CR.diet.I <- interactionMeans(CR.diet) # effect plots
plot(CR.diet.I, errorbar="ci95")

meanD.diet = lme(diet.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.diet)  
testInteractions(meanD.diet, pairwise="habitat.ordered")
plot(meanD.diet)
qqnorm(meanD.diet)
meanD.diet.I <- interactionMeans(meanD.diet) # effect plots
plot(meanD.diet.I, errorbar="ci95")

Balance.diet = lme(diet.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
summary(Balance.diet)  
testInteractions(Balance.diet, pairwise="habitat.ordered")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
plot(Balance.diet)
qqnorm(Balance.diet)
Balance.diet.I <- interactionMeans(Balance.diet) # effect plots
plot(Balance.diet.I, errorbar="ci95")

# Final figure (effect +/- standard error)

a <- ggplot(QE.diet.I , aes(x= habitat.ordered, y=QE.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.diet.I [,2]-QE.diet.I [,3], ymax=QE.diet.I [,2]+QE.diet.I [,3]), width=.2) +
  geom_point(data=QE.diet.I , mapping=aes(x=habitat.ordered, y=QE.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Phylogenetic diversity (QE), diet", cex=16) +
  geom_text(aes(label= c("a","a","a","a","b")))

b <- ggplot(CR.diet.I , aes(x= habitat.ordered, y=CR.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.diet.I [,2]-CR.diet.I [,3], ymax=CR.diet.I [,2]+CR.diet.I [,3]), width=.2) +
  geom_point(data=CR.diet.I , mapping=aes(x=habitat.ordered, y=CR.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, diet ", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

c <- ggplot(meanD.diet.I, aes(x= habitat.ordered, y=meanD.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.diet.I [,2]-meanD.diet.I [,3], ymax=meanD.diet.I [,2]+meanD.diet.I [,3]), width=.2) +
  geom_point(data=meanD.diet.I , mapping=aes(x=habitat.ordered, y=meanD.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, diet ", cex=16) +
  geom_text(aes(label= c("a","ab","ab","ab","b")))

d <- ggplot(Balance.diet.I , aes(x= habitat.ordered, y=Balance.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.diet.I [,2]-Balance.diet.I [,3], ymax=Balance.diet.I [,2]+Balance.diet.I [,3]), width=.2) +
  geom_point(data=Balance.diet.I , mapping=aes(x=habitat.ordered, y=Balance.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, diet ", cex=16) +
  geom_text(aes(label= c("ab","ab","ab","a","b")))

tiff(paste0(GoogleFigs,"/plot_FD_diet_occurrences.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()

}


########### 7. Analyses for diet diversity ocurrence natives ###############
#############################################################################

{
## Import biodiversity metrics for communities

x<-read.table(paste0(workingData,"/Diet diversity metrics for communities ocurrences natives.txt"))

# We restrict the analyses to studies where there is information inside and outside the city
x <- subset(x, used.urban.nonurban=="yes")
x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)


# Three ways to code habitats

# 1. All habitats separated
# levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")

# 2. All urban habitats separated, all non-urban habitat together
levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
x <- cbind(x,habitat.ordered)

# 3. All urban and non-urban habitats pooled together
# levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")



## Tests of the effect of urbanization on biodiversity*
################################################################

# * Note that the best random structure has been previously evaluated with the method=REML
# * We need to include confounds (city age, human density, coordinates)

QE.diet = lme(QE.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(QE.diet)
testInteractions(QE.diet, pairwise="habitat.ordered")
plot(QE.diet)
qqnorm(QE.diet)
QE.diet.I <- interactionMeans(QE.diet) # effect plots
plot(QE.diet.I, errorbar="ci95")

CR.diet = lme(CR.diet ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(CR.diet)  
testInteractions(CR.diet, pairwise="habitat.ordered")
plot(CR.diet)
qqnorm(CR.diet)
CR.diet.I <- interactionMeans(CR.diet) # effect plots
plot(CR.diet.I, errorbar="ci95")

meanD.diet = lme(diet.meanD ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
summary(meanD.diet)  
testInteractions(meanD.diet, pairwise="habitat.ordered")
plot(meanD.diet)
qqnorm(meanD.diet)
meanD.diet.I <- interactionMeans(meanD.diet) # effect plots
plot(meanD.diet.I, errorbar="ci95")

Balance.diet = lme(diet.Balance ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
summary(Balance.diet)  
testInteractions(Balance.diet, pairwise="habitat.ordered")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
plot(Balance.diet)
qqnorm(Balance.diet)
Balance.diet.I <- interactionMeans(Balance.diet) # effect plots
plot(Balance.diet.I, errorbar="ci95")

# Final figure (effect +/- standard error)

a <- ggplot(QE.diet.I , aes(x= habitat.ordered, y=QE.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=QE.diet.I [,2]-QE.diet.I [,3], ymax=QE.diet.I [,2]+QE.diet.I [,3]), width=.2) +
  geom_point(data=QE.diet.I , mapping=aes(x=habitat.ordered, y=QE.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional diversity (QE), diet", cex=16) +
  geom_text(aes(label= c("a","a","a","a","b")))

b <- ggplot(CR.diet.I , aes(x= habitat.ordered, y=CR.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.diet.I [,2]-CR.diet.I [,3], ymax=CR.diet.I [,2]+CR.diet.I [,3]), width=.2) +
  geom_point(data=CR.diet.I , mapping=aes(x=habitat.ordered, y=CR.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Functional redundancies, diet ", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

c <- ggplot(meanD.diet.I, aes(x= habitat.ordered, y=meanD.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.diet.I [,2]-meanD.diet.I [,3], ymax=meanD.diet.I [,2]+meanD.diet.I [,3]), width=.2) +
  geom_point(data=meanD.diet.I , mapping=aes(x=habitat.ordered, y=meanD.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Average species dissimilarity, diet ", cex=16) +
  geom_text(aes(label= c("a","ab","ab","ab","b")))

d <- ggplot(Balance.diet.I , aes(x= habitat.ordered, y=Balance.diet.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.diet.I [,2]-Balance.diet.I [,3], ymax=Balance.diet.I [,2]+Balance.diet.I [,3]), width=.2) +
  geom_point(data=Balance.diet.I , mapping=aes(x=habitat.ordered, y=Balance.diet.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance component, diet ", cex=16) +
  geom_text(aes(label= c("ab","ab","ab","a","b")))

tiff(paste0(GoogleFigs,"/plot_FD_diet_ocurrence_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()
}




########### 8. Analyses for morphological diversity natives, ocurrence data ######
##################################################################################

### need to estimate input file

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities ocurrences natives.txt"))
  
  # We restrict the analyses to studies where there is information inside and outside the city
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  
  # Three ways to code habitats
  
  # 1. All habitats separated
  # levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")
  
  # 2. All urban habitats separated, all non-urban habitat together
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  # 3. All urban and non-urban habitats pooled together
  # levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")
  
  
  
  ## Tests of the effect of urbanization on biodiversity*
  ################################################################
  
  # * Note that the best random structure has been previously evaluated with the method=REML
  # * We need to include confounds (city age, human density, coordinates)
  
  
  
  # Species richness
  
  spp.richness = lme(Species.richness ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(spp.richness)
  testInteractions(spp.richness, pairwise="habitat.ordered")
  plot(spp.richness)
  qqnorm(spp.richness)
  spp.richness.I <- interactionMeans(spp.richness) # effect plots
  plot(spp.richness.I, errorbar="ci95")
  
  
  
  # 8 functional traits 
  
  QE.all.morph = lme(QE.all.morph ~ habitat.ordered, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.all.morph)
  testInteractions(QE.all.morph, pairwise="habitat.ordered")
  plot(QE.all.morph)
  qqnorm(QE.all.morph)
  QE.all.morph.I <- interactionMeans(QE.all.morph) # effect plots
  plot(QE.all.morph.I, errorbar="ci95")
  
  
  
  a <- ggplot(spp.richness.I, aes(x= habitat.ordered, y=spp.richness.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=spp.richness.I[,2]-spp.richness.I[,3], ymax=spp.richness.I[,2]+spp.richness.I[,3]), width=.2) +
    geom_point(data=spp.richness.I, mapping=aes(x=habitat.ordered, y=spp.richness.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Species richness", cex=16) +
    geom_text(aes(label= c("a","a","a","a","b")))
  
  b <- ggplot(QE.all.morph.I, aes(x= habitat.ordered, y=QE.all.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.all.morph.I[,2]-QE.all.morph.I[,3], ymax=QE.all.morph.I[,2]+QE.all.morph.I[,3]), width=.2) +
    geom_point(data=QE.all.morph.I, mapping=aes(x=habitat.ordered, y=QE.all.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional diversity (QE), all traits", cex=16) +
    geom_text(aes(label= c("a","b","ab","ab","ab")))
  
  tiff(paste0(GoogleFigs,"/plot_Biodiversity_loss_Occurrences_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b)
  dev.off()
}




############# 9. Analyses for morphology-diet diversity, only natives ######
############################################################################

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Morphology-Diet diversity metrics for communities natives.txt"))
  
  # We restrict the analyses to studies where there is information inside and outside the city
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  
  # Three ways to code habitats
  
  # 1. All habitats separated
  # levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")
  
  # 2. All urban habitats separated, all non-urban habitat together
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  # 3. All urban and non-urban habitats pooled together
  # levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")
  
  
  
  ## Insectivory ##
  #################
  
    # 8 morphology diversity within insectivorous (>=40% diet insects) 
  
  x1 <- subset(x,QE.insectiv>0)
  
  QE.insectiv = lme(QE.insectiv ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x1), method="ML")
  summary(QE.insectiv)
  testInteractions(QE.insectiv, pairwise="habitat.ordered")
  plot(QE.insectiv)
  qqnorm(QE.insectiv)
  QE.insectiv.I <- interactionMeans(QE.insectiv) # effect plots
  plot(QE.insectiv.I, errorbar="ci95")
  
  CR.insectiv = lme(CR.insectiv ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x1), method="ML")
  summary(CR.insectiv)  
  testInteractions(CR.insectiv, pairwise="habitat.ordered")
  plot(CR.insectiv)
  qqnorm(CR.insectiv)
  CR.insectiv.I <- interactionMeans(CR.insectiv) # effect plots
  plot(CR.insectiv.I, errorbar="ci95")
  
  meanD.insectiv = lme(insectiv.meanD ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x1), method="ML")
  summary(meanD.insectiv)  
  testInteractions(meanD.insectiv, pairwise="habitat.ordered")
  plot(meanD.insectiv)
  qqnorm(meanD.insectiv)
  meanD.insectiv.I <- interactionMeans(meanD.insectiv) # effect plots
  plot(meanD.insectiv.I, errorbar="ci95")
  
  Balance.insectiv = lme(insectiv.Balance ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x1), method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
  summary(Balance.insectiv)  
  testInteractions(Balance.insectiv, pairwise="habitat.ordered")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
  plot(Balance.insectiv)
  qqnorm(Balance.insectiv)
  Balance.insectiv.I <- interactionMeans(Balance.insectiv) # effect plots
  plot(Balance.insectiv.I, errorbar="ci95")
  
  
  
  a <- ggplot(QE.insectiv.I, aes(x= habitat.ordered, y=QE.insectiv.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.insectiv.I[,2]-QE.insectiv.I[,3], ymax=QE.insectiv.I[,2]+QE.insectiv.I[,3]), width=.2) +
    geom_point(data=QE.insectiv.I, mapping=aes(x=habitat.ordered, y=QE.insectiv.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "QE morphology within insectivorous", cex=16) +
    geom_text(aes(label= c("a","a","a","a","b")))
  
  b <- ggplot(CR.insectiv.I , aes(x= habitat.ordered, y=CR.insectiv.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.insectiv.I [,2]-CR.insectiv.I [,3], ymax=CR.insectiv.I [,2]+CR.insectiv.I [,3]), width=.2) +
    geom_point(data=CR.insectiv.I , mapping=aes(x=habitat.ordered, y=CR.insectiv.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, morphology within insectivorous ", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  c <- ggplot(meanD.insectiv.I, aes(x= habitat.ordered, y=meanD.insectiv.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.insectiv.I [,2]-meanD.insectiv.I [,3], ymax=meanD.insectiv.I [,2]+meanD.insectiv.I [,3]), width=.2) +
    geom_point(data=meanD.insectiv.I , mapping=aes(x=habitat.ordered, y=meanD.insectiv.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, morphology within insectivorous ", cex=16) +
    geom_text(aes(label= c("a","ab","ab","ab","b")))
  
  d <- ggplot(Balance.insectiv.I , aes(x= habitat.ordered, y=Balance.insectiv.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.insectiv.I [,2]-Balance.insectiv.I [,3], ymax=Balance.insectiv.I [,2]+Balance.insectiv.I [,3]), width=.2) +
    geom_point(data=Balance.insectiv.I , mapping=aes(x=habitat.ordered, y=Balance.insectiv.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, morphology within insectivorous ", cex=16) +
    geom_text(aes(label= c("ab","ab","ab","a","b")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_morphology_insectivory_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
 
  
  
  
  ## Granivory ##
  #################
  
  x2 <- subset(x,QE.seeds>0)
  
  # 8 morphology diversity within seedorous (>=40% diet insects) 
  
  QE.seeds = lme(QE.seeds ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x2), method="ML")
  summary(QE.seeds)
  testInteractions(QE.seeds, pairwise="habitat.ordered")
  plot(QE.seeds)
  qqnorm(QE.seeds)
  QE.seeds.I <- interactionMeans(QE.seeds) # effect plots
  plot(QE.seeds.I, errorbar="ci95")
  
  CR.seeds = lme(CR.seeds ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x2), method="ML")
  summary(CR.seeds)  
  testInteractions(CR.seeds, pairwise="habitat.ordered")
  plot(CR.seeds)
  qqnorm(CR.seeds)
  CR.seeds.I <- interactionMeans(CR.seeds) # effect plots
  plot(CR.seeds.I, errorbar="ci95")
  
  meanD.seeds = lme(seeds.meanD ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x2), method="ML")
  summary(meanD.seeds)  
  testInteractions(meanD.seeds, pairwise="habitat.ordered")
  plot(meanD.seeds)
  qqnorm(meanD.seeds)
  meanD.seeds.I <- interactionMeans(meanD.seeds) # effect plots
  plot(meanD.seeds.I, errorbar="ci95")
  
  Balance.seeds = lme(seeds.Balance ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x2), method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
  summary(Balance.seeds)  
  testInteractions(Balance.seeds, pairwise="habitat.ordered")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
  plot(Balance.seeds)
  qqnorm(Balance.seeds)
  Balance.seeds.I <- interactionMeans(Balance.seeds) # effect plots
  plot(Balance.seeds.I, errorbar="ci95")
  
  
  
  a <- ggplot(QE.seeds.I, aes(x= habitat.ordered, y=QE.seeds.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.seeds.I[,2]-QE.seeds.I[,3], ymax=QE.seeds.I[,2]+QE.seeds.I[,3]), width=.2) +
    geom_point(data=QE.seeds.I, mapping=aes(x=habitat.ordered, y=QE.seeds.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "QE morphology within granivorous", cex=16) +
    geom_text(aes(label= c("a","a","a","a","b")))
  
  b <- ggplot(CR.seeds.I , aes(x= habitat.ordered, y=CR.seeds.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.seeds.I [,2]-CR.seeds.I [,3], ymax=CR.seeds.I [,2]+CR.seeds.I [,3]), width=.2) +
    geom_point(data=CR.seeds.I , mapping=aes(x=habitat.ordered, y=CR.seeds.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, morphology within granivorous ", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  c <- ggplot(meanD.seeds.I, aes(x= habitat.ordered, y=meanD.seeds.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.seeds.I [,2]-meanD.seeds.I [,3], ymax=meanD.seeds.I [,2]+meanD.seeds.I [,3]), width=.2) +
    geom_point(data=meanD.seeds.I , mapping=aes(x=habitat.ordered, y=meanD.seeds.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, morphology within granivorous ", cex=16) +
    geom_text(aes(label= c("a","ab","ab","ab","b")))
  
  d <- ggplot(Balance.seeds.I , aes(x= habitat.ordered, y=Balance.seeds.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.seeds.I [,2]-Balance.seeds.I [,3], ymax=Balance.seeds.I [,2]+Balance.seeds.I [,3]), width=.2) +
    geom_point(data=Balance.seeds.I , mapping=aes(x=habitat.ordered, y=Balance.seeds.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, morphology within granivorous ", cex=16) +
    geom_text(aes(label= c("ab","ab","ab","a","b")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_morphology_seeds_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
  
  
  
  
  
  
  ## Frugivory ##
  #################
  
  x3 <- subset(x,QE.fruits>0)
  
  # 8 morphology diversity within seedorous (>=40% diet insects) 
  
  QE.fruits = lme(QE.fruits ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x2), method="ML")
  summary(QE.fruits)
  testInteractions(QE.fruits, pairwise="habitat.ordered")
  plot(QE.fruits)
  qqnorm(QE.fruits)
  QE.fruits.I <- interactionMeans(QE.fruits) # effect plots
  plot(QE.fruits.I, errorbar="ci95")
  
  CR.fruits = lme(CR.fruits ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x2), method="ML")
  summary(CR.fruits)  
  testInteractions(CR.fruits, pairwise="habitat.ordered")
  plot(CR.fruits)
  qqnorm(CR.fruits)
  CR.fruits.I <- interactionMeans(CR.fruits) # effect plots
  plot(CR.fruits.I, errorbar="ci95")
  
  meanD.fruits = lme(fruits.meanD ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x2), method="ML")
  summary(meanD.fruits)  
  testInteractions(meanD.fruits, pairwise="habitat.ordered")
  plot(meanD.fruits)
  qqnorm(meanD.fruits)
  meanD.fruits.I <- interactionMeans(meanD.fruits) # effect plots
  plot(meanD.fruits.I, errorbar="ci95")
  
  Balance.fruits = lme(fruits.Balance ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x2), method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
  summary(Balance.fruits)  
  testInteractions(Balance.fruits, pairwise="habitat.ordered")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
  plot(Balance.fruits)
  qqnorm(Balance.fruits)
  Balance.fruits.I <- interactionMeans(Balance.fruits) # effect plots
  plot(Balance.fruits.I, errorbar="ci95")
  
  
  
  a <- ggplot(QE.fruits.I, aes(x= habitat.ordered, y=QE.fruits.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.fruits.I[,2]-QE.fruits.I[,3], ymax=QE.fruits.I[,2]+QE.fruits.I[,3]), width=.2) +
    geom_point(data=QE.fruits.I, mapping=aes(x=habitat.ordered, y=QE.fruits.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "QE morphology within granivorous", cex=16) +
    geom_text(aes(label= c("a","a","a","a","b")))
  
  b <- ggplot(CR.fruits.I , aes(x= habitat.ordered, y=CR.fruits.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.fruits.I [,2]-CR.fruits.I [,3], ymax=CR.fruits.I [,2]+CR.fruits.I [,3]), width=.2) +
    geom_point(data=CR.fruits.I , mapping=aes(x=habitat.ordered, y=CR.fruits.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, morphology within granivorous ", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  c <- ggplot(meanD.fruits.I, aes(x= habitat.ordered, y=meanD.fruits.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.fruits.I [,2]-meanD.fruits.I [,3], ymax=meanD.fruits.I [,2]+meanD.fruits.I [,3]), width=.2) +
    geom_point(data=meanD.fruits.I , mapping=aes(x=habitat.ordered, y=meanD.fruits.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, morphology within granivorous ", cex=16) +
    geom_text(aes(label= c("a","ab","ab","ab","b")))
  
  d <- ggplot(Balance.fruits.I , aes(x= habitat.ordered, y=Balance.fruits.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.fruits.I [,2]-Balance.fruits.I [,3], ymax=Balance.fruits.I [,2]+Balance.fruits.I [,3]), width=.2) +
    geom_point(data=Balance.fruits.I , mapping=aes(x=habitat.ordered, y=Balance.fruits.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, morphology within granivorous ", cex=16) +
    geom_text(aes(label= c("ab","ab","ab","a","b")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_morphology_fruits_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
  
  
}




############# 9. Analyses for diet axes diversity, only natives ######
############################################################################

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Diet diversity metrics for communities natives.txt"))
  
  # We restrict the analyses to studies where there is information inside and outside the city
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  
  # Three ways to code habitats
  
  # 1. All habitats separated
  # levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")
  
  # 2. All urban habitats separated, all non-urban habitat together
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  # 3. All urban and non-urban habitats pooled together
  # levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")
  
  
  
  ## PCoA 1: Insectivory vs frugiv/granivourous ##
  ################################################
  
  {
  #  The first PCoA is negative for insectivorous and positive for graniv/frugiv
  
  QE.PCoA1 = lme(QE.PCoA1 ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")
  summary(QE.PCoA1)
  testInteractions(QE.PCoA1, pairwise="habitat.ordered")
  plot(QE.PCoA1)
  qqnorm(QE.PCoA1)
  QE.PCoA1.I <- interactionMeans(QE.PCoA1) # effect plots
  plot(QE.PCoA1.I, errorbar="ci95")
  
  CR.PCoA1 = lme(CR.PCoA1 ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")
  summary(CR.PCoA1)  
  testInteractions(CR.PCoA1, pairwise="habitat.ordered")
  plot(CR.PCoA1)
  qqnorm(CR.PCoA1)
  CR.PCoA1.I <- interactionMeans(CR.PCoA1) # effect plots
  plot(CR.PCoA1.I, errorbar="ci95")
  
  meanD.PCoA1 = lme(PCoA1.meanD ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")
  summary(meanD.PCoA1)  
  testInteractions(meanD.PCoA1, pairwise="habitat.ordered")
  plot(meanD.PCoA1)
  qqnorm(meanD.PCoA1)
  meanD.PCoA1.I <- interactionMeans(meanD.PCoA1) # effect plots
  plot(meanD.PCoA1.I, errorbar="ci95")
  
  Balance.PCoA1 = lme(PCoA1.Balance ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
  summary(Balance.PCoA1)  
  testInteractions(Balance.PCoA1, pairwise="habitat.ordered")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
  plot(Balance.PCoA1)
  qqnorm(Balance.PCoA1)
  Balance.PCoA1.I <- interactionMeans(Balance.PCoA1) # effect plots
  plot(Balance.PCoA1.I, errorbar="ci95")
  
  
  a <- ggplot(QE.PCoA1.I, aes(x= habitat.ordered, y=QE.PCoA1.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=QE.PCoA1.I[,2]-QE.PCoA1.I[,3], ymax=QE.PCoA1.I[,2]+QE.PCoA1.I[,3]), width=.2) +
    geom_point(data=QE.PCoA1.I, mapping=aes(x=habitat.ordered, y=QE.PCoA1.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "QE diet PCoA1", cex=16) +
    geom_text(aes(label= c("a","a","a","a","b")))
  
  b <- ggplot(CR.PCoA1.I , aes(x= habitat.ordered, y=CR.PCoA1.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.PCoA1.I [,2]-CR.PCoA1.I [,3], ymax=CR.PCoA1.I [,2]+CR.PCoA1.I [,3]), width=.2) +
    geom_point(data=CR.PCoA1.I , mapping=aes(x=habitat.ordered, y=CR.PCoA1.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Functional redundancies, diet PCoA1 ", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  c <- ggplot(meanD.PCoA1.I, aes(x= habitat.ordered, y=meanD.PCoA1.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.PCoA1.I [,2]-meanD.PCoA1.I [,3], ymax=meanD.PCoA1.I [,2]+meanD.PCoA1.I [,3]), width=.2) +
    geom_point(data=meanD.PCoA1.I , mapping=aes(x=habitat.ordered, y=meanD.PCoA1.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Average species dissimilarity, diet PCoA1", cex=16) +
    geom_text(aes(label= c("a","ab","ab","ab","b")))
  
  d <- ggplot(Balance.PCoA1.I , aes(x= habitat.ordered, y=Balance.PCoA1.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.PCoA1.I [,2]-Balance.PCoA1.I [,3], ymax=Balance.PCoA1.I [,2]+Balance.PCoA1.I [,3]), width=.2) +
    geom_point(data=Balance.PCoA1.I , mapping=aes(x=habitat.ordered, y=Balance.PCoA1.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance component, diet PCoA1", cex=16) +
    geom_text(aes(label= c("ab","ab","ab","a","b")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_morphology_PCoA1_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
  
  }
  
  

  ## PCoA 2: Insectivory vs granivourous ##
  #########################################
  
  {
    #  The first PCoA is negative for insectivorous and positive for graniv/frugiv
    
    QE.PCoA2 = lme(QE.PCoA2 ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")
    summary(QE.PCoA2)
    testInteractions(QE.PCoA2, pairwise="habitat.ordered")
    plot(QE.PCoA2)
    qqnorm(QE.PCoA2)
    QE.PCoA2.I <- interactionMeans(QE.PCoA2) # effect plots
    plot(QE.PCoA2.I, errorbar="ci95")
    
    CR.PCoA2 = lme(CR.PCoA2 ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")
    summary(CR.PCoA2)  
    testInteractions(CR.PCoA2, pairwise="habitat.ordered")
    plot(CR.PCoA2)
    qqnorm(CR.PCoA2)
    CR.PCoA2.I <- interactionMeans(CR.PCoA2) # effect plots
    plot(CR.PCoA2.I, errorbar="ci95")
    
    meanD.PCoA2 = lme(PCoA2.meanD ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")
    summary(meanD.PCoA2)  
    testInteractions(meanD.PCoA2, pairwise="habitat.ordered")
    plot(meanD.PCoA2)
    qqnorm(meanD.PCoA2)
    meanD.PCoA2.I <- interactionMeans(meanD.PCoA2) # effect plots
    plot(meanD.PCoA2.I, errorbar="ci95")
    
    Balance.PCoA2 = lme(PCoA2.Balance ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
    summary(Balance.PCoA2)  
    testInteractions(Balance.PCoA2, pairwise="habitat.ordered")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
    plot(Balance.PCoA2)
    qqnorm(Balance.PCoA2)
    Balance.PCoA2.I <- interactionMeans(Balance.PCoA2) # effect plots
    plot(Balance.PCoA2.I, errorbar="ci95")
    
    
    a <- ggplot(QE.PCoA2.I, aes(x= habitat.ordered, y=QE.PCoA2.I[,2])) + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
      geom_errorbar(aes(ymin=QE.PCoA2.I[,2]-QE.PCoA2.I[,3], ymax=QE.PCoA2.I[,2]+QE.PCoA2.I[,3]), width=.2) +
      geom_point(data=QE.PCoA2.I, mapping=aes(x=habitat.ordered, y=QE.PCoA2.I[,2]), size=8, shape=21, fill="white") +
      labs(x = "", y = "QE diet PCoA2", cex=16) +
      geom_text(aes(label= c("a","a","a","a","b")))
    
    b <- ggplot(CR.PCoA2.I , aes(x= habitat.ordered, y=CR.PCoA2.I [,2])) + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
      geom_errorbar(aes(ymin=CR.PCoA2.I [,2]-CR.PCoA2.I [,3], ymax=CR.PCoA2.I [,2]+CR.PCoA2.I [,3]), width=.2) +
      geom_point(data=CR.PCoA2.I , mapping=aes(x=habitat.ordered, y=CR.PCoA2.I [,2]), size=8, shape=21, fill="white") +
      labs(x = "", y = "Functional redundancies, diet PCoA2 ", cex=16) +
      geom_text(aes(label= c("a","a","a","a","a")))
    
    c <- ggplot(meanD.PCoA2.I, aes(x= habitat.ordered, y=meanD.PCoA2.I [,2])) + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
      geom_errorbar(aes(ymin=meanD.PCoA2.I [,2]-meanD.PCoA2.I [,3], ymax=meanD.PCoA2.I [,2]+meanD.PCoA2.I [,3]), width=.2) +
      geom_point(data=meanD.PCoA2.I , mapping=aes(x=habitat.ordered, y=meanD.PCoA2.I [,2]), size=8, shape=21, fill="white") +
      labs(x = "", y = "Average species dissimilarity, diet PCoA2", cex=16) +
      geom_text(aes(label= c("a","ab","ab","ab","b")))
    
    d <- ggplot(Balance.PCoA2.I , aes(x= habitat.ordered, y=Balance.PCoA2.I [,2])) + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
      geom_errorbar(aes(ymin=Balance.PCoA2.I [,2]-Balance.PCoA2.I [,3], ymax=Balance.PCoA2.I [,2]+Balance.PCoA2.I [,3]), width=.2) +
      geom_point(data=Balance.PCoA2.I , mapping=aes(x=habitat.ordered, y=Balance.PCoA2.I [,2]), size=8, shape=21, fill="white") +
      labs(x = "", y = "Balance component, diet PCoA2", cex=16) +
      geom_text(aes(label= c("ab","ab","ab","a","b")))
    
    tiff(paste0(GoogleFigs,"/plot_FD_morphology_PCoA2_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
    ggplot2.multiplot(a,b,c,d)
    dev.off()
    
  }
  

  
  ## PCoA 3: granivourous ##
  ##########################
  
  {
    #  The first PCoA is negative for insectivorous and positive for graniv/frugiv
    
    QE.PCoA3 = lme(QE.PCoA3 ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")
    summary(QE.PCoA3)
    testInteractions(QE.PCoA3, pairwise="habitat.ordered")
    plot(QE.PCoA3)
    qqnorm(QE.PCoA3)
    QE.PCoA3.I <- interactionMeans(QE.PCoA3) # effect plots
    plot(QE.PCoA3.I, errorbar="ci95")
    
    CR.PCoA3 = lme(CR.PCoA3 ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")
    summary(CR.PCoA3)  
    testInteractions(CR.PCoA3, pairwise="habitat.ordered")
    plot(CR.PCoA3)
    qqnorm(CR.PCoA3)
    CR.PCoA3.I <- interactionMeans(CR.PCoA3) # effect plots
    plot(CR.PCoA3.I, errorbar="ci95")
    
    meanD.PCoA3 = lme(PCoA3.meanD ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")
    summary(meanD.PCoA3)  
    testInteractions(meanD.PCoA3, pairwise="habitat.ordered")
    plot(meanD.PCoA3)
    qqnorm(meanD.PCoA3)
    meanD.PCoA3.I <- interactionMeans(meanD.PCoA3) # effect plots
    plot(meanD.PCoA3.I, errorbar="ci95")
    
    Balance.PCoA3 = lme(PCoA3.Balance ~ habitat.ordered, random = ~ 1|country/city, data = na.omit(x), method="ML")  # if the most abundant species tended to be dissimilar = positive Balance component
    summary(Balance.PCoA3)  
    testInteractions(Balance.PCoA3, pairwise="habitat.ordered")                                            # a lack of association between abundance and dissimilarity = Balance component close to zero
    plot(Balance.PCoA3)
    qqnorm(Balance.PCoA3)
    Balance.PCoA3.I <- interactionMeans(Balance.PCoA3) # effect plots
    plot(Balance.PCoA3.I, errorbar="ci95")
    
    
    a <- ggplot(QE.PCoA3.I, aes(x= habitat.ordered, y=QE.PCoA3.I[,2])) + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
      geom_errorbar(aes(ymin=QE.PCoA3.I[,2]-QE.PCoA3.I[,3], ymax=QE.PCoA3.I[,2]+QE.PCoA3.I[,3]), width=.2) +
      geom_point(data=QE.PCoA3.I, mapping=aes(x=habitat.ordered, y=QE.PCoA3.I[,2]), size=8, shape=21, fill="white") +
      labs(x = "", y = "QE diet PCoA3", cex=16) +
      geom_text(aes(label= c("a","a","a","a","b")))
    
    b <- ggplot(CR.PCoA3.I , aes(x= habitat.ordered, y=CR.PCoA3.I [,2])) + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
      geom_errorbar(aes(ymin=CR.PCoA3.I [,2]-CR.PCoA3.I [,3], ymax=CR.PCoA3.I [,2]+CR.PCoA3.I [,3]), width=.2) +
      geom_point(data=CR.PCoA3.I , mapping=aes(x=habitat.ordered, y=CR.PCoA3.I [,2]), size=8, shape=21, fill="white") +
      labs(x = "", y = "Functional redundancies, diet PCoA3 ", cex=16) +
      geom_text(aes(label= c("a","a","a","a","a")))
    
    c <- ggplot(meanD.PCoA3.I, aes(x= habitat.ordered, y=meanD.PCoA3.I [,2])) + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
      geom_errorbar(aes(ymin=meanD.PCoA3.I [,2]-meanD.PCoA3.I [,3], ymax=meanD.PCoA3.I [,2]+meanD.PCoA3.I [,3]), width=.2) +
      geom_point(data=meanD.PCoA3.I , mapping=aes(x=habitat.ordered, y=meanD.PCoA3.I [,2]), size=8, shape=21, fill="white") +
      labs(x = "", y = "Average species dissimilarity, diet PCoA3", cex=16) +
      geom_text(aes(label= c("a","ab","ab","ab","b")))
    
    d <- ggplot(Balance.PCoA3.I , aes(x= habitat.ordered, y=Balance.PCoA3.I [,2])) + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
      geom_errorbar(aes(ymin=Balance.PCoA3.I [,2]-Balance.PCoA3.I [,3], ymax=Balance.PCoA3.I [,2]+Balance.PCoA3.I [,3]), width=.2) +
      geom_point(data=Balance.PCoA3.I , mapping=aes(x=habitat.ordered, y=Balance.PCoA3.I [,2]), size=8, shape=21, fill="white") +
      labs(x = "", y = "Balance component, diet PCoA3", cex=16) +
      geom_text(aes(label= c("ab","ab","ab","a","b")))
    
    tiff(paste0(GoogleFigs,"/plot_FD_morphology_PCoA3_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
    ggplot2.multiplot(a,b,c,d)
    dev.off()
    
  }
  
  
  
  
}
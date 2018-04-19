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


### ANALYSES: ###

## 1. Analyses for morphological diversity
## 2. Analyses for morphological diversity, only natives




########### 1. Analyses for morphological diversity ##############
##################################################################

{## Import biodiversity metrics for communities

<<<<<<< HEAD
# x <- read.table("/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis/Morphological diversity metrics for communities.txt") # metrics estimated including all species
# x <- read.table("/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis/Morphological diversity metrics for communities natives.txt")  # metrics estimated excluding exotics
=======
#x <- read.table("/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis/Morphological diversity metrics for communities.txt") # metrics estimated including all species
x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities.txt"))
# x <- read.table("Morphological diversity metrics for communities natives.txt")  # metrics estimated excluding exotics
>>>>>>> 1483dfcd4dc55a894b3bfccc615633c05db0fe09


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
  geom_text(aes(label= c("a","a","ab","b","c")))

c <- ggplot(AEveness.I, aes(x= habitat.ordered, y=AEveness.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=AEveness.I[,2]-AEveness.I[,3], ymax=AEveness.I[,2]+AEveness.I[,3]), width=.2) +
  geom_point(data=AEveness.I, mapping=aes(x=habitat.ordered, y=AEveness.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Abundance evenness index", cex=16) +
  geom_text(aes(label= c("a","a","a","b","b")))

tiff(paste0(GoogleFigs,"/plot_taxonomic_diversity.tiff"), width = 9, height = 8, units = 'in', res = 200)
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
  labs(x = "", y = "QE all traits", cex=16) +
  geom_text(aes(label= c("a","b","ab","b","b")))

b <- ggplot(CR.all.morph.I, aes(x= habitat.ordered, y=CR.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.all.morph.I[,2]-CR.all.morph.I[,3], ymax=CR.all.morph.I[,2]+CR.all.morph.I[,3]), width=.2) +
  geom_point(data=CR.all.morph.I, mapping=aes(x=habitat.ordered, y=CR.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CR all traits", cex=16) +
  geom_text(aes(label= c("ab","c","abc","ac","b")))

c <- ggplot(meanD.all.morph.I, aes(x= habitat.ordered, y=meanD.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.all.morph.I[,2]-meanD.all.morph.I[,3], ymax=meanD.all.morph.I[,2]+meanD.all.morph.I[,3]), width=.2) +
  geom_point(data=meanD.all.morph.I, mapping=aes(x=habitat.ordered, y=meanD.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "meanD all traits", cex=16) +
  geom_text(aes(label= c("a","b","ab","ab","a")))
  
d <- ggplot(Balance.all.morph.I, aes(x= habitat.ordered, y=Balance.all.morph.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.all.morph.I[,2]-Balance.all.morph.I[,3], ymax=Balance.all.morph.I[,2]+Balance.all.morph.I[,3]), width=.2) +
  geom_point(data=Balance.all.morph.I, mapping=aes(x=habitat.ordered, y=Balance.all.morph.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance all traits", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

tiff(paste0(GoogleFigs,"plot_FD_all_traits.tiff"), width = 11, height = 8, units = 'in', res = 200)
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
  labs(x = "", y = "QE PCA3", cex=16) +
  geom_text(aes(label= c("a","b","ab","b","b")))

b <- ggplot(CR.PCA3.I, aes(x= habitat.ordered, y=CR.PCA3.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.PCA3.I[,2]-CR.PCA3.I[,3], ymax=CR.PCA3.I[,2]+CR.PCA3.I[,3]), width=.2) +
  geom_point(data=CR.PCA3.I, mapping=aes(x=habitat.ordered, y=CR.PCA3.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CR PCA3", cex=16) +
  geom_text(aes(label= c("a","b","ab","ab","a")))

c <- ggplot(meanD.PCA3.I, aes(x= habitat.ordered, y=meanD.PCA3.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.PCA3.I[,2]-meanD.PCA3.I[,3], ymax=meanD.PCA3.I[,2]+meanD.PCA3.I[,3]), width=.2) +
  geom_point(data=meanD.PCA3.I, mapping=aes(x=habitat.ordered, y=meanD.PCA3.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "meanD PCA3", cex=16) +
  geom_text(aes(label= c("a","b","ab","ab","a")))

d <- ggplot(Balance.PCA3.I, aes(x= habitat.ordered, y=Balance.PCA3.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.PCA3.I[,2]-Balance.PCA3.I[,3], ymax=Balance.PCA3.I[,2]+Balance.PCA3.I[,3]), width=.2) +
  geom_point(data=Balance.PCA3.I, mapping=aes(x=habitat.ordered, y=Balance.PCA3.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance PCA3", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

#tiff(paste0(GoogleFigs,"/plot_FD_PCA3.tiff"), width = 11, height = 8, units = 'in', res = 200)
tiff(paste0(GoogleFigs,"plot_FD_PCA3.tiff"), width = 11, height = 8, units = 'in', res = 200)
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
  labs(x = "", y = "QE beak", cex=16) +
  geom_text(aes(label= c("a","a","ab","b","b")))

b <- ggplot(CR.beak.I, aes(x= habitat.ordered, y=CR.beak.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.beak.I[,2]-CR.beak.I[,3], ymax=CR.beak.I[,2]+CR.beak.I[,3]), width=.2) +
  geom_point(data=CR.beak.I, mapping=aes(x=habitat.ordered, y=CR.beak.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CR beak", cex=16) +
  geom_text(aes(label= c("a","a","a","b","b")))

c <- ggplot(meanD.beak.I, aes(x= habitat.ordered, y=meanD.beak.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.beak.I[,2]-meanD.beak.I[,3], ymax=meanD.beak.I[,2]+meanD.beak.I[,3]), width=.2) +
  geom_point(data=meanD.beak.I, mapping=aes(x=habitat.ordered, y=meanD.beak.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "meanD beak", cex=16) +
  geom_text(aes(label= c("a","a","ab","b","c")))

d <- ggplot(Balance.beak.I, aes(x= habitat.ordered, y=Balance.beak.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.beak.I[,2]-Balance.beak.I[,3], ymax=Balance.beak.I[,2]+Balance.beak.I[,3]), width=.2) +
  geom_point(data=Balance.beak.I, mapping=aes(x=habitat.ordered, y=Balance.beak.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance beak", cex=16) +
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
  labs(x = "", y = "QE locomotory system", cex=16) +
  geom_text(aes(label= c("a","ab","ab","ab","b")))

b <- ggplot(CR.locom.I , aes(x= habitat.ordered, y=CR.locom.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.locom.I [,2]-CR.locom.I [,3], ymax=CR.locom.I [,2]+CR.locom.I [,3]), width=.2) +
  geom_point(data=CR.locom.I , mapping=aes(x=habitat.ordered, y=CR.locom.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CR locomotory system", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

c <- ggplot(meanD.locom.I , aes(x= habitat.ordered, y=meanD.locom.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.locom.I [,2]-meanD.locom.I [,3], ymax=meanD.locom.I [,2]+meanD.locom.I [,3]), width=.2) +
  geom_point(data=meanD.locom.I , mapping=aes(x=habitat.ordered, y=meanD.locom.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "meanD locomotory system", cex=16) +
  geom_text(aes(label= c("ab","ab","a","b","b")))

d <- ggplot(Balance.locom.I , aes(x= habitat.ordered, y=Balance.locom.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.locom.I [,2]-Balance.locom.I [,3], ymax=Balance.locom.I [,2]+Balance.locom.I [,3]), width=.2) +
  geom_point(data=Balance.locom.I , mapping=aes(x=habitat.ordered, y=Balance.locom.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance locomotory system", cex=16) +
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
  labs(x = "", y = "QE body size", cex=16) +
  geom_text(aes(label= c("a","b","a","b","b")))

b <- ggplot(CR.size.I , aes(x= habitat.ordered, y=CR.size.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.size.I [,2]-CR.size.I [,3], ymax=CR.size.I [,2]+CR.size.I [,3]), width=.2) +
  geom_point(data=CR.size.I , mapping=aes(x=habitat.ordered, y=CR.size.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CR body size", cex=16) +
  geom_text(aes(label= c("a","b","ab","b","a")))

c <- ggplot(meanD.size.I , aes(x= habitat.ordered, y=meanD.size.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.size.I [,2]-meanD.size.I [,3], ymax=meanD.size.I [,2]+meanD.size.I [,3]), width=.2) +
  geom_point(data=meanD.size.I , mapping=aes(x=habitat.ordered, y=meanD.size.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "meanD body size", cex=16) +
  geom_text(aes(label= c("a","b","abc","abc","ac")))

d <- ggplot(Balance.size.I , aes(x= habitat.ordered, y=Balance.size.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.size.I [,2]-Balance.size.I [,3], ymax=Balance.size.I [,2]+Balance.size.I [,3]), width=.2) +
  geom_point(data=Balance.size.I , mapping=aes(x=habitat.ordered, y=Balance.size.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance body size", cex=16) +
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
  labs(x = "", y = "QE wing hand", cex=16) +
  geom_text(aes(label= c("a","ab","ab","ab","b")))

b <- ggplot(CR.winghand.I , aes(x= habitat.ordered, y=CR.winghand.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.winghand.I [,2]-CR.winghand.I [,3], ymax=CR.winghand.I [,2]+CR.winghand.I [,3]), width=.2) +
  geom_point(data=CR.winghand.I , mapping=aes(x=habitat.ordered, y=CR.winghand.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CR wing hand ", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

c <- ggplot(meanD.winghand.I , aes(x= habitat.ordered, y=meanD.winghand.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.winghand.I [,2]-meanD.winghand.I [,3], ymax=meanD.winghand.I [,2]+meanD.winghand.I [,3]), width=.2) +
  geom_point(data=meanD.winghand.I , mapping=aes(x=habitat.ordered, y=meanD.winghand.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "meanD wing hand ", cex=16) +
  geom_text(aes(label= c("a","a","a","a","a")))

d <- ggplot(Balance.winghand.I , aes(x= habitat.ordered, y=Balance.winghand.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.winghand.I [,2]-Balance.winghand.I [,3], ymax=Balance.winghand.I [,2]+Balance.winghand.I [,3]), width=.2) +
  geom_point(data=Balance.winghand.I , mapping=aes(x=habitat.ordered, y=Balance.winghand.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance wing hand ", cex=16) +
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
  labs(x = "", y = "QE Erickson’s phylogeny", cex=16) +
  geom_text(aes(label= c("ac","b","ac","abc","bc")))

b <- ggplot(CR.phyE.I , aes(x= habitat.ordered, y=CR.phyE.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.phyE.I [,2]-CR.phyE.I [,3], ymax=CR.phyE.I [,2]+CR.phyE.I [,3]), width=.2) +
  geom_point(data=CR.phyE.I , mapping=aes(x=habitat.ordered, y=CR.phyE.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CR Erickson’s phylogeny ", cex=16) +
  geom_text(aes(label= c("ab","a","bc","bc","c")))

c <- ggplot(meanD.phyE.I , aes(x= habitat.ordered, y=meanD.phyE.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.phyE.I [,2]-meanD.phyE.I [,3], ymax=meanD.phyE.I [,2]+meanD.phyE.I [,3]), width=.2) +
  geom_point(data=meanD.phyE.I , mapping=aes(x=habitat.ordered, y=meanD.phyE.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "meanD Erickson’s phylogeny ", cex=16) +
  geom_text(aes(label= c("ab","a","b","ab","ab")))

d <- ggplot(Balance.phyE.I , aes(x= habitat.ordered, y=Balance.phyE.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.phyE.I [,2]-Balance.phyE.I [,3], ymax=Balance.phyE.I [,2]+Balance.phyE.I [,3]), width=.2) +
  geom_point(data=Balance.phyE.I , mapping=aes(x=habitat.ordered, y=Balance.phyE.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance Erickson’s phylogeny ", cex=16) +
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
  labs(x = "", y = "QE Hackett’s phylogeny", cex=16) +
  geom_text(aes(label= c("a","b","ac","abc","b")))

b <- ggplot(CR.phyH.I , aes(x= habitat.ordered, y=CR.phyH.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=CR.phyH.I [,2]-CR.phyH.I [,3], ymax=CR.phyH.I [,2]+CR.phyH.I [,3]), width=.2) +
  geom_point(data=CR.phyH.I , mapping=aes(x=habitat.ordered, y=CR.phyH.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "CR Hackett’s phylogeny ", cex=16) +
  geom_text(aes(label= c("a","b","ac","abc","c")))

c <- ggplot(meanD.phyH.I, aes(x= habitat.ordered, y=meanD.phyH.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=meanD.phyH.I [,2]-meanD.phyH.I [,3], ymax=meanD.phyH.I [,2]+meanD.phyH.I [,3]), width=.2) +
  geom_point(data=meanD.phyH.I , mapping=aes(x=habitat.ordered, y=meanD.phyH.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "meanD Hackett’s phylogeny ", cex=16) +
  geom_text(aes(label= c("abc","b","c","abc","abc")))

d <- ggplot(Balance.phyH.I , aes(x= habitat.ordered, y=Balance.phyH.I [,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Balance.phyH.I [,2]-Balance.phyH.I [,3], ymax=Balance.phyH.I [,2]+Balance.phyH.I [,3]), width=.2) +
  geom_point(data=Balance.phyH.I , mapping=aes(x=habitat.ordered, y=Balance.phyH.I [,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Balance Hackett’s phylogeny ", cex=16) +
  geom_text(aes(label= c("a","a","ab","bc","c")))

tiff(paste0(GoogleFigs,"/plot_FD_Hackett_phylogeny.tiff"), width = 11, height = 8, units = 'in', res = 200)
ggplot2.multiplot(a,b,c,d)
dev.off()
}

########### 2. Analyses for morphological diversity, omly natives ##############
################################################################################

{## Import biodiversity metrics for communities
  
  # x <- read.table("/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis/Morphological diversity metrics for communities natives.txt")  # metrics estimated excluding exotics
  x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities natives.txt"))
  
  
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
    labs(x = "", y = "QE all traits", cex=16) +
    geom_text(aes(label= c("a","ab","ab","b","b")))
  
  b <- ggplot(CR.all.morph.I, aes(x= habitat.ordered, y=CR.all.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.all.morph.I[,2]-CR.all.morph.I[,3], ymax=CR.all.morph.I[,2]+CR.all.morph.I[,3]), width=.2) +
    geom_point(data=CR.all.morph.I, mapping=aes(x=habitat.ordered, y=CR.all.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "CR all traits", cex=16) +
    geom_text(aes(label= c("ab","a","ab","ab","b")))
  
  c <- ggplot(meanD.all.morph.I, aes(x= habitat.ordered, y=meanD.all.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.all.morph.I[,2]-meanD.all.morph.I[,3], ymax=meanD.all.morph.I[,2]+meanD.all.morph.I[,3]), width=.2) +
    geom_point(data=meanD.all.morph.I, mapping=aes(x=habitat.ordered, y=meanD.all.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "meanD all traits", cex=16) +
    geom_text(aes(label= c("a","b","ab","ab","a")))
  
  d <- ggplot(Balance.all.morph.I, aes(x= habitat.ordered, y=Balance.all.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.all.morph.I[,2]-Balance.all.morph.I[,3], ymax=Balance.all.morph.I[,2]+Balance.all.morph.I[,3]), width=.2) +
    geom_point(data=Balance.all.morph.I, mapping=aes(x=habitat.ordered, y=Balance.all.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance all traits", cex=16) +
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
    labs(x = "", y = "QE PCA3", cex=16) +
    geom_text(aes(label= c("a","b","ab","b","b")))
  
  b <- ggplot(CR.PCA3.I, aes(x= habitat.ordered, y=CR.PCA3.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.PCA3.I[,2]-CR.PCA3.I[,3], ymax=CR.PCA3.I[,2]+CR.PCA3.I[,3]), width=.2) +
    geom_point(data=CR.PCA3.I, mapping=aes(x=habitat.ordered, y=CR.PCA3.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "CR PCA3", cex=16) +
    geom_text(aes(label= c("ab","a","ab","ab","b")))
  
  c <- ggplot(meanD.PCA3.I, aes(x= habitat.ordered, y=meanD.PCA3.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.PCA3.I[,2]-meanD.PCA3.I[,3], ymax=meanD.PCA3.I[,2]+meanD.PCA3.I[,3]), width=.2) +
    geom_point(data=meanD.PCA3.I, mapping=aes(x=habitat.ordered, y=meanD.PCA3.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "meanD PCA3", cex=16) +
    geom_text(aes(label= c("a","b","ab","ab","a")))
  
  d <- ggplot(Balance.PCA3.I, aes(x= habitat.ordered, y=Balance.PCA3.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.PCA3.I[,2]-Balance.PCA3.I[,3], ymax=Balance.PCA3.I[,2]+Balance.PCA3.I[,3]), width=.2) +
    geom_point(data=Balance.PCA3.I, mapping=aes(x=habitat.ordered, y=Balance.PCA3.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance PCA3", cex=16) +
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
    labs(x = "", y = "QE beak", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  b <- ggplot(CR.beak.I, aes(x= habitat.ordered, y=CR.beak.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.beak.I[,2]-CR.beak.I[,3], ymax=CR.beak.I[,2]+CR.beak.I[,3]), width=.2) +
    geom_point(data=CR.beak.I, mapping=aes(x=habitat.ordered, y=CR.beak.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "CR beak", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b")))
  
  c <- ggplot(meanD.beak.I, aes(x= habitat.ordered, y=meanD.beak.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.beak.I[,2]-meanD.beak.I[,3], ymax=meanD.beak.I[,2]+meanD.beak.I[,3]), width=.2) +
    geom_point(data=meanD.beak.I, mapping=aes(x=habitat.ordered, y=meanD.beak.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "meanD beak", cex=16) +
    geom_text(aes(label= c("a","ab","ab","ab","bc")))
  
  d <- ggplot(Balance.beak.I, aes(x= habitat.ordered, y=Balance.beak.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.beak.I[,2]-Balance.beak.I[,3], ymax=Balance.beak.I[,2]+Balance.beak.I[,3]), width=.2) +
    geom_point(data=Balance.beak.I, mapping=aes(x=habitat.ordered, y=Balance.beak.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance beak", cex=16) +
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
    labs(x = "", y = "QE locomotory system", cex=16) +
    geom_text(aes(label= c("a","ab","ab","ab","b")))
  
  b <- ggplot(CR.locom.I , aes(x= habitat.ordered, y=CR.locom.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.locom.I [,2]-CR.locom.I [,3], ymax=CR.locom.I [,2]+CR.locom.I [,3]), width=.2) +
    geom_point(data=CR.locom.I , mapping=aes(x=habitat.ordered, y=CR.locom.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "CR locomotory system", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  c <- ggplot(meanD.locom.I , aes(x= habitat.ordered, y=meanD.locom.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.locom.I [,2]-meanD.locom.I [,3], ymax=meanD.locom.I [,2]+meanD.locom.I [,3]), width=.2) +
    geom_point(data=meanD.locom.I , mapping=aes(x=habitat.ordered, y=meanD.locom.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "meanD locomotory system", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  d <- ggplot(Balance.locom.I , aes(x= habitat.ordered, y=Balance.locom.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.locom.I [,2]-Balance.locom.I [,3], ymax=Balance.locom.I [,2]+Balance.locom.I [,3]), width=.2) +
    geom_point(data=Balance.locom.I , mapping=aes(x=habitat.ordered, y=Balance.locom.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance locomotory system", cex=16) +
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
    labs(x = "", y = "QE body size", cex=16) +
    geom_text(aes(label= c("a","ab","ab","b","b")))
  
  b <- ggplot(CR.size.I , aes(x= habitat.ordered, y=CR.size.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.size.I [,2]-CR.size.I [,3], ymax=CR.size.I [,2]+CR.size.I [,3]), width=.2) +
    geom_point(data=CR.size.I , mapping=aes(x=habitat.ordered, y=CR.size.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "CR body size", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  c <- ggplot(meanD.size.I , aes(x= habitat.ordered, y=meanD.size.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.size.I [,2]-meanD.size.I [,3], ymax=meanD.size.I [,2]+meanD.size.I [,3]), width=.2) +
    geom_point(data=meanD.size.I , mapping=aes(x=habitat.ordered, y=meanD.size.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "meanD body size", cex=16) +
    geom_text(aes(label= c("a","b","ab","ab","a")))
  
  d <- ggplot(Balance.size.I , aes(x= habitat.ordered, y=Balance.size.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.size.I [,2]-Balance.size.I [,3], ymax=Balance.size.I [,2]+Balance.size.I [,3]), width=.2) +
    geom_point(data=Balance.size.I , mapping=aes(x=habitat.ordered, y=Balance.size.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance body size", cex=16) +
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
    labs(x = "", y = "QE wing hand", cex=16) +
    geom_text(aes(label= c("a","ab","ab","ab","b")))
  
  b <- ggplot(CR.winghand.I , aes(x= habitat.ordered, y=CR.winghand.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.winghand.I [,2]-CR.winghand.I [,3], ymax=CR.winghand.I [,2]+CR.winghand.I [,3]), width=.2) +
    geom_point(data=CR.winghand.I , mapping=aes(x=habitat.ordered, y=CR.winghand.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "CR wing hand ", cex=16) +
    geom_text(aes(label= c("ab","ab","a","ab","b")))
  
  c <- ggplot(meanD.winghand.I , aes(x= habitat.ordered, y=meanD.winghand.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.winghand.I [,2]-meanD.winghand.I [,3], ymax=meanD.winghand.I [,2]+meanD.winghand.I [,3]), width=.2) +
    geom_point(data=meanD.winghand.I , mapping=aes(x=habitat.ordered, y=meanD.winghand.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "meanD wing hand ", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  d <- ggplot(Balance.winghand.I , aes(x= habitat.ordered, y=Balance.winghand.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.winghand.I [,2]-Balance.winghand.I [,3], ymax=Balance.winghand.I [,2]+Balance.winghand.I [,3]), width=.2) +
    geom_point(data=Balance.winghand.I , mapping=aes(x=habitat.ordered, y=Balance.winghand.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance wing hand ", cex=16) +
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
    labs(x = "", y = "QE Erickson’s phylogeny", cex=16) +
    geom_text(aes(label= c("a","b","a","ab","b")))
  
  b <- ggplot(CR.phyE.I , aes(x= habitat.ordered, y=CR.phyE.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.phyE.I [,2]-CR.phyE.I [,3], ymax=CR.phyE.I [,2]+CR.phyE.I [,3]), width=.2) +
    geom_point(data=CR.phyE.I , mapping=aes(x=habitat.ordered, y=CR.phyE.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "CR Erickson’s phylogeny ", cex=16) +
    geom_text(aes(label= c("ab","b","a","a","a")))
  
  c <- ggplot(meanD.phyE.I , aes(x= habitat.ordered, y=meanD.phyE.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.phyE.I [,2]-meanD.phyE.I [,3], ymax=meanD.phyE.I [,2]+meanD.phyE.I [,3]), width=.2) +
    geom_point(data=meanD.phyE.I , mapping=aes(x=habitat.ordered, y=meanD.phyE.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "meanD Erickson’s phylogeny ", cex=16) +
    geom_text(aes(label= c("ab","a","b","ab","ab")))
  
  d <- ggplot(Balance.phyE.I , aes(x= habitat.ordered, y=Balance.phyE.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.phyE.I [,2]-Balance.phyE.I [,3], ymax=Balance.phyE.I [,2]+Balance.phyE.I [,3]), width=.2) +
    geom_point(data=Balance.phyE.I , mapping=aes(x=habitat.ordered, y=Balance.phyE.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance Erickson’s phylogeny ", cex=16) +
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
    labs(x = "", y = "QE Hackett’s phylogeny", cex=16) +
    geom_text(aes(label= c("a","b","ac","abc","b")))
  
  b <- ggplot(CR.phyH.I , aes(x= habitat.ordered, y=CR.phyH.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=CR.phyH.I [,2]-CR.phyH.I [,3], ymax=CR.phyH.I [,2]+CR.phyH.I [,3]), width=.2) +
    geom_point(data=CR.phyH.I , mapping=aes(x=habitat.ordered, y=CR.phyH.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "CR Hackett’s phylogeny ", cex=16) +
    geom_text(aes(label= c("a","a","b","b","b")))  # check significances
  
  c <- ggplot(meanD.phyH.I, aes(x= habitat.ordered, y=meanD.phyH.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=meanD.phyH.I [,2]-meanD.phyH.I [,3], ymax=meanD.phyH.I [,2]+meanD.phyH.I [,3]), width=.2) +
    geom_point(data=meanD.phyH.I , mapping=aes(x=habitat.ordered, y=meanD.phyH.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "meanD Hackett’s phylogeny ", cex=16) +
    geom_text(aes(label= c("a","a","a","a","a")))
  
  d <- ggplot(Balance.phyH.I , aes(x= habitat.ordered, y=Balance.phyH.I [,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Balance.phyH.I [,2]-Balance.phyH.I [,3], ymax=Balance.phyH.I [,2]+Balance.phyH.I [,3]), width=.2) +
    geom_point(data=Balance.phyH.I , mapping=aes(x=habitat.ordered, y=Balance.phyH.I [,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Balance Hackett’s phylogeny ", cex=16) +
    geom_text(aes(label= c("a","a","ab","b","b")))
  
  tiff(paste0(GoogleFigs,"/plot_FD_Hackett_phylogeny_natives.tiff"), width = 11, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,b,c,d)
  dev.off()
}

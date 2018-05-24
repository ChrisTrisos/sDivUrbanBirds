###########################################################################################################
################ 17. Analyses of species-level functional uniqueness and vulnerability ####################
###########################################################################################################

dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)
dat$community <- factor(dat$community)
diet.Kbar <- read.table(paste0(workingData,"/vulnerability_species_diet_natives.txt"))
morphology.Kbar <- read.table(paste0(workingData,"/vulnerability_species_morphology_natives.txt"))

# Amravati

C1 <- subset(morphology.Kbar, community=="1")
hist(C1$Kbar.morphology)
C4 <- subset(morphology.Kbar, community=="4")
hist(C4$Kbar.morphology)


# Barcelona

C5 <- subset(morphology.Kbar, community=="5")
hist(C5$Kbar.morphology)
C6 <- subset(morphology.Kbar, community=="6")
hist(C6$Kbar.morphology)


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



## Tests of the effect of urbanization
################################################################

library(nlme)

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


a <- ggplot(Kbar.diet.skew.I, aes(x= habitat.ordered, y=Kbar.diet.skew.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Kbar.diet.skew.I[,2]-Kbar.diet.skew.I[,3], ymax=Kbar.diet.skew.I[,2]+Kbar.diet.skew.I[,3]), width=.2) +
  geom_point(data=Kbar.diet.skew.I, mapping=aes(x=habitat.ordered, y=Kbar.diet.skew.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Kbar skeweness, diet", cex=16) +
  geom_text(aes(label= c("a","a","a","b","b")))

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


b <- ggplot(Kbar.morph.skew.I, aes(x= habitat.ordered, y=Kbar.morph.skew.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Kbar.morph.skew.I[,2]-Kbar.morph.skew.I[,3], ymax=Kbar.morph.skew.I[,2]+Kbar.morph.skew.I[,3]), width=.2) +
  geom_point(data=Kbar.morph.skew.I, mapping=aes(x=habitat.ordered, y=Kbar.morph.skew.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Kbar skeweness, morphology", cex=16) +
  geom_text(aes(label= c("a","a","ab","a","b")))


tiff(paste0(GoogleFigs,"/plot_skewness.natives.tiff"), width = 9, height = 11, units = 'in', res = 200)
ggplot2.multiplot(a,b, cols=1)
dev.off()
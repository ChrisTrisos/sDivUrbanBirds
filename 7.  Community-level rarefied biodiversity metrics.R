
######################################################################################
### 7.	Community level rarefied estimations of taxonomic and phylogenetic metrics ###
######################################################################################


### We use rarefications to:

  # 1) asssess whether species richness is well estimated by the sampling effort of the studies (package iNEXT)
  # 2a) estimate rarefied measures of species richness, Simpson's index and phylogenetic diversity
  # 2b) and test for differences across habitats (library BAT)


### Inputs:
  # A list of abundance (number of individuals) vectors per community (iNEXT)
  # A list of species*community for cities, with abundances (number of individuals)




#### Data preparation ###
#########################

dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)

# dat <- subset(dat,status=="native")  # if we want to exclude exotics
# dat[] <- lapply(dat0, function(x) if(is.factor(x)) factor(x) else x)

dat01 <- subset(dat, used.urban.nonurban=="yes") # if we focus on communities along urbanisation gradients
dat01[] <- lapply(dat01, function(x) if(is.factor(x)) factor(x) else x)
dat01 <- na.omit(dat01[,-c(14,19)])

dat02 <- dat01[dat01$Individuals>0,] # we exclude species absent
dat02[] <- lapply(dat02, function(x) if(is.factor(x)) factor(x) else x)




#### 1. sample-size-based rarefaction curves  ###
#################################################

## General estimations


dl <-split(dat02[,20], dat02$community)    # we create a list of vectors of abundance for all communities

out <- iNEXT(dl, q=0, datatype="abundance") 
out$DataInfo

options(max.print=999999)
out$AsyEst

sampling.accuracy <- out$AsyEst$Observed-out$AsyEst$Estimator


## Examples urban assemblages to draw figures

dat.Amrav <- subset(dat02, city=="Amravati")
dat.Amrav[] <- lapply(dat.Amrav, function(x) if(is.factor(x)) factor(x) else x)


dl.Amrav <-split(dat.Amrav[,20], dat.Amrav$community)    # we create a list of vectors of abundance for all communities

out.Amrav <- iNEXT(dl.Amrav, q=0, datatype="abundance") 
rownames(out.Amrav$DataInfo) <- c("Wildland", "Rural", "Suburban", "Urban")
out.urb$DataInfo
out.Amrav

ggiNEXT(out.Amrav, type=1, facet.var="none") # sample-size-based rarefaction curve (type = 1)



## Test to assess differences in accuracy between habitats

comm.inf <-read.table(paste0(workingData,"/Community information.txt"), header=TRUE)  # Metrics previously estimated in the section above

levels(comm.inf$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(comm.inf$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
comm.inf <- cbind(comm.inf,habitat.ordered)

Spp.richness.accuracy <- log(comm.inf$Estimated.Spp.Richness)-log(comm.inf$Observed.Spp.Richness)
Simpson.accuracy <- log(comm.inf$Estimated.Simpson)-log(comm.inf$Observed.Simpson)

comm.inf <- cbind(comm.inf,Spp.richness.accuracy,Simpson.accuracy)
comm.inf <- na.omit(comm.inf)

Spp.richness.bias = lme(Spp.richness.accuracy ~ habitat.ordered, random = ~ 1|country/city, data = comm.inf, method="ML")
summary(Spp.richness.bias)  
testInteractions(Spp.richness.bias, pairwise="habitat.ordered")
plot(Spp.richness.bias)
qqnorm(Spp.richness.bias)
Spp.richness.bias.I <- interactionMeans(Spp.richness.bias) # effect plots
plot(Spp.richness.bias.I, errorbar="ci95")


Simpson.bias = lme(Simpson.accuracy ~ habitat.ordered, random = ~ 1|country/city, data = comm.inf, method="ML")
summary(Simpson.bias)  
testInteractions(Simpson.bias, pairwise="habitat.ordered")
plot(Simpson.bias)
qqnorm(Simpson.bias)
Simpson.bias.I <- interactionMeans(Simpson.bias) # effect plots
plot(Simpson.bias.I, errorbar="ci95")

a <- ggplot(Spp.richness.bias.I, aes(x= habitat.ordered, y=Spp.richness.bias.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Spp.richness.bias.I[,2]-Spp.richness.bias.I[,3], ymax=Spp.richness.bias.I[,2]+Spp.richness.bias.I[,3]), width=.2) +
  geom_point(data=Spp.richness.bias.I, mapping=aes(x=habitat.ordered, y=Spp.richness.bias.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Species richness accuracy", cex=16) +
  geom_text(aes(label= c("a","a","b","a","a")))

b <- ggplot(Simpson.bias.I, aes(x= habitat.ordered, y=Simpson.bias.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=Simpson.bias.I[,2]-Simpson.bias.I[,3], ymax=Simpson.bias.I[,2]+Simpson.bias.I[,3]), width=.2) +
  geom_point(data=Simpson.bias.I, mapping=aes(x=habitat.ordered, y=Simpson.bias.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Simpson's index accuracy", cex=16) +
  geom_text(aes(label= c("ab","a","ab","ab","b")))


tiff(paste0(GoogleFigs,"/plot_sampling_accuracy_across_habitats.tiff"), width = 9, height = 11, units = 'in', res = 200)
ggplot2.multiplot(a,b, cols=1)
dev.off()




#### 2. Rarefied estimations of species richness  ###
#####################################################

{

ctree1 <- read.nexus(paste0(workingData,"/AllBirdsEricson1_summary.tre"))    # This is Ericson concensus tree
ctree2 <- read.nexus(paste0(workingData,"/AllBirdsHackett1_summary.tre"))    # This is Hackett concensus tree


dl<-split(dat02, dat02$city)    # we create a list with all communities of each region
ml <- lapply(dl, function(x) xtabs(Individuals~community + animal, droplevels(x)))

PDrarefied = list()
SRrarefied = list()


for (i in 1:length(ml)){
  combined <- match.phylo.comm(ctree2, ml[[i]])
  ctree.i <- combined$phy
  comm.i <- combined$comm
  PDrarefied[[i]] <- alpha(comm.i, ctree.i, raref = 1, runs = 99) 
  SRrarefied[[i]] <- alpha(comm.i, raref = 1, runs = 99) 
}

names(PDrarefied)<- names(ml)
PDrarefied

names(SRrarefied)<- names(ml)
SRrarefied

SRraref <- as.data.frame(do.call(rbind, SRrarefied), h=TRUE) 
PDraref <- as.data.frame(do.call(rbind, PDrarefied), h=TRUE) 

### Outputs

raref <- as.data.frame(cbind(rownames(SRraref), SRraref[,1], PDraref[,1]), h=TRUE)
colnames(raref) <- c("community", "SRraref", "PDraref")


library(plyr)

# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
cdata <- ddply(dat, c("country", "city", "community", "habitat", "used.urban.nonurban"), summarise,
               regional.spp.pool  = length(relative.abundance))

write.table(merge(raref,cdata, by="community"), paste0(workingData,"/Rarefied.communities.txt"))





}



#### Tests for differences across habitats ###
##############################################

{
dat<-read.table(paste0(workingData,"/Rarefied.communities.txt"), header=TRUE)
levels(dat$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(dat$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
dat <- cbind(dat,habitat.ordered)


# Species richness rarefied

spp.richness = lme(SRraref ~ habitat.ordered, random = ~ 1|country/city, data = dat, method="ML")
summary(spp.richness)
testInteractions(spp.richness, pairwise="habitat.ordered")
plot(spp.richness)
qqnorm(spp.richness)
spp.richness.I <- interactionMeans(spp.richness) # effect plots
plot(spp.richness.I, errorbar="ci95")


# Phylogenetic diversity

PDraref = lme(PDraref ~ habitat.ordered, random = ~ 1|country/city, data = dat, method="ML")
summary(PDraref)
testInteractions(PDraref, pairwise="habitat.ordered")
plot(PDraref)
qqnorm(PDraref)
PDraref.I <- interactionMeans(PDraref) # effect plots
plot(PDraref.I, errorbar="ci95")




# Final figure (effect +/- standard error)

a <- ggplot(spp.richness.I, aes(x= habitat.ordered, y=spp.richness.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=spp.richness.I[,2]-spp.richness.I[,3], ymax=spp.richness.I[,2]+spp.richness.I[,3]), width=.2) +
  geom_point(data=spp.richness.I, mapping=aes(x=habitat.ordered, y=spp.richness.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Species richness rarefied", cex=16) +
  geom_text(aes(label= c("a","ab","ab","bc","c")))

b <- ggplot(PDraref.I, aes(x= habitat.ordered, y=PDraref.I[,2])) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_errorbar(aes(ymin=PDraref.I[,2]-PDraref.I[,3], ymax=PDraref.I[,2]+PDraref.I[,3]), width=.2) +
  geom_point(data=PDraref.I, mapping=aes(x=habitat.ordered, y=PDraref.I[,2]), size=8, shape=21, fill="white") +
  labs(x = "", y = "Phylogenetic diversity rarefied", cex=16) +
  geom_text(aes(label= c("a","a","ab","b","c")))


tiff(paste0(GoogleFigs,"/plot_taxonomic_diversity_rarefied.tiff"), width = 9, height = 11, units = 'in', res = 200)
ggplot2.multiplot(a,b, cols=1)
dev.off()





}




# Quadratic entropy

## Community data preparation

comm <- acast(na.omit(dat02[,c(5,8,20)]), community~animal, value.var="Individuals", fun.aggregate=mean)   # we transform it to a matrix species*community    # fun.aggregate=mean is used as otherwise it gives the length
# comm <- acast(na.omit(dat02[,c(5,8,15)]), community~animal, value.var="occurrence", fun.aggregate=mean)   # if we use presence/absence instead of abundances
comm[comm=="NaN"] <- 0

func<-read.table(paste0(workingData,"/morphological.axes.txt"),header=TRUE)
names(func)
funcdat<-func[,c(4,6,7,15, 8:14, 16)]
names(funcdat)
rownames(funcdat) <- func$animal
funcdat <- funcdat[order(rownames(funcdat)),] # we need to order the species
funcdat <- funcdat[labels(comm[1,]),]
all <-funcdat[,5:12]  # dataset for all traits when analysed separately
rownames(all) <- rownames(funcdat)
distallmorphology <- distance(all, "euclidean")   # all 8 traits
distallmorphology <- distallmorphology/max(distallmorphology)  # we standardize to 0-1 range

rare.Rao = list()


for (i in 1:length(ml)){
  rare.Rao[[i]] <- rare_Rao(comm.i, distallmorphology, sim = TRUE, resampling = 3)
}

rare.Rao[[i]] <- rare_Rao(ml[1], distallmorphology, sim = TRUE, resampling = 3)






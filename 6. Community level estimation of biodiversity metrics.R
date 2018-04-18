###########################################################################
#####  Estimatation of functional alpha diversity across communities ######
###########################################################################

# warnings: check the two species that are lacking for diet and foraging behaviour and the NAs of one species

### GOALS: ### 

# Estimate for each community:

# Quadratic Entropy (QE) for all morphological and behavioral traits + phylogeny
# Uniqueness* for all morphological and behavioral traits + phylogeny
# Redundancies (CR = 1-Uniqueness) for all morphological and behavioral traits + phylogeny
# Simpson index (QE taxonomy or HGS)
# Species richness
# The meanD
# The balance factor

# *A community containning species that are functionally different will achieve a high Uniqueness value, and this will increase if the relative abundance of the species is even (the absolute abundance has no effect)
# *conversely, the same community will exhibit low redundancies

## We will make the estimations for

# all species
# all native species (excluding exotics)


### INPUTS: ###

# Full species*community matrix (comm) of relative abundances
# 9 Morphological traits
# Morphological axes derived from the traits defining body size, beak shape, locmotory (tarsus) shape and wing shape
# 33 foraging behavioural traits
# 8 diet categories
# Two full phylogenies


### OUTPUT: ### 

# Data.frame containning all metrics ("Morphological diversity metrics for communities.txt")
# Data.frame containning all metrics for natives only ("Morphological diversity metrics for communities natives.txt")


### ANALYSES START HERE ### 


## Community data preparation

#dat <- read.table("/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis/Urban global data April 6 2018 for R.txt", header=TRUE)
dat<-read.table(paste0(workingData,"/Urban global data April 6 2018 for R.txt"), header=TRUE)
# dat0 <- subset(dat,status=="native")  # if we want to exclude exotics
# dat0[] <- lapply(dat0, function(x) if(is.factor(x)) factor(x) else x)

dat01 <- subset(dat, urban.analyses=="yes") # if we focus on communities within urbanisation gradients
dat01[] <- lapply(dat01, function(x) if(is.factor(x)) factor(x) else x)

# dat01 <- subset(dat0,urban.analyses=="yes") # if we focus on communities within urbanisation gradients
# dat01[] <- lapply(dat01, function(x) if(is.factor(x)) factor(x) else x)

dat02 <- dat01[dat01$relative.abundance>0,] # we exclude species absent
dat02[] <- lapply(dat02, function(x) if(is.factor(x)) factor(x) else x)

comm <- acast(na.omit(dat02[,c(5,8,17)]), community~animal, value.var="relative.abundance", fun.aggregate=mean)   # we transform it to a matrix species*community    # fun.aggregate=mean is used as otherwise it gives the length
comm[comm=="NaN"] <- 0




########### Analyses for morphological data ##############
##########################################################


## Functional data preparation

#func<-read.table("/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis/morphological.axes.txt", header=TRUE)
func<-read.table(paste0(workingData,"/morphological.axes.txt"),header=TRUE)
names(func)

funcdat<-func[,c(4,6,7,15, 8:14, 16)]
names(funcdat)
rownames(funcdat) <- func$animal

funcdat <- funcdat[order(rownames(funcdat)),] # we need to order the species
funcdat <- funcdat[labels(comm[1,]),]

beakshape <- as.data.frame(funcdat$beak.shape)  # dataset for beak shape
rownames(beakshape) <- rownames(funcdat)

locomshape <- as.data.frame(funcdat$locom.shape)  # dataset for locomotory system
rownames(locomshape) <- rownames(funcdat)

bodysize <- as.data.frame(funcdat$body.size)  # dataset for body size
rownames(bodysize) <- rownames(funcdat)

hand.wing <- as.data.frame(funcdat$hand_wing)  # dataset for wing shape
rownames(hand.wing) <- rownames(funcdat)

all <-funcdat[,5:12]  # dataset for all traits when analysed separately
rownames(all) <- rownames(funcdat)


## estimation of euclidean distances among morphological traits

distallmorphology <- distance(all, "euclidean")   # all 8 traits
distallmorphology <- distallmorphology/max(distallmorphology)  # we standardize to 0-1 range

distallmorphology3PCAs <- distance(funcdat, "euclidean")   # all three PCAs
distallmorphology3PCAs <- distallmorphology3PCAs/max(distallmorphology3PCAs)  

distbeak <- distance(beakshape, "euclidean")    # beak shape
distbeak <- distbeak/max(distbeak)   

distlocom <- distance(locomshape, "euclidean")   # locomotory shape
distlocom <- distlocom/max(distlocom)

distsize <- distance(bodysize, "euclidean")    # body size PCA all 8 traits
distsize <- distsize/max(distsize)

distwinghand<- distance(hand.wing, "euclidean")   # wing hand index
distwinghand<- distwinghand/max(distwinghand)   


## estimation of phylogenetic distances among species

ctree1 <- read.nexus(paste0(workingData,"/AllBirdsEricson1_summary.tre"))    # This is Ericson concensus tree
ctree2 <- read.nexus(paste0(worlingData,"/AllBirdsHackett1_summary.tre"))    # This is Hackett concensus tree

combined <- match.phylo.comm(ctree1, comm)
ctree.Eric <- combined$phy
comm.Eric <- combined$comm

combined <- match.phylo.comm(ctree2, comm)
ctree.Hack <- combined$phy
comm.Hack <- combined$comm

phydisE <- as.dist(cophenetic(ctree.Eric))
phydisE <- phydisE/max(phydisE)

phydisH <- as.dist(cophenetic(ctree.Hack))
phydisH <- phydisH/max(phydisH)




## We will first use adiv to estimate:

  # N (species richness, used for mistakes control)
  # Q (quadratic diversity)
  # D (Simpson diversity, i.e taxonomic Q)
  # Community uniqueness: Q/D
  # Community redundancy (CR = 1-uniqueness)

all.morph <- uniqueness(comm, distallmorphology, abundance = TRUE)
PCA3 <- uniqueness(comm, distallmorphology3PCAs, abundance = TRUE)
beak <- uniqueness(comm, distbeak, abundance = TRUE)
locom <- uniqueness(comm, distlocom, abundance = TRUE)
size <- uniqueness(comm, distsize, abundance = TRUE)
winghand <- uniqueness(comm, distwinghand, abundance = TRUE)
phyE <- uniqueness(comm.Eric, dis = phydisE, abundance = TRUE)
phyH <- uniqueness(comm.Hack, dis = phydisH, abundance = TRUE)


## We will next use the function QEpart.R to estimate:

  # The meanD
  # The balance factor

all.morph.2 <- QEpart(comm, distallmorphology)
PCA3.2 <- QEpart(comm, distallmorphology3PCAs)
beak.2 <- QEpart(comm, distbeak)
locom.2 <- QEpart(comm, distlocom)
size.2 <- QEpart(comm, distsize)
winghand.2 <- QEpart(comm, distwinghand)
phyE.2 <- QEpart(comm.Eric, dis = phydisE)
phyH.2 <- QEpart(comm.Hack, dis = phydisH)



## Preparing data for subsequent analyses

FDmorphology<-as.data.frame(cbind(labels(comm[,2]),all.morph$red$N,all.morph$red$D,all.morph$red$Q,all.morph$red$U,1-all.morph$red$U,PCA3$red$Q,PCA3$red$U,1-PCA3$red$U,beak$red$Q,beak$red$U,1-beak$red$U,locom$red$Q,locom$red$U,1-locom$red$U,size$red$Q,size$red$U,1-size$red$U,winghand$red$Q,winghand$red$U,1-winghand$red$U,phyE$red$Q,phyE$red$U,1-phyE$red$U,phyH$red$Q,phyH$red$U,1-phyH$red$U,all.morph.2$meanD,PCA3.2$meanD,beak.2$meanD,locom.2$meanD,size.2$meanD,winghand.2$meanD,phyE.2$meanD,phyH.2$meanD,all.morph.2$Balance,PCA3.2$Balance,beak.2$Balance,locom.2$Balance,size.2$Balance,winghand.2$Balance,phyE.2$Balance,phyH.2$Balance))

colnames(FDmorphology)<-c("community","Species.richness","QE.taxonomy","QE.all.morph","Uniqueness.all.morph","CR.all.morph","QE.PCA3","Uniqueness.PCA3","CR.PCA3","QE.beak","Uniqueness.beak","CR.beak","QE.locom","Uniqueness.locom","CR.locom","QE.size","Uniqueness.size","CR.size","QE.winghand","Uniqueness.winghand","CR.winghand","QE.phyE","Uniqueness.phyE","CR.phyE","QE.phyH","Uniqueness.phyH","CR.phyH","all.morph.meanD","PCA3.meanD","beak.meanD","locom.meanD","size.meanD","winghand.meanD","phyE.meanD","phyH.meanD","all.morph.Balance","PCA3.Balance","beak.Balance","locom.Balance","size.Balance","winghand.Balance","phyE.Balance","phyH.Balance")


# We add habitat and study site information

tmp <- ddply(dat0, c("country", "city", "community", "habitat"), summarise,
               N    = length(relative.abundance))

tmp2 <- merge(FDmorphology,tmp[,-5], by="community")
      
write.table(tmp2, paste0(workingData,"/Morphological diversity metrics for communities.txt"))
# write.table(tmp2, "/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis/Morphological diversity metrics for communities natives.txt")







## Not yet finished, waiting for Joe's data

################ Analyses for diet  ######################
##########################################################


### Functional diet data

setwd("~/ownCloud2/Science/Research/Urbanisation/Functional diversity and urbanization/Data and analyses")

# Missing data Thalasseus_eurygnatha, Saxicola_torquatus

diet <-read.table("Diet urban birds 8 April 2018 for R.txt", header=TRUE)
names(diet)
diet <- diet[order(diet$Scientific),] # wee need to order the species

# write.table(dietdat,"borrar.txt")

dietdat<-na.omit(diet[,c(7:17)])
names(dietdat)
rownames(dietdat) <- dietdat$Scientific
dietdat <- dietdat[-1127,-1]

dietdat <- dietdat[labels(comm[1,]),]  # we take only species in communities

dietdat1 <-  read.table("borrra1.txt", h=TRUE)
distdiet1 <- gowdis(dietdat1, ord="classic")
#distdiet <- distance(dietdat, "gower")    # all traits
distdiet2 <- gowdis(dietdat1[1:5], ord="classic")
distdiet3 <- gowdis(dietdat1[6:9], ord="classic")


dat <- read.table("Urban global data April 6 2018 for R.txt", header=TRUE)
# dat0 <- subset(dat,status=="native")  # if we want to exclude exotics
dat0 <- subset(dat,urban.analyses=="yes") # if we focus on communities within urbanisation gradients
# dat0 <- subset(dat,location=="Amravati" | location=="BCN.B")

dat01 <- dat0[dat0$animal!="Saxicola_torquatus",]

dat01 <- dat01[dat01$relative.abundance>0,] # we exclude species absent
comm01 <- acast(na.omit(dat01[,c(5,8,17)]), community~animal, value.var="relative.abundance", fun.aggregate=mean)   # we transform it to a matrix species*community    # fun.aggregate=mean is used as otherwise it gives the length
comm01[comm01=="NaN"] <- 0

Rao.diet <- QE(comm01, distdiet1)
Rao.diet.animal <- QE(comm01, distdiet2)
Rao.diet.vegetal <- QE(comm01, distdiet3)



# Prepare the data for subsequent analyses

### Outputs

Rao <- cbind(labels(comm[,2]), Rao.all, Rao.3PCAs, Rao.beak, Rao.locom, Rao.size, Rao.winghand, Rao.spp, Rao.phy, Rao.diet, Rao.diet.animal, Rao.diet.vegetal)
colnames(Rao) <- c("community", "QE.all", "QE.3PCAs", "QE.beak","QE.locom","QE.size", "QE.winghand", "QE.spp", "QE.phyl", "QE.diet", "QE.diet.anim", "QE.diet.veg")




library(plyr)

# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
cdata <- ddply(dat0, c("country", "city", "community", "habitat"), summarise,
               N    = length(relative.abundance))

write.table(merge(Rao,cdata[,-5], by="community"), "QE.morphol.comm.txt")



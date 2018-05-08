#########################################################################
#####  Estimation of functional alpha diversity across communities ######
#########################################################################

# warnings: check the two species that are lacking for diet and foraging behaviour and the NAs of one species
# warnings: adiv has some incopatibilities with other packages; if it does not work, start again, upload the package and use only the packages subsequently requested 

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


## Community data preparation (if you run this part, you get communities*species abundances of natives)

{

dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)
dat$community <- factor(dat$community)
dat <- subset(dat,status=="native")  # if we want to exclude exotics
dat[] <- lapply(dat, function(x) if(is.factor(x)) factor(x) else x)

dat01 <- subset(dat, urban.analyses=="yes") # if we focus on communities within cities or along urbanisation gradients
dat01[] <- lapply(dat01, function(x) if(is.factor(x)) factor(x) else x)

dat02 <- dat01[dat01$relative.abundance>0,] # we exclude species absent
# dat02 <- dat01  # if we want to include all occurrence data
dat02[] <- lapply(dat02, function(x) if(is.factor(x)) factor(x) else x)


comm <- acast(na.omit(dat02[,c(5,8,17)]), community~animal, value.var="relative.abundance", fun.aggregate=mean)   # we transform it to a matrix species*community    # fun.aggregate=mean is used as otherwise it gives the length
# comm <- acast(na.omit(dat02[,c(5,8,15)]), community~animal, value.var="occurrence", fun.aggregate=mean)   # if we use presence/absence instead of abundances
comm[comm=="NaN"] <- 0


}

########### Analyses for morphological data ##############
##########################################################


## Functional data preparation

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
ctree2 <- read.nexus(paste0(workingData,"/AllBirdsHackett1_summary.tre"))    # This is Hackett concensus tree

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

all.morph <- uniqueness(comm, distallmorphology)
PCA3 <- uniqueness(comm, distallmorphology3PCAs)
beak <- uniqueness(comm, distbeak)
locom <- uniqueness(comm, distlocom)
size <- uniqueness(comm, distsize)
winghand <- uniqueness(comm, distwinghand)
phyE <- uniqueness(comm.Eric, dis = phydisE)
phyH <- uniqueness(comm.Hack, dis = phydisH)


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

tmp <- ddply(dat, c("country", "city", "community", "habitat", "used.urban.nonurban"), summarise,
               Regional.spp.richness = length(relative.abundance))

tmp2 <- merge(FDmorphology,tmp, by="community")
      
write.table(tmp2, paste0(workingData,"/Morphological diversity metrics for communities.txt"))
# write.table(tmp2, paste0(workingData,"/Morphological diversity metrics for communities natives.txt"))
# write.table(tmp2, paste0(workingData,"/Morphological diversity metrics for communities ocurrences.txt"))
# write.table(tmp2, paste0(workingData,"/Morphological diversity metrics for communities ocurrences natives.txt"))




################ Analyses for diet  ######################
##########################################################


### Functional diet data

## Need to estimate first "comm" running the firts part of this code

diet <-read.table(paste0(workingData,"/Diet urban birds 28 April 2018 for R.txt"), header=TRUE)
# Missing data Thalasseus_eurygnatha, Saxicola_torquatus, already removed from diet

names(diet)
diet <- diet[order(diet$animal),] # wee need to order the species
diet<-diet[,c(5,6:16)]
names(diet)
rownames(diet) <- diet$animal
diet <- diet[,-c(1:2)]


## Estimation of distances

dietdat <- na.omit(diet[labels(comm[1,]),])  # we take only species in communities

comm1 <- comm[,rownames(dietdat)]  # as one species is missing, we need to update community

distdiet<- distance(dietdat, "euclidean") 
distdiet <- distdiet/max(distdiet) 

all.diet <- QEpart(comm1, distdiet)

# redundancies are estimated as 1-QE/Simpson (Ricotta et al. 2016)

redundancy <- 1-(all.diet$QE/all.diet$Simpson)

## Preparing data for subsequent analyses

FDdiet <-as.data.frame(cbind(labels(comm1[,2]), all.diet$QE, all.diet$meanD, all.diet$Balance, redundancy))

colnames(FDdiet)<-c("community","QE.diet","diet.meanD","diet.Balance","CR.diet")



# We add habitat and study site information

tmp <- ddply(dat, c("country", "city", "community", "habitat", "used.urban.nonurban"), summarise,
             Regional.spp.richness = length(relative.abundance))

tmp2 <- merge(FDdiet,tmp, by="community")

# write.table(tmp2, paste0(workingData,"/Diet diversity metrics for communities.txt"))
write.table(tmp2, paste0(workingData,"/Diet diversity metrics for communities natives.txt"))
# write.table(tmp2, paste0(workingData,"/Diet diversity metrics for communities ocurrences.txt"))
# write.table(tmp2, paste0(workingData,"/Diet diversity metrics for communities ocurrences natives.txt"))






################ Analyses for morphology by diet for natives ####################
#################################################################################


### We will examine FD for all traits, subsetted by insectivory, nectarivory, ...

## We need to estimate first comm for natives running the beginning of the code


diet <-read.table(paste0(workingData,"/Diet urban birds 28 April 2018 for R.txt"), header=TRUE)
names(diet)
diet<-diet[,c(5,6:16)]

func<-read.table(paste0(workingData,"/morphological.axes.txt"),header=TRUE)
names(func)
all <-func[,c(1,9:16)]  # dataset for all traits when analysed separately

diet.morph <- merge(diet,all, by="animal")

diet.morph <- diet.morph[order(diet.morph$animal),] # wee need to order the species
names(diet.morph)
rownames(diet.morph) <- diet.morph$animal
diet.morph <- diet.morph[,-c(1,2)]

diet.morph <- na.omit(diet.morph[labels(comm[1,]),])  # we take only species in communities


## Insectivorous, insects > 30% in diet

insect.50 <- diet.morph[diet.morph$Diet.Inv>50,]
comm.insect.50 <- comm[,rownames(insect.50)]

distinsect.50 <- distance(insect.50[,c(11:14,15,16,17)], "euclidean")   # all 8 traits



morph.insect.50 <- QEpart(comm.insect.50, distinsect.50)

# redundancies are estimated as 1-QE/Simpson (Ricotta et al. 2016)

morph.insect.50.redundancy <- 1-(morph.insect.50$QE/morph.insect.50$Simpson)



## Preparing data for subsequent analyses

FD.morphol.diet <-as.data.frame(cbind(labels(comm.insect.50[,2]), morph.insect.50$QE, morph.insect.50$meanD, morph.insect.50$Balance, morph.insect.50.redundancy))

colnames(FD.morphol.diet)<-c("community","QE.insectiv","insectiv.meanD","insectiv.Balance","CR.insectiv")



# We add habitat and study site information

tmp <- ddply(dat, c("country", "city", "community", "habitat", "used.urban.nonurban"), summarise,
             Regional.spp.richness = length(relative.abundance))

tmp2 <- merge(FD.morphol.diet,tmp, by="community")

write.table(tmp2, paste0(workingData,"/Morphology-Diet diversity metrics for communities natives.txt"))




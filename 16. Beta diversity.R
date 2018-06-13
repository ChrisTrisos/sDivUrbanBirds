#########################################################################
#####  Estimation of functional beta diversity across communities ######
#########################################################################

# warnings: check the two species that are lacking for diet and foraging behaviour and the NAs of one species
# warnings: adiv has some incopatibilities with other packages; if it does not work, start again, upload the package and use only the packages subsequently requested 

### GOALS: ### 



####################
### All habitats ### 
####################

## Community data preparation (if you run this part, you get communities*species abundances of natives)

{
  
  dat<-read.table(paste0(workingData,"/Urban global data May 14 2018 for R.txt"), header=TRUE)
  dat$community <- factor(dat$community)
  dat <- subset(dat,status=="native")  # if we want to exclude exotics
  dat[] <- lapply(dat, function(x) if(is.factor(x)) factor(x) else x)
  
  dat01 <- subset(dat, used.urban.nonurban=="yes") # if we focus on communities within cities or along urbanisation gradients
  dat01[] <- lapply(dat01, function(x) if(is.factor(x)) factor(x) else x)
  
  dat02 <- dat01[dat01$relative.abundance>0,] # we exclude species absent
  # dat02 <- dat01  # if we want to include all occurrence data
  dat02[] <- lapply(dat02, function(x) if(is.factor(x)) factor(x) else x)
  
  
  comm <- acast(na.omit(dat02[,c(5,8,17)]), community~animal, value.var="relative.abundance", fun.aggregate=mean)   # we transform it to a matrix species*community    # fun.aggregate=mean is used as otherwise it gives the length
  # comm <- acast(na.omit(dat02[,c(5,8,15)]), community~animal, value.var="occurrence", fun.aggregate=mean)   # if we use presence/absence instead of abundances
  comm[comm=="NaN"] <- 0
  
  
}


## Functional data preparation

{
func<-read.table(paste0(workingData,"/morphological.axes.txt"),header=TRUE)
names(func)
rownames(func) <- func$animal

funcdat<-func[levels(dat02$animal),c(4,6,7,15, 8:14, 16)]
funcdat[] <- lapply(funcdat, function(x) if(is.factor(x)) factor(x) else x)

funcdat <- funcdat[order(rownames(funcdat)),] # we need to order the species

beakshape <- as.data.frame(funcdat$beak.shape)  # dataset for beak shape
rownames(beakshape) <- rownames(funcdat)

locomshape <- as.data.frame(funcdat$locom.shape)  # dataset for locomotory system
rownames(locomshape) <- rownames(funcdat)

bodysize <- as.data.frame(funcdat$body.size)  # dataset for body size
rownames(bodysize) <- rownames(funcdat)

hand.wing <- as.data.frame(funcdat$hand_wing)  # dataset for wing shape
rownames(hand.wing) <- rownames(funcdat)

all <-funcdat[,5:12]  # dataset for all traits when analysed separately



## estimation of euclidean distances among morphological traits

distallmorphology <- distance(na.omit(all), "euclidean")   # all 8 traits
distallmorphology <- distallmorphology/max(distallmorphology)  # we standardize to 0-1 range

distbeak <- distance(beakshape, "euclidean")    # beak shape
distbeak <- distbeak/max(distbeak)   

distlocom <- distance(locomshape, "euclidean")   # locomotory shape
distlocom <- distlocom/max(distlocom)

distsize <- distance(bodysize, "euclidean")    # body size PCA all 8 traits
distsize <- distsize/max(distsize)

distwinghand<- distance(hand.wing, "euclidean")   # wing hand index
distwinghand<- distwinghand/max(distwinghand)   

}


## Functional beta diversity with BetQmult


########################################################################################################################################
# betaQmult: R function to compute functional beta-diversity based on the multiplicative decomposition of the Rao's quadratic entropy  # 
# Four contrasted examples are provided at the end of the script
# Author: Sébastien Villéger, sebastien.villeger@univ-tlse3.fr
#
# INPUTS:
#   - "abundances": matrix (C x S) of abundances of the S species in the C communities
#   - "functdist": matrix (S x S) or dist object with pairwise functional distances between the S species
#   - species names in "abundances" and in "functdist" must be identical
#      NA are automatically replaced by 0 in 'abundances'
#      NA are not allowed in 'functdist'
#
# OUTPUTS:  a list of 7 elements
#   - "nbc" : number of communities
#   - "nbsp" (vector, length C) : number of species in each community
#   - "weight" (vector, length C) : contribution of each community to the regional abundances
#   - "eqmQalpha" : equivalent number of species of the abundance-weighted mean of local Rao's quadratic entropy values(Qalpha)
#   - "eqQgamma" : equivalent number of species of the regional Rao's quadratic entropy 
#   - "Qbeta" : eqQgamma/eqmQalpha, raw value of functional beta-diversity
#   - "Qbetast" : standardized functional beta-diversity, =(Qbeta-1)/(C-1)
#####################################################################################################################################


## We estimate beta diversity with betaQmult for each each city

dl<-split(dat02, dat02$city)    # we create a list with all communities of each region
ml <- lapply(dl, function(x) xtabs(effort.corrected.abundance~community + animal, droplevels(x)))

Beta_FD = list()

dist <- as.matrix(distallmorphology)

for (i in 1:length(ml)){
  Xi <- labels(ml[[i]])[2]$animal
  disti <- dist[Xi,Xi]
  Beta_FD[[i]] <- betaQmult(disti,ml[[i]]) 
  
  Beta_FD
}

names(Beta_FD)<- names(ml) 

res <- do.call(rbind, Beta_FD)
resfinal <- res[,-c(2,3)]



## We estimate beta diversity with betaQmult for each pair of communities within each city


dl<-split(dat02, dat02$city)    # we create a list with all communities of each region
ml <- lapply(dl, function(x) xtabs(effort.corrected.abundance~community + animal, droplevels(x)))  # we use effort corrected abundances


Beta_FD = list()

dist <- as.matrix(distallmorphology)

for (i in 1:length(ml)){
  Xi <- labels(ml[[i]])[2]$animal
  disti <- dist[Xi,Xi]
  
  ml[[i]]
  comb<- combn(rownames(ml[[i]]), 2)
  resCity<- data.frame(city=names(ml)[i], comm1=comb[1,], comm2=comb[2,], Qbetast=numeric(ncol(comb)), Qbeta=numeric(ncol(comb)))
  for (j in 1:ncol(comb)){
    tmp<- ml[[i]][c(comb[1, j], comb[2, j]), ]
    resCity[j, c("Qbetast", "Qbeta")]<- as.numeric(betaQmult(disti, tmp)[c("Qbetast", "Qbeta")])
    
  }
  
  Beta_FD[[i]] <- resCity
}
names(Beta_FD)<- names(ml)

res <- do.call(rbind, Beta_FD)


levels(dat02$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
habitat.ordered  = factor(dat02$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
dat03 <- cbind(dat02,habitat.ordered)

comnames <- unique(dat03[,c(5,5,7,7)])
colnames(comnames) <- c("comm1","comm2","hab.comm1", "hab.comm2")

tmp3 <- merge(res,comnames[,c(1,3)], by="comm1")
beta.div.comm <- merge(tmp3,comnames[,c(2,4)], by="comm2")
beta.div.comm[ is.na(beta.div.comm) ] <- NA

beta.div.comm <- na.omit(transform(beta.div.comm, hab.comp=as.factor(paste(hab.comm1, hab.comm2, sep="_"))))

# "Rural_Rural” , “Rural_Suburban” , “Rural_Urban” , “Rural_Urban_Park"   ,"Rural_Wildland” , “Suburban_Rural” , “Suburban_Suburban” , “Suburban_Urban"       ,"Suburban_Urban_Park” , “Suburban_Wildland” , “Urban_Park_Rural” , “Urban_Park_Suburban" , "Urban_Park_Urban” , “Urban_Park_Urban_Park” , “Urban_Park_Wildland” , “Urban_Rural" , "Urban_Suburban” , “Urban_Urban” , “Urban_Urban_Park” , “Urban_Wildland" , "Wildland_Rural” , “Wildland_Suburban” , “Wildland_Urban” , “Wildland_Urban_Park",  "Wildland_Wildland"

levels(beta.div.comm$hab.comp) <- c("Rural_Rural" , "Suburban_Rural" , "Urban_Rural", "Urban_Park_Rural","Rural_Wildland", "Suburban_Rural", "Suburban_Suburban", "Urban_Suburban", "Suburban_Urban_Park", "Suburban_Wildland", "Urban_Park_Rural", "Suburban_Urban_Park", "Urban_Urban_Park", "Urban_Park_Urban_Park", "Urban_Park_Wildland", "Urban_Rural", "Urban_Suburban", "Urban_Urban", "Urban_Urban_Park", "Urban_Wildland", "Rural_Wildland", "Suburban_Wildland", "Urban_Wildland", "Urban_Park_Wildland",  "Wildland_Wildland")

# write.table(beta.div.comm, paste0(workingData,"/Rao beta diversity.txt"))






































## Estimation of beta diversity

# we need to create variable to define the cities and habitats within cities for each community, including only the species we have in the communities

stru <- dat02[,c(2,5,7)]
stru[] <- lapply(stru, function(x) if(is.factor(x)) factor(x) else x)
stru <- na.omit(unique(stru[,1:3]))

stru <- stru[order(stru[,2]),] # we need to order the communities

levels(stru$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")

stru$new_col <- do.call(paste, c(stru[c("habitat", "city")], sep=" "))
stru1 <- as.data.frame(stru[,4])
rownames(stru1) <- stru[,2]

EqRao(comm, , stru1)
r2_GS <- rtestEqRao(comm, , stru2, level=1, nrep=9, option="normed2")  # level 1 tests for dissimilarities among sites within regions

# EqRao(comm, distallmorphology, stru1)
r2_GS <- rtestEqRao(comm, distallmorphology, stru2, level=1, nrep=3, option="normed2") 



## We repeat the analysis with betaQmult

dl<-split(dat02, dat02$city)    # we create a list with all communities of each region
ml <- lapply(dl, function(x) xtabs(abundance~community + animal, droplevels(x)))

Beta_FD = list()

dist <- as.matrix(distallmorphology)

for (i in 1:length(ml)){
  Xi <- labels(ml[[i]])[2]$animal
  disti <- dist[Xi,Xi]
  Beta_FD[[i]] <- betaQmult(disti,ml[[i]]) 

  Beta_FD
}

names(Beta_FD)<- names(ml) 



### Analyses with betapart ###
##############################


# Based on convex hulls

dl<-split(dat02, dat02$city)    # we create a list with all communities of each region
ml <- lapply(dl, function(x) xtabs(abundance~community + animal, droplevels(x)))

Beta_FD = list()

comm.test <- ml[[1]]
traits.test <- all[labels(ml[[1]])[[2]],]


test.pair<-functional.beta.pair(x=comm.test, traits=traits.test, index.family = "jaccard" )



























############################################# 
###############  not used ###################
############################################# 
  
  
  
  
##################################
### Urban vs wildland habitats ### 
##################################

## Community data preparation (if you run this part, you get communities*species abundances of natives)

{
  
  dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)
  levels(dat$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  
  dat <- subset(dat,status=="native" & habitat=="Urban" | habitat=="Wildland")  # we want to exclude exotics anf focus on the two habitats
    dat <- subset(dat,city!="Toronto" & city!="WASH" & city!="Santiago" & city!="Saskatoon"& city!="Rennes" & city!="Orebro"& city!="Montreal" & city!="Mindanao"& city!="Jerusalem" & city!="Bristol" & city!="Jamshedpur" & city!="Phoenix" & city!="Phalaborwa" & city!="Nanchong" & city!="Brisbane")# we also exclude cities wichi lack any of the two habitats

  dat[] <- lapply(dat, function(x) if(is.factor(x)) factor(x) else x)
  
  
  dat01 <- subset(dat, used.urban.nonurban=="yes") # if we focus on communities within cities or along urbanisation gradients
  dat01[] <- lapply(dat01, function(x) if(is.factor(x)) factor(x) else x)
  
  dat02 <- dat01[dat01$relative.abundance>0,] # we exclude species absent
  # dat02 <- dat01  # if we want to include all occurrence data
  dat02[] <- lapply(dat02, function(x) if(is.factor(x)) factor(x) else x)
  
  comm <- acast(na.omit(dat02[,c(5,8,17)]), community~animal, value.var="relative.abundance", fun.aggregate=mean)   # we transform it to a matrix species*community    # fun.aggregate=mean is used as otherwise it gives the length
  # comm <- acast(na.omit(dat02[,c(5,8,15)]), community~animal, value.var="occurrence", fun.aggregate=mean)   # if we use presence/absence instead of abundances
  comm[comm=="NaN"] <- 0
  
}


## Functional data preparation

func<-read.table(paste0(workingData,"/morphological.axes.txt"),header=TRUE)
names(func)
rownames(func) <- func$animal

funcdat<-func[levels(dat02$animal),c(4,6,7,15, 8:14, 16)]
funcdat[] <- lapply(funcdat, function(x) if(is.factor(x)) factor(x) else x)

funcdat <- funcdat[order(rownames(funcdat)),] # we need to order the species

beakshape <- as.data.frame(funcdat$beak.shape)  # dataset for beak shape
rownames(beakshape) <- rownames(funcdat)

locomshape <- as.data.frame(funcdat$locom.shape)  # dataset for locomotory system
rownames(locomshape) <- rownames(funcdat)

bodysize <- as.data.frame(funcdat$body.size)  # dataset for body size
rownames(bodysize) <- rownames(funcdat)

hand.wing <- as.data.frame(funcdat$hand_wing)  # dataset for wing shape
rownames(hand.wing) <- rownames(funcdat)

all <-funcdat[,5:12]  # dataset for all traits when analysed separately



## estimation of euclidean distances among morphological traits

distallmorphology <- distance(na.omit(all), "euclidean")   # all 8 traits
distallmorphology <- distallmorphology/max(distallmorphology)  # we standardize to 0-1 range

distbeak <- distance(beakshape, "euclidean")    # beak shape
distbeak <- distbeak/max(distbeak)   

distlocom <- distance(locomshape, "euclidean")   # locomotory shape
distlocom <- distlocom/max(distlocom)

distsize <- distance(bodysize, "euclidean")    # body size PCA all 8 traits
distsize <- distsize/max(distsize)

distwinghand<- distance(hand.wing, "euclidean")   # wing hand index
distwinghand<- distwinghand/max(distwinghand)   



## Now we want to estimate beta diversity for each pair of communities within a city

# all urban and wildland communities to estimate gamma diversity

dl<-split(dat02, dat02$city)    # we create a list with all communities of each region
ml <- lapply(dl, function(x) xtabs(abundance~community + animal, droplevels(x)))

Beta_FD = list()

dist <- as.matrix(distallmorphology)

for (i in 1:length(ml)){
  Xi <- labels(ml[[i]])[2]$animal
  disti <- dist[Xi,Xi]
  Beta_FD[[i]] <- betaQmult(disti,ml[[i]]) 
  
  #names(Beta_FD)<- names(ml) 
  Beta_FD
}


# urban and wildland separated to estimate gamma diversity for each habitat 

dat02$new_col <- do.call(paste, c(dat02[c("habitat", "city")], sep=" "))

dl<-split(dat02, dat02$new_col)    # we create a list with all communities of each region
ml <- lapply(dl, function(x) xtabs(abundance~community + animal, droplevels(x)))

Beta_FD = list()

dist <- as.matrix(distallmorphology)

for (i in 1:length(ml)){
  Xi <- labels(ml[[i]])[2]$animal
  disti <- dist[Xi,Xi]
  Beta_FD[[i]] <- betaQmult(disti,ml[[i]]) 
  
  names(Beta_FD)<- names(ml) 
  Beta_FD
}






Beta_FD = list()

dist <- as.matrix(distallmorphology)

for (i in 1:length(ml)){
  Xi <- labels(ml[[i]])[2]$animal
  disti <- dist[Xi,Xi]
  
  ml[[i]]
  comb<- combn(rownames(ml[[i]]), 2)
  resCity<- data.frame(city=names(ml)[i], comm1=comb[1,], comm2=comb[2,], Qbetast=numeric(ncol(comb)), Qbeta=numeric(ncol(comb)))
  for (j in 1:ncol(comb)){
    tmp<- ml[[i]][c(comb[1, j], comb[2, j]), ]
    resCity[j, c("Qbetast", "Qbeta")]<- as.numeric(betaQmult(disti, tmp)[c("Qbetast", "Qbeta")])
    
  }
  
  Beta_FD[[i]] <- resCity
}
names(Beta_FD)<- names(ml)

do.call(rbind, Beta_FD)





















# Imagine that structures is a data frame with sites as rows and only one column representing how sites are distributed among regions.  
# two levels of beta diversity: 
# beta1 diversity represents the (functional or phylogenetic) dissimilarities among sites within regions; 
# beta2 diversity represents the (functional or phylogenetic) dissimilarities among regions. If level = 1 then functions rtestEqRS, rtestEqRSintra, rtestEqRao, and rtestwapqe will test for the significance of the dissimilarities among sites within regions (beta1 diversity); in contrast, if level = 2 functions rtestEqRS, rtestEqRSintra, rtestEqRao, and rtestwapqe will test for the significance of the dissimilarities among regions (beta2 diversity). As there is only one column in parameter structures and thus only two levels of diversity, level cannot be higher than 2.


# we will use EqRao and rtestEqRao and level = 1 to test for the significance of the dissimilarities among sites within regions. 


## Estimation of beta diversity with "adiv"

# we need to create variable to define the cities and habitats within cities for each community, including only the species we have in the communities

stru <- dat02[,c(2,5,7)]
stru[] <- lapply(stru, function(x) if(is.factor(x)) factor(x) else x)
stru <- na.omit(unique(stru[,1:3]))

stru <- stru[order(stru[,2]),] # we need to order the communities

stru$new_col <- do.call(paste, c(stru[c("habitat", "city")], sep=" "))
stru1 <- as.data.frame(stru[,4])
rownames(stru1) <- stru[,2]

EqRao(comm, , stru1)
r2_GS <- rtestEqRao(comm, , stru2, level=1, nrep=9, option="normed2")  # level 1 tests for dissimilarities among sites within regions

r2_GS

EqRao(comm, distallmorphology, stru2)
# r2_GS <- rtestEqRao(comm, distallmorphology, stru1, level=1, nrep=3, option="normed2") 


## We repeat the analysis with betaQmult

dl<-split(dat02, dat02$city)    # we create a list with all communities of each region
ml <- lapply(dl, function(x) xtabs(abundance~community + animal, droplevels(x)))

Beta_FD = list()

dist <- as.matrix(distallmorphology)

for (i in 1:length(ml)){
  Xi <- labels(ml[[i]])[2]$animal
  disti <- dist[Xi,Xi]
  Beta_FD[[i]] <- betaQmult(disti,ml[[i]]) 
  
  #names(Beta_FD)<- names(ml) 
  Beta_FD
}


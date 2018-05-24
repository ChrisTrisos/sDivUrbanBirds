#################################################
#####  Modelling tolerance to urbanisation ######
#################################################


### GOALS: ###

    # 1. Use ZI Poisson GLMMS to model corrected abundances as a function of habitat to test whether abundance outside the city is a good predictor of abundance within the city
    # 2. Estimate urban tolerance as the log-log difference between corrected abundances within and outside the city
    # 3. Use Gaussian mixed models to assess whether species exhibit consistent differences in tolerance, including mean abundance outside the city as offset
    # 4. Use Gaussian mixed models to assess whether morphological niche traits are related to tolerances, including mean abundance outside the city as offset
    # 5. Investigate the correlation between morphological niche traits and response traits previously identified to define urban tolerance


### 1.  ZI Poisson GLMMS to model corrected abundances as a function of habitat ###
###################################################################################  

## Data:

# we will first use effort.corrected.abundance, which is the appropriate way to do so

dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)
# dat0 <- subset(dat,status=="native")  # if we want to exclude exotics
# dat0[] <- lapply(dat0, function(x) if(is.factor(x)) factor(x) else x)

dat01 <- subset(dat, urban.analyses=="yes") # if we focus on communities within urbanisation gradients
dat01[] <- lapply(dat01, function(x) if(is.factor(x)) factor(x) else x)

levels(dat01$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")

dat01$link <- as.factor(paste(dat01$city,dat01$animal, dat01$habitat,sep = ""))  # we create a new variable to then link with subsequent data



## We estimate mean corrected relative abundance in the wildland

tmp <- ddply(dat01, c("city", "habitat", "animal"), summarise,
             propagule  = sum(effort.corrected.abundance)/sum(community.size))
tmp$link <- as.factor(paste(tmp$city,tmp$animal, tmp$habitat,sep = ""))  # we create a new variable to then link with subsequent data

tmp <- tmp[,-c(1:3)]

dat03 <- merge(dat01[,-19],tmp, by="link")



## Morphological data:

func<-read.table(paste0(workingData,"/morphological.axes.txt"),header=TRUE)
func <- func[,c(1, 4,6,7,15, 8:14, 16)]
names(func)


## Merged urban and morphological data:

dat.toler <-merge(dat03[,-c(1,2,7,10,11,12,16,25)], func, by="animal")
dat.toler$species <- dat.toler$animal
str(dat.toler)


## Phylogenies

ctree1 <- read.nexus(paste0(workingData,"/AllBirdsEricson1_summary.tre"))    # This is Ericson concensus tree
ctree2 <- read.nexus(paste0(workingData,"/AllBirdsHackett1_summary.tre"))    # This is Hackett concensus tree


## MCMCglmm 

zi.prior <-  list(R = list(V = diag(2), n = 0.002, fix = 2),
                  G = list(G1 = list(V = 1, n = 0.002),
                           G2 = list(V = 1, n = 0.002),
                           G3 = list(V = 1, n = 0.002)))

# residual variance for the zero-inflation is identified with rcov=idh(trait):units and fixing the [2,2] element of the covariance matrix at something (e.g  1) in the prior (i.e. add fix=2 to your R prior).


## Test of traits influencing presence/abundance in cities

dat.toler.01 <- subset(dat.toler, habitat=="Urban")
dat.toler.urb <- na.omit(dat.toler.01[,c(1,2,5,12,14,18:32)])

# How many cities and species? 
length(unique(dat.toler.urb$city)) #27 cities
length(unique(dat.toler.urb$animal)) #1068 cities



m1 <- MCMCglmm(as.integer(effort.corrected.abundance) ~ trait-1 + at.level(trait,1):log(propagule+0.1),    # trait-1 rather than -1 as we do not want the intercept for the Poisson process and the zero-inflation to be the same
               random = ~idh(at.level(trait,1)):city +  # at.level separates predictors in the Poisson and Binomial model
                 idh(at.level(trait,1)):animal + idh(at.level(trait,1)):species,    #2. If you have many observations per species then I would put a non-phylogenetic species effect in too. If there are few (at the limit, only one) then it may be hard to separate the phylogenetic from the non-phylogenetic.
               rcov = ~ idh(trait):units,  #trat1 refers to the Poisson component
               prior = zi.prior,
               pedigree=ctree1,
               nitt = 2550, burnin = 50, thin = 5,
               data = dat.toler.urb, 
               family = "zipoisson", 
               verbose = TRUE, 
               pr = FALSE, 
               pl = FALSE)   #test also zapoisson
summary(m1)

summary(m1$Sol)
HPDintervals(m1$Sol)
autocorr(m1$Sol)
plot(m1$Sol)






### 2.  Estimate urban tolerance ###
####################################

#### Warning: To avoid conflicts, download plyr first !!!#### 

library(plyr) 
library(tidyverse)
library(stringr)

## Data preparation ## 

dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)

dat0 <- subset(dat, used.urban.nonurban=="yes") # if we focus on communities along urbanisation gradients
dat0[] <- lapply(dat0, function(x) if(is.factor(x)) factor(x) else x)

dat0 <- dat0[,c(8,2,3,5,7,13,16,17,18,20,23,26)]  # select key variables

dat1<- subset(dat0,dat0$hab!="garden") #Remove records made in gardens
dat1$hab<-factor(dat1$hab)

dat.abundance <- dat1[dat1$effort.corrected.abundance!="NA" & dat1$animal!="NA",] # remove cities that report only ocurrences and a few species from Pauri that were not identified
dat.abundance[] <- lapply(dat.abundance, function(x) if(is.factor(x)) factor(x) else x)


# Change the levels names 
levels(dat1$hab)
dat1$hab1<-revalue(dat1$hab, c("urb"="urban","sub"="suburban","rural"="rural","open_wild"="low-altered","closed_wild"="low-altered","little_urbanised"="suburban","pasture"="rural","plantation"="rural","rural_wild"="rural","wild_mosaic"="low-altered"))
dat.abundance$hab1<-revalue(dat.abundance$hab, c("urb"="urban","sub"="suburban","rural"="rural","open_wild"="low-altered","closed_wild"="low-altered","little_urbanised"="low-altered","pasture"="rural","plantation"="rural","rural_wild"="rural","wild_mosaic"="low-altered"))

levels(dat.abundance$hab1)



#### Data that we will use to estimate the urban tolerance index (UTI)


urban_abundance<-as_tibble(dat.abundance) # These data include cities with absolute abundance and relative abundance
length(unique(urban_abundance$animal)) # 1093 species
length(unique(urban_abundance$city)) # 28 cities




#creating urban tolerance index (UTI)

# We subset urban and non-urban abundances per species (averaging all abundances per species, habitat and city) and then estimate the log-log difference

urban_abundance<-dat.abundance%>%
  filter(status=="native", hab1=="urban")%>%
  group_by(animal, city)%>%
  summarise(
    count=n(),
    abundance.urban=mean(effort.corrected.abundance,na.rm=TRUE))

urban_abundance<-urban_abundance%>%
  unite(match.code, animal,city, sep = "_", remove = FALSE)


surrounding_abundance<-dat.abundance%>%
  filter(status=="native", hab1=="low-altered")%>%
  group_by(animal, city)%>%
  summarise(abundance.surrounding=mean(effort.corrected.abundance,na.rm=TRUE))

surrounding_abundance<-surrounding_abundance%>%
  unite(match.code, animal,city, sep = "_", remove = FALSE)

UTI_within_species <-urban_abundance%>%
  mutate(
    abundance.surrounding=surrounding_abundance$abundance.surrounding[match(match.code,surrounding_abundance$match.code)])%>%###### There are several species not found in the surrounding (closed+open) depicted as NAs, but present in the urban areas.
  mutate_if(is.numeric, funs(replace(., is.na(.), 0)))%>%     # Here we replace NA´s with 0s and then remove lines with 0s in both urban and surrounding
  filter(abundance.urban!="0" | abundance.surrounding!="0")%>%    # Eliminate cases where abundance within and outside the city are both zero
  mutate(UTI=log(abundance.urban+0.1)-log(abundance.surrounding+0.1))

UTI_within_species$species <- UTI_within_species$animal
UTI_within_species[] <- lapply(UTI_within_species, function(x) if(is.factor(x)) factor(x) else x)


write.table(UTI_within_species, paste0(workingData,"/UTI_within_species.txt"))


#Summarize UTI data at species level

UTI_species <-UTI_within_species%>%
  group_by(animal)%>%
  summarise(UTI_min=min(UTI), UTI_max=max(UTI), meanUTI=mean(UTI), n=length(UTI), abundance.surrounding=mean(abundance.surrounding), abundance.urban=mean(abundance.urban))
UTI_species[] <- lapply(UTI_species, function(x) if(is.factor(x)) factor(x) else x)

write.table(UTI_species, paste0(workingData,"/UTI_species.txt"))


### 3.  Use Gaussian mixed models to assess whether species exhibit consistent differences in tolerance ###
###########################################################################################################


UTI_within_species<-read.table(paste0(workingData,"/UTI_within_species.txt"), header=TRUE)

isTree <- drop.tip(ctree1,tip= ctree1$tip.label[!ctree1$tip.label %in% unique(UTI_within_species$animal)])  # we prune the tree
isTree$tip.label # 959 species

prior1 <- list(R=list(V=1, nu=0),
               G=list(G1=list(V=1, nu=0, alpha.mu=0, alpha.V=1e3),
                      G2=list(V=1, nu=0, alpha.mu=0, alpha.V=1e3),
                      G3=list(V=1, nu=0, alpha.mu=0, alpha.V=1e3)))

fit <- MCMCglmm(
       fixed = UTI ~ log(abundance.surrounding+0.1),
       random = ~ animal + species + city,
       data = UTI_within_species,
       family = "gaussian",
       pedigree=isTree,
       prior = prior1,
       nitt = 1010, thin = 10, burnin = 10,
       verbose = FALSE
       )




### 4.  Use Gaussian mixed models to assess whether tolerance depends on morphology ###
#######################################################################################

# We will ask whether morphology is related to UTI and abundance in the city, and whether this reflects that species with particular traits are too rare in the surroundings of the city to colonise it

# Data used

UTI_within_species<-read.table(paste0(workingData,"/UTI_within_species.txt"), header=TRUE)
UTI_species<-read.table(paste0(workingData,"/UTI_species.txt"), header=TRUE)

isTree <- drop.tip(ctree1,tip= ctree1$tip.label[!ctree1$tip.label %in% unique(UTI_within_species$animal)])  # we prune the tree
isTree$tip.label # 959 species


## Morphological data:

func<-read.table(paste0(workingData,"/morphological.axes.txt"),header=TRUE)
func <- func[,c(1, 4,6,7,15, 8:14, 16)]
names(func)


## Merged urban and morphological data:

UTI_within_species <-merge(UTI_within_species, func, by="animal")
UTI_within_species$species <- UTI_within_species$animal
str(UTI_within_species)


UTI_species <-merge(UTI_species, func, by="animal")
UTI_species$species <- UTI_species$animal
str(UTI_species)


# This firt set of models asks whether 
    # urban tolerant species tend to be morphologically different than the species from the regional pool
    # morphology is related to abundance inside and outside the city

prior1 <- list(R=list(V=1, nu=0),
               G=list(G1=list(V=1, nu=0, alpha.mu=0, alpha.V=1e3),
                      G2=list(V=1, nu=0, alpha.mu=0, alpha.V=1e3),
                      G3=list(V=1, nu=0, alpha.mu=0, alpha.V=1e3)))

tolerance.diff.morph.comm <- MCMCglmm(
  fixed = UTI ~ log(abundance.surrounding+0.1)+beak.shape+locom.shape+poly(body.size,2)+hand_wing,
  random = ~ animal + species + city,
  data = UTI_within_species,
  family = "gaussian",
  pedigree=isTree,
  prior = prior1,
  nitt = 1010, thin = 10, burnin = 10,
  verbose = FALSE
)
summary(tolerance.diff.morph.comm)


inside.abundance.morph.comm <- MCMCglmm(
  fixed = log(abundance.urban+0.1) ~ beak.shape+locom.shape+poly(body.size,2)+hand_wing,
  random = ~ animal + species + city,
  data = UTI_within_species,
  family = "gaussian",
  pedigree=isTree,
  prior = prior1,
  nitt = 1010, thin = 10, burnin = 10,
  verbose = FALSE
)

summary(inside.abundance.morph.comm)


outside.abundance.morph.comm <- MCMCglmm(
  fixed = log(abundance.surrounding+0.1) ~ beak.shape+locom.shape+poly(body.size,2)+hand_wing,
  random = ~ animal + species + city,
  data = UTI_within_species,
  family = "gaussian",
  pedigree=isTree,
  prior = prior1,
  nitt = 1010, thin = 10, burnin = 10,
  verbose = FALSE
)

summary(outside.abundance.morph.comm)



# This second model asks whether urban tolerant species tend to be morphologically different than urban avoiders

prior2 <- list(R=list(V=1, nu=0),
               G=list(G1=list(V=1, nu=0, alpha.mu=0, alpha.V=1e3)))

tolerance.diff.morph <- MCMCglmm(
  fixed = meanUTI ~ log(abundance.surrounding+0.1)+beak.shape+locom.shape+body.size+hand_wing,
  random = ~ animal,
  data = UTI_species,
  family = "gaussian",
  pedigree=isTree,
  prior = prior2,
  nitt = 1010, thin = 10, burnin = 10,
  verbose = FALSE
)

summary(tolerance.diff.morph)


inside.abundance.morph <- MCMCglmm(
  fixed = log(abundance.surrounding+0.1) ~ beak.shape+locom.shape+body.size+hand_wing,
  random = ~ animal,
  data = UTI_species,
  family = "gaussian",
  pedigree=isTree,
  prior = prior2,
  nitt = 1010, thin = 10, burnin = 10,
  verbose = FALSE
)

summary(inside.abundance.morph)


outside.abundance.morph <- MCMCglmm(
  fixed = log(abundance.surrounding+0.1) ~ beak.shape+locom.shape+body.size+hand_wing,
  random = ~ animal,
  data = UTI_species,
  family = "gaussian",
  pedigree=isTree,
  prior = prior2,
  nitt = 1010, thin = 10, burnin = 10,
  verbose = FALSE
)

summary(outside.abundance.morph)








### Area in preparation ##############

## Random forests to assess the relaibility of morphological traits as functional traits

library(randomForest)

fit1 <- randomForest(effort.corrected.abundance ~ propagule+beak.shape+body.size+hand_wing+locom.shape+community, data = dat.toler.urb, do.trace=10, ntree=100, labelVar=TRUE, importance=TRUE)

print(fit1)
plot(fit1)
round(importance(fit1), 2)
varImpPlot(fit1)

fit1$test


library(rpart)

fit <- rpart(effort.corrected.abundance ~ propagule+beak.shape+body.size+hand_wing+locom.shape, data = dat.toler.urb, method="poisson")
print(fit)








###############################
# Estimation of urban tolerance
###############################


### This script estimates the difference in abundance between habitats ###

# We first use relative abundance

dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)
dat <- subset(dat,used.urban.nonurban=="yes" & relative.abundance>0)
dat[] <- lapply(dat, function(x) if(is.factor(x)) factor(x) else x)

d <- dat[,c(2,3,5,7,8,17)]
levels(d$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")

d <- subset(d, city=="Amravati" | city=="Barcelona"| city=="Cameron_Highlands"| city=="Cayenne"| city=="Ciudad_de_Mexico")  
d[] <- lapply(d, function(x) if(is.factor(x)) factor(x) else x)




abundanceDiff<- by(d, d[, c("city", "animal")], function(x){
  diffs<- expand.grid(x$community, x$community)
  diffs<- diffs[diffs[,1] != diffs[,2],]
  res<- data.frame()
  for (i in 1:nrow(diffs)){
    tmp<- x[x$community %in% diffs[i,],]
    resTmp<- as.numeric(as.character(tmp$relative.abundance[2])) - as.numeric(as.character(tmp$relative.abundance[1]))
    resTmp<- data.frame(city=x$city[1], animal=x$animal[1], community1=diffs[i,1], community2=diffs[i,2],
                        habitat1=tmp$habitat[1], Rel.Abund.habitat1=as.numeric(as.character(tmp$relative.abundance[1])), habitat2=tmp$habitat[2], tolerance_hab2_relative_hab1=resTmp, stringsAsFactors = FALSE)
    res<- rbind(res, resTmp)
  }
  return(res)
})

abundanceDiff<- do.call(rbind, abundanceDiff)  

Urb.tol <- subset(abundanceDiff, habitat2=="urb")


# We eliminate duplicates

Urb.tol <- ddply(Urb.tol, c("location", "species", "habitat1"), summarise,
             urban.tolerance  = mean(tolerance_hab2_relative_hab1))

colnames(Urb.tol) <- c("location", "animal", "habitat.outside", "urban_tolerance") 

## Merged urban and morphological data:

dat.toler <-merge(Urb.tol, func, by="animal")
dat.toler$species <- dat.toler$animal
str(dat.toler)

## Random forests to assess the relaibility of morphological traits as functional traits

library(randomForest)

fit1 <- randomForest(urban_tolerance ~ propagule+beak.shape+body.size+hand_wing+locom.shape+community, data = dat.toler.urb, do.trace=10, ntree=100, labelVar=TRUE, importance=TRUE)

print(fit1)
plot(fit1)
round(importance(fit1), 2)
varImpPlot(fit1)

fit1$test


library(rpart)

fit <- rpart(effort.corrected.abundance ~ propagule+beak.shape+body.size+hand_wing+locom.shape, data = dat.toler.urb, method="poisson")
print(fit)

#################################################
#####  Modelling tolerance to urbanisation ######
#################################################

## Urban data:

# we will first use effort.corrected.abundance, which is the appropriate way to do so

dat<-read.table(paste0(workingData,"/Urban global data April 6 2018 for R.txt"), header=TRUE)
# dat0 <- subset(dat,status=="native")  # if we want to exclude exotics
# dat0[] <- lapply(dat0, function(x) if(is.factor(x)) factor(x) else x)

dat01 <- subset(dat, urban.analyses=="yes") # if we focus on communities within urbanisation gradients
dat02[] <- lapply(dat01, function(x) if(is.factor(x)) factor(x) else x)

levels(dat02$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")

dat02$link <- as.factor(paste(dat02$city,dat02$animal, dat02$habitat,sep = ""))  # we create a new variable to then link with subsequent data



## We estimate mean corrected relative abundance in the wildland

tmp <- ddply(dat02, c("city", "habitat", "animal"), summarise,
             propagule  = sum(effort.corrected.abundance)/sum(community.size))
tmp$link <- as.factor(paste(tmp$city,tmp$animal, tmp$habitat,sep = ""))  # we create a new variable to then link with subsequent data

tmp <- tmp[,-c(1:3)]

dat03 <- merge(dat02[,-19],tmp, by="link")



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

library(MCMCglmm)

zi.prior <-  list(R = list(V = diag(2), n = 0.002, fix = 2),
                  G = list(G1 = list(V = 1, n = 0.002),
                           G2 = list(V = 1, n = 0.002),
                           G3 = list(V = 1, n = 0.002)))

# residual variance for the zero-inflation is identified with rcov=idh(trait):units and fixing the [2,2] element of the covariance matrix at something (e.g  1) in the prior (i.e. add fix=2 to your R prior).
# Check whether using trait-1 or simply -1
# we may use zapoisson with a zero-altered model


## Test of traits influencing presence/abundance in cities

dat.toler.01 <- subset(dat.toler, habitat="urban")
dat.toler.urb <- na.omit(dat.toler.01[,c(1,5,12,14,18:31)])


m1 <- MCMCglmm(as.integer(effort.corrected.abundance) ~ trait-1 + at.level(trait,1):log(propagule+0.5) + at.level(trait,1):log(propagule+0.5) + at.level(trait,1):body.size + at.level(trait,1):beak.shape + at.level(trait,1):locom.shape + at.level(trait,1):hand_wing,    # trait-1 rather than -1 as we do not want the intercept for the Poisson process and the zero-inflation to be the same
               random = ~idh(at.level(trait,1)):community +   # at.level separates predictors in the Poisson and Binomial model
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






###############################
# Estimation of urban tolerance
###############################


### This script estimates the difference in abundance between habitats ###

dat<-read.table(paste0(workingData,"/Urban global data April 6 2018 for R.txt"), header=TRUE)
d <- dat[,c(3,5,7,9,17)]
d <- subset(d, location=="Amravati" | location=="NTL")  
d[] <- lapply(d, function(x) if(is.factor(x)) factor(x) else x)



abundanceDiff<- by(d, d[, c("location", "original_sp")], function(x){
  diffs<- expand.grid(x$community, x$community)
  diffs<- diffs[diffs[,1] != diffs[,2],]
  res<- data.frame()
  for (i in 1:nrow(diffs)){
    tmp<- x[x$community %in% diffs[i,],]
    resTmp<- as.numeric(as.character(tmp$relative.abundance[2])) - as.numeric(as.character(tmp$relative.abundance[1]))
    resTmp<- data.frame(location=x$location[1], species=x$original_sp[1], community1=diffs[i,1], community2=diffs[i,2],
                        habitat1=tmp$habitat[1], habitat2=tmp$habitat[2], tolerance_hab2_relative_hab1=resTmp, stringsAsFactors = FALSE)
    res<- rbind(res, resTmp)
  }
  return(res)
})

abundanceDiff<- do.call(rbind, abundanceDiff)  



## Estimation of urban tolerance


Urb.tol <- subset(abundanceDiff, habitat2=="urb")

# We eliminate duplicates

Urb.tol <- ddply(Urb.tol, c("location", "species", "habitat1"), summarise,
             urban.tolerance  = mean(tolerance_hab2_relative_hab1))

colnames(Urb.tol) <- c("location", "species", "habitat.outside", "urban_tolerance") 





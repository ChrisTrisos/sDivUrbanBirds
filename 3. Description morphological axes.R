###################################################
### 3.	Description of morphological niche axes ###
###################################################

# 1.1. PCA ALL 11 TRAITS
# 1.2. PCA BEAK TRAITS
# 1.3. PCA LOCOMOTOR TRAITS
# 1.4. PCA BODY SIZE
# 1.5. OUTPUT


################
### Libraries
################

library(ade4)
source("/Users/d.sol/Documents/Science/Research/Urbanisation/Functional diversity and urbanization/Data and analyses/Hypervolumes/utils.R")


###########
### Data                              
##########

setwd("~/Documents/Science/Research/Urbanisation/Functional diversity and urbanization/Data and analyses")

morph0 <- read.table("Morphological traits urban birds 24 Feb 2018 for R.txt", h=TRUE)
# morph0<- subset(morph0, animal!="Struthio_camelus")

morph <- morph0[,c(8,10,12:18)]
names(morph) <- c("mass", "bill_length", "bill_width", "bill_depth", "tarsus", "second","wing","hand_wing","tail")

# Log10-transformed and standardized to zero mean and unit variance (z-transformation).

morph$mass <- log(morph$mass)
morph$bill_length <- log(morph$bill_length)
morph$bill_width <- log(morph$bill_width)
morph$bill_depth <- log(morph$bill_depth)
morph$tarsus <- log(morph$tarsus)
morph$second <- log(morph$second)
morph$wing <- log(morph$wing)
morph$hand_wing <- log(morph$hand_wing)
morph$tail <- log(morph$tail)

morph <- as.data.frame(scale(morph))   # we scale the data


########################################################################
###                1. PCA TO ALL 11 TRAITS                           ###                 
########################################################################


###################################################
###                 Apply PCA
###################################################

pca11 <- dudi.pca(morph[-8], scannf = F, nf = 3)
scatter(pca11, clab.row=0, clab.col = 1.5, xax = 1, yax = 2, posieig="bottomright")

pca.11b <- princomp(morph[-8], cor=FALSE)  # as variables are standardised, we can use either the cor or var-covar mattrix
summary(pca.11b)
loadings(pca.11b)
body.size <- -pca.11b$scores[,1]
plot(body.size,morph$mass)


###################################################
###      Test the significance of PCA axes
###################################################
morph.rnd <- cbind(morph, rnd = rnorm(nrow(morph)))  
pca11.rnd <- dudi.pca(morph.rnd, scannf = F, nf = 3)
test11.rnd <- testdim(pca11.rnd, nrepet = 99)  # 999 for final analyses
test11.rnd


###################################################
###    Test the significance of PCA loadings
###################################################

# boot11 <- netoboot(morph.rnd, scannf = F, nf = 3)
# plotboot(pca11.rnd$c1, boot11)


###################################################
###    Analysis of the 11D volume
###################################################

# pc <- scale(as.matrix(morph))


###################################################
###  Compute observed volume (all and 95 % of species)
###################################################
# pc95 <- subselect.data(pc, 0.95)
#blob95 <- convhulln(pc95,"FA")
#obs.vol95 <- blob95$vol






########################################################################
###                2. PCA TO BEAK TRAITS                           ###                 
########################################################################

beak <- morph[,c(2:4)]
pca.beak <- dudi.pca(beak, scannf = F, nf = 3)
scatter(pca.beak, clab.row=0, clab.col = 1.5, xax = 1, yax = 2, posieig="bottomright")

pca.beak2 <- princomp(beak, cor=TRUE)
summary(pca.beak2)
loadings(pca.beak2)
beak.size <- pca.beak2$scores[,1]
plot(beak.size,body.size)

beak.shape <- pca.beak2$scores[,2]
plot(beak.shape,body.size)

############################################################################
###                3. PCA TO LOCOMOTOR TRAITS                           ###                 
############################################################################

locom <- morph[,c(5:7,9)]
pca.locom <- dudi.pca(locom, scannf = F, nf = 3)
scatter(pca.locom, clab.row=0, clab.col = 1.5, xax = 1, yax = 2, posieig="bottomright")

pca.locom2 <- princomp(locom, cor=TRUE)
summary(pca.locom2)
loadings(pca.locom2)

locom.size <- pca.locom2$scores[,1]
plot(locom.size,morph$mass)

locom.shape <- pca.locom2$scores[,2]
plot(locom.shape,morph$mass)



#################################################################
###             4. PCA BODY SIZE TWO STEPS                    ###                 
#################################################################

# does not work so well 

size <- as.matrix(beak.size, locom.size, morph$mass)

pca.size2 <- princomp(size, cor=TRUE)
summary(pca.size2)
loadings(pca.size2)
body.size2 <- pca.size2$scores[,1]

plot(pca.size2$scores[,1],morph$mass)
plot(pca.size2$scores[,1],locom.size)
plot(pca.size2$scores[,1],beak.size)



#################################################################
###                5. OUTPUT                                 ###                 
#################################################################


morphological.axes <- data.frame(morph0$OrigNam, morph0$animal, beak.size, beak.shape, locom.size, locom.shape, body.size, morph$mass, morph$bill_length, morph$bill_width, morph$bill_depth, morph$tarsus, morph$second, morph$wing, morph$hand_wing, morph$tail)

colnames(morphological.axes) <- c("animal", "species", "beak.size", "beak.shape", "locom.size", "locom.shape", "body.size", "mass", "bill_length", "bill_width", "bill_depth", "tarsus", "second", "wing", "hand_wing", "tail" )

write.table(morphological.axes,"morphological.axes.txt")

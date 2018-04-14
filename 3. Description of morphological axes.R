###################################################################
#################### 1. MORPHOLOGICAL ANALYSES ####################
###################################################################

# GOALS: 1) Estimate ecologically meaningful traits based on 9 morpohological traits

# INPUT: Species-level information on 9 traits

# OUTPUT: A file for subsequent analyses ("morphological.axes.txt"), which contains the 9 traits (including the hand-wing index) plus PCA components describing body size, beak (shape and size) and locomotor system (shape and size)




###########
### Data                              
##########

morph0 <- read.table("/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis/Morphological traits urban birds 24 Feb 2018 for R.txt", h=TRUE)
# read.table(paste0(workingData,"Morphological traits urban birds 24 Feb 2018 for R.txt"))

morph <- morph0[,c(8,10,12:18)] # we exclude redundant traits
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
###      PCA TO ALL 9 TRAITS TO ESTIMATE BODY SIZE                ###                 
########################################################################


pca.9 <- dudi.pca(morph[-8], scannf = F, nf = 3)  # We exclude hand_wing index
scatter(pca.9, clab.row=0, clab.col = 1.5, xax = 1, yax = 2, posieig="bottomright")

pca.9b <- princomp(morph[-8], cor=FALSE)  # as variables are standardised, we can use either the cor or var-covar mattrix
summary(pca.9b)
loadings(pca.9b)
body.size <- -pca.9b$scores[,1]
plot(body.size,morph$mass)




########################################################################
###                2. PCA TO BEAK TRAITS                            ###                 
########################################################################

beak <- morph[,c(2:4)]
pca.beak <- dudi.pca(beak, scannf = F, nf = 3)
scatter(pca.beak, clab.row=0, clab.col = 1.5, xax = 1, yax = 2, posieig="bottomright")

pca.beak2 <- princomp(beak, cor=TRUE)
summary(pca.beak2)
loadings(pca.beak2)
beak.size <- -pca.beak2$scores[,1]
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

locom.size <- -pca.locom2$scores[,1]
plot(locom.size,morph$mass)

locom.shape <- pca.locom2$scores[,2]  # Tarsus vs tail
plot(locom.shape,morph$mass)



#################################################################
###                5. OUTPUT                                 ###                 
#################################################################


morphological.axes <- data.frame(morph0$OrigNam, morph0$animal, beak.size, beak.shape, locom.size, locom.shape, body.size, morph$mass, morph$bill_length, morph$bill_width, morph$bill_depth, morph$tarsus, morph$second, morph$wing, morph$hand_wing, morph$tail)

colnames(morphological.axes) <- c("animal", "species", "beak.size", "beak.shape", "locom.size", "locom.shape", "body.size", "mass", "bill_length", "bill_width", "bill_depth", "tarsus", "second", "wing", "hand_wing", "tail" )

write.table(morphological.axes,"/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis/morphological.axes.txt")

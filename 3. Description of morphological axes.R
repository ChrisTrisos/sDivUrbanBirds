###################################################################
#################### MORPHOLOGICAL ANALYSES ####################
###################################################################

# GOALS: 1) Estimate ecologically meaningful traits based on 9 morpohological traits

# INPUT: Species-level information on 9 traits

# OUTPUT: 
    # A file for subsequent analyses ("morphological.axes.txt"), which contains the 9 traits (including the hand-wing index) plus PCA components describing body size, beak (shape and size) and locomotor system (shape and size)
    # An analysis of the relationship between morphology and diet

# SUMMARY:
# PCA TO ALL 9 TRAITS TO ESTIMATE BODY SIZE   
# PCA TO DESCRIBE BEAK SHAPE                                   
# PCA TO DESCRIBE HINDLIMB MORPHOLOGY 
# ANALYSIS OF DIET
# ANALYSIS OF DIET VS MORPHOLOGY
# CONVEX HULLS TO EXAMINE NICHE OVERLAP             
# CONVEX HULLS TO EXAMINE NICHE OVERLAP, FOR ABUNDANT SPECIES (>10%)


###########
### Data                              
##########

morph0<- read.table(paste0(workingData,"/Morphological traits urban birds 24 Feb 2018 for R.txt"), h=TRUE)

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
###      PCA TO ALL 9 TRAITS TO ESTIMATE BODY SIZE                   ###                 
########################################################################


pca.9 <- dudi.pca(morph[-8], scannf = F, nf = 3)  # We exclude hand_wing index
scatter(pca.9, clab.row=0, clab.col = 1.5, xax = 1, yax = 2, posieig="bottomright")

pca.9b <- princomp(morph[-8], cor=FALSE)  # as variables are standardised, we can use either the cor or var-covar mattrix
summary(pca.9b)
loadings(pca.9b)
body.size <- -pca.9b$scores[,1]
plot(body.size,morph$mass)




########################################################################
###           PCA TO DESCRIBE BEAK SHAPE                             ###                 
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
###           PCA TO DESCRIBE HINDLIMB MORPHOLOGY                        ###                 
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
###                    OUTPUT                                 ###                 
#################################################################


morphological.axes <- data.frame(morph0$OrigNam, morph0$animal, beak.size, beak.shape, locom.size, locom.shape, body.size, morph$mass, morph$bill_length, morph$bill_width, morph$bill_depth, morph$tarsus, morph$second, morph$wing, morph$hand_wing, morph$tail)

colnames(morphological.axes) <- c("animal", "species", "beak.size", "beak.shape", "locom.size", "locom.shape", "body.size", "mass", "bill_length", "bill_width", "bill_depth", "tarsus", "second", "wing", "hand_wing", "tail" )

write.table(morphological.axes,paste0(workingData,"/morphological.axes.txt"))




########################################################################
###               ANALYSIS OF DIET VS MORPHOLOGY                     ###                 
########################################################################

## We first merge diet and morphology

diet <-read.table(paste0(workingData,"/Diet urban birds 28 April 2018 for R.txt"), header=TRUE)
names(diet)
diet<-diet[,c(5,6:16)]

func<-read.table(paste0(workingData,"/morphological.axes.txt"),header=TRUE)
names(func)
all <-func[,c(1,3:16)]  # dataset for all traits when analysed separately

diet.morph <- merge(diet,all, by="animal")

diet.morph <- diet.morph[order(diet.morph$animal),] # wee need to order the species
names(diet.morph)
rownames(diet.morph) <- diet.morph$animal
diet.morph <- diet.morph[,-c(1,2)]

diet.morph <- na.omit(diet.morph[labels(comm[1,]),])  # we take only species in communities


## we create a file adding community information to know which diet types are sufficiently represented in each study area

dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)

diet.morph.tmp <- diet.morph
diet.morph.tmp$animal <- rownames(diet.morph)

morph.diet.comm <- merge(dat, diet.morph.tmp, by="animal")  # we add morphology to community abundance data

## we obtain mean values for species for descriptive purposes

tmp2 <- ddply(tmp, c("country", "city", "used.urban.nonurban", "animal"), summarise,
             Diet.Inv.Spp = mean(Diet.Inv),
             Diet.Vend.Spp = mean(Diet.Vend),
             Diet.Vect.Spp = mean(Diet.Vect),
             Diet.Vfish.Spp = mean(Diet.Vfish),
             Diet.Vunk.Spp = mean(Diet.Vunk),
             Diet.Scav.Spp = mean(Diet.Scav),
             Diet.Fruit.Spp = mean(Diet.Fruit),
             Diet.Nect.Spp = mean(Diet.Nect),
             Diet.Seed.Spp = mean(Diet.Seed),
             Diet.PlantO.Spp = mean(Diet.PlantO))

table(tmp2$city,tmp2$Diet.Inv)
table(tmp2$city,tmp2$Diet.Vend) 
table(tmp2$city,tmp2$Diet.Vect)  
table(tmp2$city,tmp2$Diet.Vfish) 
table(tmp2$city,tmp2$Diet.Vunk)
table(tmp2$city,tmp2$Diet.Scav) 
table(tmp2$city,tmp2$Diet.Fruit)
table(tmp2$city,tmp2$Diet.Nect)
table(tmp2$city,tmp2$Diet.Seed)
table(tmp2$city,tmp2$Diet.PlantO)




########################################################################
###            CONVEX HULLS TO EXAMINE NICHE OVERLAP                ###                 
########################################################################

## we draw some convex hulls for native species in particular cities, based on PCAs

# Newcastle

Newcastle <- subset(morph.diet.comm, morph.diet.comm$city=="Newcastle" & morph.diet.comm$relative.abundance>0 & morph.diet.comm$status=="native")
dataf <- Newcastle[,c(41,40,8)]
colnames(dataf) <- c("x","y","cat")

library(geometry)

chulls <- ddply(dataf, .(cat), function(df) df[chull(df$x, df$y), ])

# plot polygons

ggplot(data=dataf, aes(x=x, y=y, color=cat)) + geom_point() +
  geom_polygon(data=chulls, aes(x=x, y=y, group=cat), fill=NA)

# if you wanna get fancy with the colors
# (will give annoying warnings if number of groups is smaller than 3, ignore)

ggplot(data=dataf, aes(x=x, y=y, color=cat)) + geom_point() +
  geom_polygon(data=chulls, aes(x=x, y=y, fill=cat), alpha=0.2) +
  scale_color_brewer(palette='Set1') +
  scale_fill_brewer(palette='Set1')




#######################################################################################
###            CONVEX HULLS TO EXAMINE NICHE OVERLAP FOR ABUNDANT SPECIES           ###                 
#######################################################################################

## we draw some convex hulls for native species in particular cities, based on PCAs

# Newcastle

Newcastle <- subset(morph.diet.comm, morph.diet.comm$city=="Newcastle" & morph.diet.comm$relative.abundance>0.001 & morph.diet.comm$status=="native")
dataf <- Newcastle[,c(41,40,8)]
colnames(dataf) <- c("x","y","cat")

library(geometry)

chulls <- ddply(dataf, .(cat), function(df) df[chull(df$x, df$y), ])

# plot polygons

ggplot(data=dataf, aes(x=x, y=y, color=cat)) + geom_point() +
  geom_polygon(data=chulls, aes(x=x, y=y, group=cat), fill=NA)

# if you wanna get fancy with the colors
# (will give annoying warnings if number of groups is smaller than 3, ignore)

ggplot(data=dataf, aes(x=x, y=y, color=cat)) + geom_point() +
  geom_polygon(data=chulls, aes(x=x, y=y, fill=cat), alpha=0.2) +
  scale_color_brewer(palette='Set1') +
  scale_fill_brewer(palette='Set1')




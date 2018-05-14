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

diet.morph <- merge(all, diet, by="animal")

diet.morph <- diet.morph[order(diet.morph$animal),] # wee need to order the species
names(diet.morph)
rownames(diet.morph) <- diet.morph$animal
diet.morph <- diet.morph[,-c(1,2)]

#diet.morph <- na.omit(diet.morph[labels(comm[1,]),])  # we take only species in communities


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



#############################################################################################
###       Principal Coordinate analysis to decribe morphological and resource niches      ###                 
#############################################################################################

library(ks)

#trait data
traits0 <- diet.morph # this data.frame has been built in the prtevious section

row.names(traits0) <- traits0$animal

# traits are coded as:
# 1. Quantitative "Q": beak.size"   "beak.shape"  "locom.size"  "locom.shape" "body.size"   "mass"        "bill_length" "bill_width"  "bill_depth"  "tarsus"      "second"      "wing"        "hand_wing"   "tail" 
# 2. Nominal "N": None
# 3. Fuzzy codes "F": "Diet.Inv"    "Diet.Vend"   "Diet.Vect"   "Diet.Vfish"  "Diet.Vunk"   "Diet.Scav"   "Diet.Fruit"  "Diet.Nect"   "Diet.Seed"  "Diet.PlantO"
# 4. Dichotomus "D": None

# We will compute the functional distances between species using the modified Gower' distance as in Pavoine et al., 2009
# all analyses are based on the function "dist.ktab"

# source("dist.ktab.R")

# preparation of the variables
# 1. Quantitative "Q":
tabQ <- traits0[,11:24]

# 2. Fuzzy codes "F": 
tabF <- traits0[, c(1:10)]
tabFp <- prep.fuzzy(tabF, 10)
# tabFp <- prep.fuzzy(tabF, c(10, 33))  # if we add foraging strategy with 33 levels


# Calculation of the global functional distances between species

ktab1 <- ktab.list.df(list(tabQ, tabFp))
distrait <- dist.ktab(ktab1, c("Q", "F"))

# ktab1 <- ktab.list.df(list(tabQ, tabN, tabFp))
# distrait <- dist.ktab(ktab1, c("Q", "N", "F", "D"))
summary(distrait)

#euclidean properties
is.euclid(distrait)

#transformation for Euclidean distances
disl <- lingoes(distrait)
summary(disl)

#calculate the number of dimensions for further analyses
# source("quality_funct_space_fromdist.R")

qual_fs_r<-quality_funct_space_fromdist(disl,  nbdim=15) 
qual_fs_r$meanSD # 7D 0.006718
K<-4 # despite the analyses showed significance at the 7D axes, we finally kept 4 dimensions after assessing the eigenvalues in the PCoA analysis


#Principal Coordiantes Analyses
pcodisl <- dudi.pco(disl, scan = F, nf=4)
barplot(pcodisl$eig) # Eigenvalue barplot
pcodisl$eig
pcodisl$eig[1:4]/sum(pcodisl$eig) # variance explained by each axes
s.label(pcodisl$li, clabel = 0.7)
summary (pcodisl)

# correlation between original trais and pco axes, saving only 2 decimals
cor.table1<-round(cor(traits0[,],pcodisl$li,method="spearman",use="pairwise.complete.obs"),2)
cor.table1
# write.table(cor.table1,"pco_cor_table_r1.txt",sep="\t") # Saving table to interpret PCoA axes

# pdf("Scatter plot.pdf")
scatter(pcodisl)
scatter(pcodisl, clab.row=0, clab.col = 1.5, xax = 1, yax = 2, posieig="topleft")
dev.off()

# kernel density estimation on the first two main axes based on Diaz et al 2016
pc12<-pcodisl$li[,1:2]
H <- Hpi(x=pc12)      # optimal bandwidth estimation
est<- kde(x=pc12, H=H, compute.cont=TRUE)     

# set contour probabilities for drawing contour levels
cl<-contourLevels(est, prob=c(0.5, 0.05, 0.001), approx=TRUE)

fit<-envfit(pc12, traits0) # use envfit for drawing arrows, can be also done using trait loadings
fit2<-fit$vectors$arrows*-1 # drawing line segments in arrow opposites direction for pretty layout

# pdf("_PCoA_Iberian.pdf", width=8, height=6)
par(mar=c(4,4,2,2))
plot(est, cont=seq(1,100,by=1), display="filled.contour2", add=FALSE, ylab="PC2", xlab="PC1", cex.axis=0.75, ylim=c(-0.5, 0.5), xlim=c(-0.5, 0.5),las=1) 
plot(est,abs.cont=cl[1], labels=c(0.5),labcex=0.75, add=TRUE, lwd=0.75, col="grey30")
plot(est,abs.cont=cl[2], labels=c(0.95),labcex=0.75, add=TRUE, lwd=0.5, col="grey60")
plot(est,abs.cont=cl[3], labels=c(0.99),labcex=0.75, add=TRUE, lwd=0.5, col="grey60")
points( pc12[,], pch=16, cex=0.25, col="black") 
dev.off()


#The variables are displayed on the factorial map by Family
s.class(pcodisl$li, traits0$BLFamilyLatin, sub = "Family")







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

# convex hull

ggplot(data=dataf, aes(x=x, y=y, color=cat)) + geom_point() +
  geom_polygon(data=chulls, aes(x=x, y=y, fill=cat), alpha=0.2) +
  scale_color_brewer(palette='Set1') +
  scale_fill_brewer(palette='Set1')


# Kernell

m <- ggplot(data=dataf, aes(x=x, y=y, colour=factor(cat))) + 
  geom_point(alpha = 0.2) + 
    xlim(-7, 10) +
    ylim(-2.5, 2.5) 
m + stat_density_2d(aes(color=cat), geom = "polygon")




g <- ggplot(data=dataf, aes(x=x, y=y, colour=factor(cat))) + 
  geom_point(alpha = 0.2) +
  stat_density_2d(aes(fill = ..cat..), alpha=.1, geom = "polygon")

dat_lims <- lapply(dataf[,c(1,2)], function(v) c(min(v), max(v)))
plot_lims <- ggplot_build(g)$panel$ranges[[1]][c("x.range", "y.range")]

g +
  scale_x_continuous(limits = dataf$x * 1.1) + 
  scale_y_continuous(limits = dataf$y * 1.1) +
  coord_cartesian(xlim(-3.439797, 9.103873), ylim(-1.654770, 1.440386))



ggplot(data=dataf, aes(x=x, y=y, colour=factor(cat))) + 
  geom_point(alpha = 0.2) + 
  xlim(-7, 10) +
  ylim(-2.5, 2.5) +
  geom_density_2d(data=dataf, aes(x=x, y=y, colour=factor(cat))) + theme(panel.background = element_rect(fill = '#ffffff'))




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




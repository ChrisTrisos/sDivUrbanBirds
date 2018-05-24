############################################################################
#################### 4.	Description of diet niche axes  ####################
############################################################################

# GOALS: 1) Estimate ecologically meaningful traits based on diet categories by means of PCoA
       # 2) Estimate degree of foraging specialization

# SUMMARY:


#############################################################################################
###       Principal Coordinate analysis to decribe morphological and resource niches      ###                 
#############################################################################################

library(ks)

#trait data

diet <-read.table(paste0(workingData,"/Diet urban birds 28 April 2018 for R.txt"), header=TRUE)
# Missing data Thalasseus_eurygnatha, Saxicola_torquatus, already removed from diet
names(diet)
diet <- diet[order(diet$animal),] # wee need to order the species
diet<-diet[,c(5,6:16)]
names(diet)
rownames(diet) <- diet$animal
diet <- diet[,-c(1:2)]


## Estimation of distances

distdiet<- dist.prop(diet, method = 1) # method 1 is Manly
distdiet <- distdiet/max(distdiet) 


#euclidean properties
is.euclid(distdiet)

#transformation for Euclidean distances
disl <- lingoes(distdiet)
summary(disl)


#Principal Coordiantes Analyses
pcodisl <- dudi.pco(disl, scan = F, nf=4)
barplot(pcodisl$eig) # Eigenvalue barplot
pcodisl$eig
pcodisl$eig[1:4]/sum(pcodisl$eig) # variance explained by each axes
s.label(pcodisl$li, clabel = 0.7)
summary (pcodisl)

# correlation between original trais and pco axes, saving only 2 decimals
cor.table1<-round(cor(diet[,],pcodisl$li,method="spearman",use="pairwise.complete.obs"),2)
cor.table1 

colnames(pcodisl$li) <- c("Inv_SeedFruit", "Inv_Seed", "Seed", "not_used")
write.table(pcodisl$li,paste0(workingData,"/diet.axes.txt"))



## Figures

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





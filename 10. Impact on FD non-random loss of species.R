###########################################################################################################################
#######  10.	Impact on FD of non-random loss of species associated with functional redundancies, abundance and both  #####
###########################################################################################################################


### GOALS: ###

# Model how fuctional diversity metrics co-vary with species richness 

### INPUTS: ### 

# Morphologically-based diversity metrics for communities ("Morphological diversity metrics for communities.txt")
# Diet-based diversity metrics for communities ("Morphological diversity metrics for communities.txt")
# Foraging behaviour data-based diversity metrics for communities ("Morphological diversity metrics for communities.txt")


### OUTPUTS: ###



### ANALYSES: ###

## 1. Analyses for morphological diversity
## 2. Analyses for morphological diversity, only natives




########### 1. Analyses for morphological diversity ##############
##################################################################

{## Import biodiversity metrics for communities
  
     x<-read.table(paste0(workingData,"/Morphological diversity metrics for communities.txt"))
     # x <- read.table("Morphological diversity metrics for communities natives.txt")  # metrics estimated excluding exotics
 
  
  ## To compare both FD and S values across studies, diversities of each land use were z-scored within each study site (i.e. subtracting study mean and dividing by the study standard deviation).
  
  cdata <- ddply(dat0, c("country", "city", "community", "habitat"), summarise,
                 N    = length(relative.abundance))
  
  tmp <- ddply(x,"city", summarise,
        QE.taxonomy.mean =mean(QE.taxonomy),
        QE.taxonomy.sd=sd(QE.taxonomy))
  
x1 <- merge(x, tmp, by="city")
  
  # Three ways to code habitats
  
  # 1. All habitats separated
  # levels(x$habitat) <- c("closed_wild", "Urban_Park", "little_urbanised", "open_wild", "pasture", "plantation", "rural", "rural_wild", "sub", "urb", "urban_mosaic", "wild_mosaic")
  
  # 2. All urban habitats separated, all non-urban habitat together
  levels(x$habitat) <- c("Wildland",       "Urban_Park",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Suburban", "Urban", "Suburban", "Wildland")
  habitat.ordered  = factor(x$habitat, levels=c("Wildland","Rural","Urban_Park","Suburban","Urban"))
  x <- cbind(x,habitat.ordered)
  
  # 3. All urban and non-urban habitats pooled together
  # levels(habitat.city) <- c("Wildland",       "Urban",     "Wildland",     "Wildland", "Rural",   "Rural",  "Rural",   "Rural", "Urban", "Urban", "Urban", "Wildland")
  
  
  
  ## Tests of the effect of species richness on FD urbanization on biodiversity ##
  #################################################################################
  
  
  
  # 8 functional traits vs Species richness

  
  ggplot(x, aes(QE.all.morph, Species.richness)) +
    geom_point() +
    theme_minimal()
  
  
  ggplot(x, aes(QE.all.morph, log(Species.richness)-log(Regional.spp.richness))) +
    geom_point() +
    theme_minimal()
  
  QE.all.morph.spp = lme(QE.all.morph ~ habitat.ordered+Species.richness, random = ~ 1|country/city, data = x, method="ML")
  QE.all.morph.inter.spp = lme(QE.all.morph ~ habitat.ordered*Species.richness, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.all.morph.spp)
  summary(QE.all.morph.inter.spp)
  anova(QE.all.morph.inter.spp,QE.all.morph.spp)
  
  testInteractions(QE.all.morph.spp)
  plot(QE.all.morph.spp)
  qqnorm(QE.all.morph.spp)
  QE.all.morph.I <- interactionMeans(QE.all.morph.spp) # effect plots
  plot(QE.all.morph.I.spp, errorbar="ci95")
  
  
  
  # 8 functional traits vs Simpson's index
  
  ggplot(x, aes(QE.all.morph, QE.taxonomy)) +  #QE.taxonomy is the Simpson index
    geom_point() +
    theme_minimal()
  
  QE.all.morph.SI = lme(QE.all.morph ~ habitat.ordered+QE.taxonomy, random = ~ 1|country/city, data = x, method="ML")
  QE.all.morph.inter.SI = lme(QE.all.morph ~ habitat.ordered*QE.taxonomy, random = ~ 1|country/city, data = x, method="ML")
  summary(QE.all.morph.SI)
  summary(QE.all.morph.inter.SI)
  anova(QE.all.morph.inter.SI,QE.all.morph.SI)
  
  testInteractions(QE.all.morph.SI)
  plot(QE.all.morph.SI)
  qqnorm(QE.all.morph.SI)
  QE.all.morph.I.SI <- interactionMeans(QE.all.morph.SI) # effect plots
  plot(QE.all.morph.I.SI, errorbar="ci95")
  
  
  
  
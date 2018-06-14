########### Test of Beta diversity, only natives ##############
###############################################################

### Morphology Beta Q test

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Rao beta diversity morphology natives.txt"))
  
  y <- subset(x, hab.comp!="Urban_Urban_Park" & hab.comp!="Urban_Park_Urban_Park" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Urban_Park_Wildland" & hab.comp!="Urban_Park_Rural" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Wildland_Wildland" & hab.comp!="Suburban_Suburban" & hab.comp!="Urban_Urban" & hab.comp!="Rural_Rural")
  y[] <- lapply(y, function(x) if(is.factor(x)) factor(x) else x)
  
  habitat.ordered  = factor(y$hab.comp, levels=c("Rural_Wildland","Suburban_Rural","Suburban_Wildland","Urban_Wildland","Urban_Rural","Urban_Suburban"))
  y <- cbind(y,habitat.ordered)
  
  habitat.ordered  = factor(y$hab.comp, levels=c("Rural_Wildland","Suburban_Rural","Suburban_Wildland","Urban_Wildland","Urban_Rural","Urban_Suburban"))
  y <- cbind(y,habitat.ordered)
  
  
  
  ## Tests of the effect of urbanization on beta biodiversity*
  ################################################################
  
  
  # Changes in functional composition
  
  Qbetast.diff.morph = lme(Qbetast ~ habitat.ordered, random = ~ 1|city, data = y, method="ML")
  summary(Qbetast.diff.morph)
  testInteractions(Qbetast.diff.morph, pairwise="habitat.ordered")
  plot(Qbetast.diff.morph)
  qqnorm(Qbetast.diff.morph)
  Qbetast.diff.morph.I <- interactionMeans(Qbetast.diff.morph) # effect plots
  plot(Qbetast.diff.morph.I, errorbar="ci95")
  
  
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(Qbetast.diff.morph.I, aes(x= habitat.ordered, y=Qbetast.diff.morph.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Qbetast.diff.morph.I[,2]-Qbetast.diff.morph.I[,3], ymax=Qbetast.diff.morph.I[,2]+Qbetast.diff.morph.I[,3]), width=.2) +
    geom_point(data=Qbetast.diff.morph.I, mapping=aes(x=habitat.ordered, y=Qbetast.diff.morph.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Beta Q morphology", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b","a")))
  
 } 
  



### Foraging Beta Q test

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Rao beta diversity foraging natives.txt"))
  
  y <- subset(x, hab.comp!="Urban_Urban_Park" & hab.comp!="Urban_Park_Urban_Park" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Urban_Park_Wildland" & hab.comp!="Urban_Park_Rural" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Wildland_Wildland" & hab.comp!="Suburban_Suburban" & hab.comp!="Urban_Urban" & hab.comp!="Rural_Rural")
  y[] <- lapply(y, function(x) if(is.factor(x)) factor(x) else x)
  
  habitat.ordered  = factor(y$hab.comp, levels=c("Rural_Wildland","Suburban_Rural","Suburban_Wildland","Urban_Wildland","Urban_Rural","Urban_Suburban"))
  y <- cbind(y,habitat.ordered)
  
  ## Tests of the effect of urbanization on beta biodiversity*
  ################################################################
  
  
  # Changes in functional composition
  
  Qbetast.diff.forag = lme(Qbetast ~ habitat.ordered, random = ~ 1|city, data = y, method="ML")
  summary(Qbetast.diff.forag)
  testInteractions(Qbetast.diff.forag, pairwise="habitat.ordered")
  plot(Qbetast.diff.forag)
  qqnorm(Qbetast.diff.forag)
  Qbetast.diff.forag.I <- interactionMeans(Qbetast.diff.forag) # effect plots
  plot(Qbetast.diff.forag.I, errorbar="ci95")
  
  
  
  # Final figure (effect +/- standard error)
  
  b <- ggplot(Qbetast.diff.forag.I, aes(x= habitat.ordered, y=Qbetast.diff.forag.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Qbetast.diff.forag.I[,2]-Qbetast.diff.forag.I[,3], ymax=Qbetast.diff.forag.I[,2]+Qbetast.diff.forag.I[,3]), width=.2) +
    geom_point(data=Qbetast.diff.forag.I, mapping=aes(x=habitat.ordered, y=Qbetast.diff.forag.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Beta Q foraging", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b","a")))
  
}
  


### Diet Beta Q test

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Rao beta diversity diet natives.txt"))
  
  y <- subset(x, hab.comp!="Urban_Urban_Park" & hab.comp!="Urban_Park_Urban_Park" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Urban_Park_Wildland" & hab.comp!="Urban_Park_Rural" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Wildland_Wildland" & hab.comp!="Suburban_Suburban" & hab.comp!="Urban_Urban" & hab.comp!="Rural_Rural")
  y[] <- lapply(y, function(x) if(is.factor(x)) factor(x) else x)
  
  habitat.ordered  = factor(y$hab.comp, levels=c("Rural_Wildland","Suburban_Rural","Suburban_Wildland","Urban_Wildland","Urban_Rural","Urban_Suburban"))
  y <- cbind(y,habitat.ordered)
  
  ## Tests of the effect of urbanization on beta biodiversity*
  ################################################################
  
  
  # Changes in functional composition
  
  Qbetast.diff.diet = lme(Qbetast ~ habitat.ordered, random = ~ 1|city, data = y, method="ML")
  summary(Qbetast.diff.diet)
  testInteractions(Qbetast.diff.diet, pairwise="habitat.ordered")
  plot(Qbetast.diff.diet)
  qqnorm(Qbetast.diff.diet)
  Qbetast.diff.diet.I <- interactionMeans(Qbetast.diff.diet) # effect plots
  plot(Qbetast.diff.diet.I, errorbar="ci95")
  
  
  
  # Final figure (effect +/- standard error)
  
  c <- ggplot(Qbetast.diff.diet.I, aes(x= habitat.ordered, y=Qbetast.diff.diet.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Qbetast.diff.diet.I[,2]-Qbetast.diff.diet.I[,3], ymax=Qbetast.diff.diet.I[,2]+Qbetast.diff.diet.I[,3]), width=.2) +
    geom_point(data=Qbetast.diff.diet.I, mapping=aes(x=habitat.ordered, y=Qbetast.diff.diet.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Beta Q diet", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b","a")))
  
}


  
  
   tiff(paste0(GoogleFigs,"/plot_Beta_diversity_natives.tiff"), width = 9, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a,c,b, cols=1)
  dev.off()
  
  
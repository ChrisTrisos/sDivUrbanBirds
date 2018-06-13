########### Test of Beta diversity, only natives ##############
###############################################################

{## Import biodiversity metrics for communities
  
  x<-read.table(paste0(workingData,"/Rao beta diversity.txt"))
  
  # We restrict the analyses to studies where there is information inside and outside the city
  x <- subset(x, used.urban.nonurban=="yes")
  x[] <- lapply(x, function(x) if(is.factor(x)) factor(x) else x)
  
  y <- subset(x, hab.comp!="Urban_Urban_Park" & hab.comp!="Urban_Park_Urban_Park" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Urban_Park_Wildland" & hab.comp!="Urban_Park_Rural" & hab.comp!="Suburban_Urban_Park" & hab.comp!="Wildland_Wildland" & hab.comp!="Suburban_Suburban" & hab.comp!="Urban_Urban" & hab.comp!="Rural_Rural")
  y[] <- lapply(y, function(x) if(is.factor(x)) factor(x) else x)
  
  habitat.ordered  = factor(y$hab.comp, levels=c("Rural_Wildland","Suburban_Rural","Suburban_Wildland","Urban_Wildland","Urban_Rural","Urban_Suburban"))
  y <- cbind(y,habitat.ordered)
  
  ## Tests of the effect of urbanization on beta biodiversity*
  ################################################################
  
  
  # Changes in functional composition
  
  Qbetast.diff = lme(Qbetast ~ habitat.ordered, random = ~ 1|city, data = y, method="ML")
  summary(Qbetast.diff)
  testInteractions(Qbetast.diff, pairwise="habitat.ordered")
  plot(Qbetast.diff)
  qqnorm(Qbetast.diff)
  Qbetast.diff.I <- interactionMeans(Qbetast.diff) # effect plots
  plot(Qbetast.diff.I, errorbar="ci95")
  
  
  
  # Final figure (effect +/- standard error)
  
  a <- ggplot(Qbetast.diff.I, aes(x= habitat.ordered, y=Qbetast.diff.I[,2])) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    geom_errorbar(aes(ymin=Qbetast.diff.I[,2]-Qbetast.diff.I[,3], ymax=Qbetast.diff.I[,2]+Qbetast.diff.I[,3]), width=.2) +
    geom_point(data=Qbetast.diff.I, mapping=aes(x=habitat.ordered, y=Qbetast.diff.I[,2]), size=8, shape=21, fill="white") +
    labs(x = "", y = "Beta functional diversity", cex=16) +
    geom_text(aes(label= c("a","a","a","b","b","a")))
  
   tiff(paste0(GoogleFigs,"/plot_Beta_diversity_natives.tiff"), width = 9, height = 8, units = 'in', res = 200)
  ggplot2.multiplot(a, cols=1)
  dev.off()
  
  
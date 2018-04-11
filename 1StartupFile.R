################
## Set up file
################

## This file sets the wd for each user
##and loads functions for the rest of the script

#load libraries
library(picante)
library(adiv)


# set wd path for each user

# Chris
if(Sys.info()['user']=='ctrisos') {
  wd <-'/Users/ctrisos/Google Drive/sDivUrbBirds/'
  workingData <-'/Users/ctrisos/Google Drive/sDivUrbBirds/Data/DataForAnalysis'
  rawData<-'/Users/ctrisos/Google Drive/sDivUrbBirds/Data/OriginalData'
  #registerDoParallel(2) # choose how many cores it can use
  # for output to google drive
  GoogleOut<-'/Users/ctrisos/Google Drive/sDivUrbBirds/AnalysisOutputs'
  GoogleFigs<-'/Users/ctrisos/Google Drive/sDivUrbBirds/Figures'
}

# Alienor
if(Sys.info()['user']=='ctrisos') {
  wd <-'/Users/ctrisos/Google Drive/sDivUrbBirds/'
  workingData <-'/Users/ctrisos/Google Drive/sDivUrbBirds/Data/DataForAnalysis'
  rawData<-'/Users/ctrisos/Google Drive/sDivUrbBirds/Data/OriginalData'
  #registerDoParallel(2) # choose how many cores it can use
  # for output to google drive
  GoogleOut<-'/Users/ctrisos/Google Drive/sDivUrbBirds/AnalysisOutputs'
  GoogleFigs<-'/Users/ctrisos/Google Drive/sDivUrbBirds/Figures'
}

#Sandrine
if(Sys.info()['user']=='ctrisos') {
  wd <-'/Users/ctrisos/Google Drive/sDivUrbBirds/'
  workingData <-'/Users/ctrisos/Google Drive/sDivUrbBirds/Data/DataForAnalysis'
  rawData<-'/Users/ctrisos/Google Drive/sDivUrbBirds/Data/OriginalData'
  #registerDoParallel(2) # choose how many cores it can use
  # for output to google drive
  GoogleOut<-'/Users/ctrisos/Google Drive/sDivUrbBirds/AnalysisOutputs'
  GoogleFigs<-'/Users/ctrisos/Google Drive/sDivUrbBirds/Figures'
}

#Daniel
if(Sys.info()['user']=='ctrisos') {
  wd <-'/Users/ctrisos/Google Drive/sDivUrbBirds/'
  workingData <-'/Users/ctrisos/Google Drive/sDivUrbBirds/Data/DataForAnalysis'
  rawData<-'/Users/ctrisos/Google Drive/sDivUrbBirds/Data/OriginalData'
  #registerDoParallel(2) # choose how many cores it can use
  # for output to google drive
  GoogleOut<-'/Users/ctrisos/Google Drive/sDivUrbBirds/AnalysisOutputs'
  GoogleFigs<-'/Users/ctrisos/Google Drive/sDivUrbBirds/Figures'
}

## Any user independent file paths can be here

#### Any custom functions we all want can be loaded here

################
## Set up file
################
# The work flow is a series of R scripts. They are numbered by the order in which they are run. See the work flow document 
# in the google drive.

# Always start by running this set up file first. It will set all of the file paths for your computer.
# It will also load any custom functions that we all need for the analysis.


#load libraries
library(picante)
library(adiv)


# set wd path for each user
# wd is the file path to the google drive
# workingData is the file path to the data for analysis
# rawData is the file path to the original data (this is not updated unless an error is corrected)
# GoogleOut is the file path you use if you want to write analysis outputs to google drive
# GoogleFigs is the file path you use if you want to print a figure to the google drive


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
if(Sys.info()['user']=='d.sol') {
  wd <-'/Users/d.sol/Google Drive/sDivUrbBirds'
  workingData <-'/Users/d.sol/Google Drive/sDivUrbBirds/Data/DataForAnalysis'
  rawData<-'/Users/d.sol/Google Drive/sDivUrbBirds/Data/OriginalData'
  #registerDoParallel(2) # choose how many cores it can use
  # for output to google drive
  GoogleOut<-'/Users/d.sol/Google Drive/sDivUrbBirds/AnalysisOutputs'
  GoogleFigs<-'/Users/d.sol/Google Drive/sDivUrbBirds/Figures'
}

## Any user independent file paths can be here

#### Any custom functions we all want can be loaded here

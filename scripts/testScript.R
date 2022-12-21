cd = strsplit(getwd(),"/")[[1]]
if(cd[length(cd)]=="BayesianNationalBorealCaribouPVA"){
  setwd(paste0(getwd(),"/DemographyWBayesianUpdating/EakerModified/CaribouDemo_v1.15"))
}
wdir <- getwd()

dir.create("figs")
dir.create("tabs")

source("CaribouDemoFns.R")

# usage
packages <- c("shiny", "R2jags","gdata","mcmcplots",
  "ggplot2", "rmarkdown", "knitr",
  "RODBC","plyr","survival","gdata","data.table")

ipak(packages)
load.module("glm")

library(caribouMetrics)
library(tidyr)
library(dplyr)

###############
#Use Eacker example data for collaring parameters
eParsIn = list()
eParsIn$cowCounts <- data.frame(Year = 1981:2023,
                                Count = 100,
                                Class = "cow")
eParsIn$freqStartsByYear <- data.frame(Year = 1981:2023,
                                       numStarts = 30)
eParsIn$collarOnTime=1
eParsIn$collarOffTime=12
eParsIn$collarNumYears=3

adjustR = T #adjust recruitment for delayed age of first reproduction or no.

##########
#Get full set of sims for comparison
simBig<-getSimsNational(adjustR=adjustR)#If called with default parameters, use saved object to speed things up.

###############
#Step 1: confirm appropriate prior variability in survival intercept using minimal (2) observed data points & 0 fire/anthro covariates. Controlled by priors on l.Saf, phi and sig.Saf.
#################
#source("CaribouDemoFns.R")
#eParsIn$collarNumYears=1


str(eParsIn)
scns=expand.grid(P=4,st=30,cmult=3,ri=2,iA=0,sQ=0.24,rQ=0.3,N0=10000)
scResults = runScnSet(scns,eParsIn,simBig,getKSDists=F)

str(scResults)$obs.all
print(plotRes(scResults$rr.summary.all, "Population growth rate",obs=scResults$obs.all,
              lowBound=0,simRange = scResults$sim.all,facetVars=c("P","sQ")))

print(plotRes(scResults$rr.summary.all, "Recruitment",obs=scResults$obs.all,
              lowBound=0,simRange = scResults$sim.all,facetVars=c("P","sQ")))

print(plotRes(scResults$rr.summary.all, "Adult female survival",obs=scResults$obs.all,
              lowBound=0.65,simRange = scResults$sim.all,facetVars=c("P","sQ")))

print(plotRes(scResults$rr.summary.all, "Female population size",obs=scResults$obs.all,
              lowBound=0,highBound=2000,facetVars=c("P","sQ")))

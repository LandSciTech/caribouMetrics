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
eParsIn = getParamsFromEacker(path="C:/Users/HughesJo/Documents/gitprojects/BayesianNationalBorealCaribouPVA/DemographyWBayesianUpdating/EakerModified/CaribouDemo_v1.15")

eParsIn$cowCounts <- data.frame(Year = 1981:2018,
                        Count = 100,
                        Class = "cow")
eParsIn$freqStartsByYear <- data.frame(Year = 1981:2018,
                               numStarts = 25)


##########
#Get full set of sims for comparison
simBig<-getSimsNational()#If called with default parameters, use saved object to speed things up.

###############
#Step 1: confirm appropriate prior variability in survival intercept using minimal (2) observed data points & 0 fire/anthro covariates. Controlled by priors on l.Saf, phi and sig.Saf.
#################
#source("CaribouDemoFns.R")
#eParsIn$collarNumYears=1

str(eParsIn)
numObsYrs=c(10);startsByYr = 15 #25
scns=expand.grid(P=numObsYrs,sQ=c(0.5),st=startsByYr)
scResults = runScnSet(scns,eParsIn,simBig)
print(plotRes(scResults$rr.summary.all, "Female population size",obs=scResults$obs.all,
              lowBound=0,highBound=2000,facetVars=c("P","sQ")))

print(plotRes(scResults$rr.summary.all, "Adult female survival",obs=scResults$obs.all,
              lowBound=0.6,simRange=scResults$sim.all,facetVars=c("P","sQ")))
print(plotRes(scResults$rr.summary.all, "Population growth rate",obs=scResults$obs.all,
              lowBound=0,simRange=scResults$sim.all,facetVars=c("P","sQ")))



numObsYrs=c(2,5,20);startsByYr = 25;lse=1
scns=expand.grid(P=numObsYrs,sQ=c(0.025,0.5,0.975),st=startsByYr,lse=lse)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"lse",lse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","sQ"))

numObsYrs=c(2,5,20);startsByYr = 25;sse=0.1
scns=expand.grid(P=numObsYrs,sQ=c(0.025,0.5,0.975),st=startsByYr,sse=sse)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"sse",sse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","sQ"))

numObsYrs=c(2,5,20);startsByYr = 25;ssv=0.08696
scns=expand.grid(P=numObsYrs,sQ=c(0.025,0.5,0.975),st=startsByYr,pse=pse)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"ssv",ssv)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","sQ"))

numObsYrs=c(2,5,20);startsByYr = 25
scns=expand.grid(P=numObsYrs,sQ=c(0.025,0.5,0.975),st=startsByYr)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","sQ"))

numObsYrs=c(2,5,20);startsByYr = 5
scns=expand.grid(P=numObsYrs,sQ=c(0.025,0.5,0.975),st=startsByYr)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","sQ"))

numObsYrs=c(2,5,20);startsByYr = 25; iA=90
scns=expand.grid(P=numObsYrs,sQ=c(0.025,0.5,0.975),iA=iA,st=startsByYr)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"iA",iA)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","sQ"))

numObsYrs=c(2,5,20);startsByYr = 25; iA=90; bse=1
scns=expand.grid(P=numObsYrs,sQ=c(0.025,0.5,0.975),iA=iA,st=startsByYr,bse=bse)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"iA",iA,"bse",bse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","sQ"))

numObsYrs=c(20);startsByYr = 25;J=3; iA=0; aS=4; bse=7;sS=c(0,1,2)
scns=expand.grid(P=numObsYrs,J=J,sQ=c(0.025,0.5,0.975),iA=iA,aS=4,
                 st=startsByYr,bse=bse,sS=sS)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"aS",aS,"bse",bse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("sS","sQ"))

numObsYrs=c(20);startsByYr = 25;J=3; iA=0; aS=4; bse=1;sS=c(0,1,2)
scns=expand.grid(P=numObsYrs,J=J,sQ=c(0.025,0.5,0.975),iA=iA,aS=4,
                 st=startsByYr,bse=bse,sS=sS)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"aS",aS,"bse",bse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("sS","sQ"))

#Results
#ssv=0.8696: too much variability when no info available
#sse=0.1: too much variability when no info available
#lse=1: too much constraint, can't adjust mean even with lots of local data
#startsByYr: reducing sample size increases variation in observations among years. More difficult to infer differences or narrow CIs even with many years of data.

###############
#Step 2: confirm appropriate prior variability in recruitment intercept using minimal (2) observed data points & 0 fire/anthro covariates. Controlled by priors on l.Saf, phi and sig.Saf.
#################
#source("CaribouDemoFns.R")
numObsYrs=c(2);startsByYr = 25
scns=expand.grid(P=numObsYrs,rQ=c(0.5),st=startsByYr)
scResults = runScnSet(scns,eParsIn,simBig)
print(plotRes(scResults$rr.summary.all, "Recruitment",obs=scResults$obsRec.all,
              lowBound=0,simRange=scResults$sim.all,facetVars=c("P","rQ")))

numObsYrs=c(2,5,20);startsByYr = 25;lre=1
scns=expand.grid(P=numObsYrs,rQ=c(0.025,0.5,0.975),st=startsByYr,lre=lre)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"lre",lre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","rQ"))

numObsYrs=c(2,5,20);startsByYr = 25;sre=0.1
scns=expand.grid(P=numObsYrs,rQ=c(0.025,0.5,0.975),st=startsByYr,sre=sre)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"sre",sre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","rQ"))

numObsYrs=c(2,5,20);startsByYr = 25
scns=expand.grid(P=numObsYrs,rQ=c(0.025,0.5,0.975),st=startsByYr)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","rQ"))

numObsYrs=c(2,5,20);startsByYr = 25; iA=90
scns=expand.grid(P=numObsYrs,rQ=c(0.025,0.5,0.975),iA=iA,st=startsByYr)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"iA",iA)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","rQ"))

numObsYrs=c(2,5,20);startsByYr = 25; iA=90; bre=1
scns=expand.grid(P=numObsYrs,rQ=c(0.025,0.5,0.975),iA=iA,st=startsByYr,bre=bre)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"iA",iA,"bre",bre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","rQ"))

numObsYrs=c(20);startsByYr = 25;J=2; iA=0; aS=4; bre=7;rS=c(0,1,2)
scns=expand.grid(P=numObsYrs,J=J,rQ=c(0.025,0.5,0.975),iA=iA,aS=4,st=startsByYr,bre=bre,rS=rS)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"aS",aS,"bre",bre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("rS","rQ"))

numObsYrs=c(20);startsByYr = 25;J=2; iA=0; aS=4; bre=1;rS=c(0,1,2)
scns=expand.grid(P=numObsYrs,J=J,rQ=c(0.025,0.5,0.975),iA=iA,aS=4,st=startsByYr,bre=bre,rS=rS)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"aS",aS,"bre",bre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("rS","rQ"))

#TO DO: extend anthro logic to fire.
#TO DO: make warning fn from example above.
#Simulating reference range of variability takes some time, so provide option in UI to skip this check.
#But default should be to complain about parameter combos that give answers outside the range of what has been observed across the country.
#"Warning: for (anthro/fire) x,x,etc fitted 95% CI for (rec/surv/lambda) does not overlap range of variation simulated from national model.
#Note: could reduce time for simulating reference range by saving results - they generally don't change.

#Results
#sre=0.1: not enough variability when no info available
#lre=1: too much constraint, can't adjust mean even with lots of local data
#startsByYr: reducing sample size increases variation in observations among years. More difficult to infer differences or narrow CIs even with many years of data.

#TO DO: put lambda=1 on lambda plots

#################
#TO DO: show effect of changing uncertainty about anthro slope when no previous anthro
#source("CaribouDemoFns.R")
anthroSlope=0 #%change in anthro per year
anthroSlopeFuture=4
numObsYrs = c(2,5,20) # number of years of data
numProjectYrs = 20 #number of years to project
startsByYr = 25; iA=0
bre=1
scns=expand.grid(P=numObsYrs,J=numProjectYrs,rQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bre=bre)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bre",bre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","rQ"))
#TO DO: this doesn't make sense. Figure out what is wrong...

bre=7
scns=expand.grid(P=numObsYrs,J=numProjectYrs,rQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bre=bre)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bre",bre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","rQ"))

bre=3
scns=expand.grid(P=numObsYrs,J=numProjectYrs,rQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bre=bre)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bre",bre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","rQ"))


anthroSlope=2 #%change in anthro per year
anthroSlopeFuture=2
numProjectYrs = 20 #number of years to project
startsByYr = 25; iA=0
numObsYrs=20
bre=3; recSlopeMultiplier=c(0,1,2)
scns=expand.grid(P=numObsYrs,J=numProjectYrs,rQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bre=bre,rS=recSlopeMultiplier)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bre",bre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("rS","rQ"))

bre=1; recSlopeMultiplier=c(0,1,2)
scns=expand.grid(P=numObsYrs,J=numProjectYrs,rQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bre=bre,rS=recSlopeMultiplier)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bre",bre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("rS","rQ"))

bre=7; recSlopeMultiplier=c(0,1,2)
scns=expand.grid(P=numObsYrs,J=numProjectYrs,rQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bre=bre,rS=recSlopeMultiplier)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("rQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bre",bre)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("rS","rQ"))

#Result: bre=3 is a reasonable compromise. Not too much variation in absence of info,
#and vague enough to allow possibility of some other relationship.

############
#TO DO: show effect of changing uncertainty about anthro slope when no previous anthro
#source("CaribouDemoFns.R")
anthroSlope=0 #%change in anthro per year
anthroSlopeFuture=4
numObsYrs = c(2,5,20) # number of years of data
numProjectYrs = 20 #number of years to project
startsByYr = 25; iA=0
bse=1
scns=expand.grid(P=numObsYrs,J=numProjectYrs,sQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bse=bse)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bse",bse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","sQ"))

bse=5
scns=expand.grid(P=numObsYrs,J=numProjectYrs,sQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bse=bse)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bse",bse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","sQ"))

bse=7
scns=expand.grid(P=numObsYrs,J=numProjectYrs,sQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bse=bse)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bse",bse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("P","sQ"))

anthroSlope=2 #%change in anthro per year
anthroSlopeFuture=2
numProjectYrs = 20 #number of years to project
startsByYr = 25; iA=0
numObsYrs=20
bse=5; sefSlopeMultiplier=c(0,1,2)
scns=expand.grid(P=numObsYrs,J=numProjectYrs,sQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bse=bse,sS=sefSlopeMultiplier)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bse",bse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("sS","sQ"))

bse=1; sefSlopeMultiplier=c(0,1,2)
scns=expand.grid(P=numObsYrs,J=numProjectYrs,sQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bse=bse,sS=sefSlopeMultiplier)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bse",bse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("sS","sQ"))

bse=7; sefSlopeMultiplier=c(0,1,2)
scns=expand.grid(P=numObsYrs,J=numProjectYrs,sQ=c(0.025,0.5,0.975),
                 iA=iA,aS=anthroSlope,aSf=anthroSlopeFuture,
                 st=startsByYr,bse=bse,sS=sefSlopeMultiplier)
scResults = runScnSet(scns,eParsIn,simBig)
addBit = paste0("sQStarts",startsByYr,"aS",anthroSlope,"aSf",anthroSlopeFuture,"bse",bse)
makeInterceptPlots(scResults,addBit=addBit,facetVars=c("sS","sQ"))

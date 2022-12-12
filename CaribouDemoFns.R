runRMModel<-function(survData="simSurvData.csv",ageRatio.herd="simAgeRatio.csv",
                     disturbance="simDisturbance.csv",betaPriors="default",
                     startYear = 1998, endYear = 2018, Nchains = 4,Niter = 15000,Nburn = 10000,Nthin = 2,N0=1000,
                     survAnalysisMethod = "KaplanMeier",
                     inpFixed=list()){
  #survData=oo$simSurvObs;ageRatio.herd=oo$ageRatioOut;disturbance=oo$simDisturbance;
  #betaPriors=betaPriors;startYear = minYr;endYear=maxYr;N0=cs$N0;survAnalysisMethod = "KaplanMeier"
  #Nchains = 2;Niter = 20000;Nburn = 10000;Nthin = 1;inpFixed=list()

  # combine defaults in function with inputs from input list
  inputArgs = c("survData","ageRatio.herd","disturbance","startYear", "endYear",
                   "Nchains","Niter","Nburn","Nthin","N0","survAnalysisMethod")
  addArgs<-inputArgs#setdiff(inputArgs,names(inp))
  inp=list()
  for (a in addArgs){
    if(is.element(a,names(inpFixed))){
      inp[[a]]=inpFixed[[a]]
    }else{
      inp[[a]]=eval(parse(text=a))
    }
  }

  if(betaPriors[[1]]=="default"){
    betaPriors=getPriors()
  }

  #Run model
  if(is.character(inp$ageRatio.herd)){
    ageRatio.herd2 <- read.csv(paste0("tabs/",inp$ageRatio.herd),header=T)
    ageRatio.herd2$X=NULL
  }else{
    ageRatio.herd2 <- ageRatio.herd
  }
  if(is.character(inp$disturbance)){
    disturbance=read.csv(paste0("tabs/",inp$disturbance))
    disturbance$X = NULL
  }else{
    disturbance=inp$disturbance
  }
  disturbance=merge(data.frame(Year=seq(inp$startYear,inp$endYear)),disturbance,all.x=T)

  if(is.character(inp$survData)){
    survData <- read.csv(paste0("tabs/",inp$survData), header=T)
    survData$X=NULL
  }else{
    survData=inp$survData
  }
  compData <- ageRatio.herd2

  data=survData

  #check that year range is within data - model will run either way
  data$Year<-as.numeric(data$Year)

  #check that year range is within data - model will run either way
  if(inp$endYear > max(data$Year) | inp$startYear < min(data$Year))
  { warning(c("requested year range does not match survival data",c(" year range:","  ",min(data$Year)," - ",max(data$Year))))}

  data<-subset(data, data$Year <= inp$endYear & data$Year >= inp$startYear)
  data$id<-factor(data$id)

  test1 <- length(c(inp$startYear:inp$endYear))
  test2 <- length(unique(data$Year))

  #check that year range is within data - model will run either way
  if(test1 > test2)
  { warning(c("missing years of survival data.","Start model at beginning of consecutive years."," ",c("Years of survival data:","  ",list(sort(unique(data$Year))))))}

  data<-data[order(data$exit),]
  list_data1 <- split(data, data$Year)
  nYears<-length(levels(as.factor(data$Year)))
  n.ind<-numeric(nYears)

  for(i in 1:nYears){
    n.ind[i]<-length(list_data1[[i]]$exit)
  }

  #check that year range is within data - model will run either way
  if(any(n.ind < 20))
  { warning(c("warning, low sample size of adult females in at least one year"))}

  #get KM estimates to use for adult female survival

  #Note: biased results from years with <12 months of observations.
  #And problems with adding animals part way through the year, so omitting those
  dSubset =subset(data,enter==0)#;dSubset=subset(dSubset,!((exit<12)&(event==0)))
  if(nrow(dSubset)==0){
    stop("Years with less than 12 months of collar data are omitted from survival analysis. Please ensure there is 12 months of collar data in at least one year.")

  }
  if(inp$survAnalysisMethod=="KaplanMeier"){
    survData<-getKMSurvivalEstimates(dSubset)
    #omitting years with less than 12 months of observations of collared animals
    sumDat <-dSubset %>%
      group_by(Year) %>%
      summarise(minEnter = min(enter),maxExit=max(exit))
    includeYrs = subset(sumDat,minEnter==0&maxExit==12)
    survData$Year = as.numeric(gsub("as.factor(Year)=","",as.character(survData$Var1),fixed=T))
    survData <- merge(survData,includeYrs)
    if(nrow(survData)==0){
      stop("Years with less than 12 months of collar data are omitted from survival analysis. Please ensure there is 12 months of collar data in at least one year.")
    }

    if(any(survData$surv==1)){
      # which years does survival equal 1
      survOne=which(unlist(lapply(split(dSubset, dSubset$Year), function(x) sum(x$event)))==0)
      yearsOne=as.numeric(names(survOne))
      data.sub=data[data$Year %in% yearsOne,]
      nriskYears=data.frame(with(data.sub, table(Year)))

      #Specify model
      sink("binLik.txt")
      cat("
	   model{
	     for(i in 1:N){
	      lived[i] ~ dbin(s[i], atrisk[i])
	      s[i] ~ dbeta(1,1) # vague prior
	      }
	     }
	      ",fill = TRUE)
      sink()

      data1=list(N=nrow(nriskYears), lived=nriskYears$Freq, atrisk=nriskYears$Freq)
      params=c("s")
      inits=function(){
        list(s=runif(nrow(nriskYears), 0.80, 0.99))
      }

      # run model in JAGS
      out1=jags(data=data1, inits=inits, parameters.to.save=params,
                model.file="binLik.txt",n.chains=2, n.iter=5000,
                n.burnin=1000, n.thin=1)
      # create standard deviation variable from survData$tau above
      survData$tau <- 1/(survData$se^2)
      survData$tau[survOne]=1/(out1$BUGSoutput$sd$s^2)
    }else

      survData$tau <- 1/(survData$se^2)

    surv_id = which(survData$surv!=1)
    nSurv = length(surv_id)
    survData$Var1=as.character(survData$Var1)
  }else{
    if(inp$survAnalysisMethod=="exp"){
      #parametric exponential survival model
      survData = dSubset
      survData$t.to.death = survData$exit/12
      survData$t.to.death[!survData$event]=NA
      survData$t.cen=survData$exit/12
      survData$t.cen[survData$event]=0
    }else{
      dExpand <- apply(subset(dSubset,select=c(id,Year,event,enter,exit)),1,expandSurvivalRecord)
      survData <- do.call(rbind,dExpand)
    }
  }

  #split data into calf and cow recruitment data
  calf.cnt<-subset(compData, compData$Class=="calf")
  calf.cnt$Class<-factor(calf.cnt$Class)
  cow.cnt<-subset(compData, compData$Class=="cow")
  cow.cnt$Class<-factor(cow.cnt$Class)

  #deal with missing years of data between year ranges
  Years2<-data.frame(sort(unique(data$Year)))
  names(Years2) <- "Year"
  y1<-merge(Years2, calf.cnt, by="Year", all=TRUE)
  data3<-y1[,1:3]

  y2<-merge(Years2, cow.cnt, by="Year", all=TRUE)
  data4<-y2[,1:3]
  data3$Name <- rep(unique(compData$Name), length(data3$Count))
  data4$Name  <- rep(unique(compData$Name), length(data4$Count))

  data3$Count <- ifelse(data3$Count>data4$Count,NA,data3$Count)
  data4$Count <- ifelse(data3$Count>data4$Count,NA,data4$Count)

  xCalf <- which(is.na(data3$Count)==T)
  xCow <- which(is.na(data4$Count)==T)
  Years4 <- levels(as.factor(data$Year))[xCalf]

  if(any(is.na(data3$Count)==T))
  {
    warning("missing composition data; missing years:"," ",list(Years4))}

  t.pred = max(inp$endYear-max(data3$Year),0)

  if(t.pred>0){
    survAddBit = survData[1,]
    if(inp$survAnalysisMethod=="KaplanMeier"){
      survAddBit[1,]=NA;survAddBit$Var1=NULL;survAddBit=merge(survAddBit,data.frame(Var1=seq(1:t.pred)))
    }else{
      survAddBit[1:ncol(survAddBit)]=NA;survAddBit$Year=NULL;survAddBit = merge(survAddBit,data.frame(Year=(max(survData$Year)+seq(1:t.pred))))
    }

    survDatat=rbind(survData,survAddBit)

    dat3Bit = data3[1,]
    dat3Bit[,2:3]=NA;dat3Bit$Year=NULL;dat3Bit=merge(dat3Bit,data.frame(Year=seq(max(data3$Year)+1,max(data3$Year)+t.pred)))
    data3t=rbind(data3,dat3Bit)

    data4t=rbind(data4,dat3Bit)

  }else{
    survDatat=survData
    data3t=data3
    data4t=data4
  }

  #for beta model, remove tau, adjust Surv: (survDatat$surv*45+0.5)/46
  #survDatat$surv,tau=survDatat$tau
  #               phi.Saf.Prior1=betaPriors$phi.Saf.Prior1,phi.Saf.Prior2=betaPriors$phi.Saf.Prior2,
  if(inp$survAnalysisMethod=="KaplanMeier"){
    sp.data=list(Surv=survDatat$surv,tau=survDatat$tau,anthro=disturbance$Anthro,fire=disturbance$fire_excl_anthro,
                 beta.Saf.Prior1=betaPriors$beta.Saf.Prior1,beta.Saf.Prior2=betaPriors$beta.Saf.Prior2,
                 beta.Rec.anthro.Prior1=betaPriors$beta.Rec.anthro.Prior1,beta.Rec.anthro.Prior2=betaPriors$beta.Rec.anthro.Prior2,
                 beta.Rec.fire.Prior1=betaPriors$beta.Rec.fire.Prior1,beta.Rec.fire.Prior2=betaPriors$beta.Rec.fire.Prior2,
                 l.Saf.Prior1 = betaPriors$l.Saf.Prior1,l.Saf.Prior2 = betaPriors$l.Saf.Prior2,
                 l.R.Prior1 = betaPriors$l.R.Prior1,l.R.Prior2 = betaPriors$l.R.Prior2,
                 sig.Saf.Prior1=betaPriors$sig.Saf.Prior1,sig.Saf.Prior2 = betaPriors$sig.Saf.Prior2,
                 sig.R.Prior1 = betaPriors$sig.R.Prior1,sig.R.Prior2=betaPriors$sig.R.Prior2,
                 Ninit=inp$N0,
                 nCounts=length(which(is.na(data3t$Count)==FALSE)),count_id=which(is.na(data3t$Count)==FALSE),
                 nSurvs=length(which(is.na(survDatat$surv)==FALSE)),surv_id=which(is.na(survDatat$surv)==FALSE),
                 nYears=nYears+t.pred,calves=round(data3t$Count),CountAntlerless=round(data4t$Count))

    sp.params <- c("S.annual.KM" ,"R", "Rfemale", "pop.growth","fpop.size","var.R.real","l.R","l.Saf","beta.Rec.anthro","beta.Rec.fire","beta.Saf")

    rr.surv <- jags(data=sp.data, parameters.to.save=sp.params,
                    model.file="INF_saf_INF_r_JAGS_extentTime.txt",
                    n.chains=inp$Nchains, n.iter=inp$Niter, n.burnin=inp$Nburn, n.thin=inp$Nthin)
  }else{
    if(inp$survAnalysisMethod=="exp"){
      sp.data=list(t.to.death=survDatat$t.to.death,t.cen=survDatat$t.cen,survYr = survDatat$Year-inp$startYear,anthro=disturbance$Anthro,fire=disturbance$fire_excl_anthro,
                   beta.Saf.Prior1=betaPriors$beta.Saf.Prior1,beta.Saf.Prior2=betaPriors$beta.Saf.Prior2,
                   beta.Rec.anthro.Prior1=betaPriors$beta.Rec.anthro.Prior1,beta.Rec.anthro.Prior2=betaPriors$beta.Rec.anthro.Prior2,
                   beta.Rec.fire.Prior1=betaPriors$beta.Rec.fire.Prior1,beta.Rec.fire.Prior2=betaPriors$beta.Rec.fire.Prior2,
                   l.Saf.Prior1 = betaPriors$l.Saf.Prior1,l.Saf.Prior2 = betaPriors$l.Saf.Prior2,
                   l.R.Prior1 = betaPriors$l.R.Prior1,l.R.Prior2 = betaPriors$l.R.Prior2,
                   sig.Saf.Prior1=betaPriors$sig.Saf.Prior1,sig.Saf.Prior2 = betaPriors$sig.Saf.Prior2,
                   sig.R.Prior1 = betaPriors$sig.R.Prior1,sig.R.Prior2=betaPriors$sig.R.Prior2,
                   Ninit=inp$N0,
                   nCounts=length(which(is.na(data3t$Count)==FALSE)),count_id=which(is.na(data3t$Count)==FALSE),
                   nSurvs=length(which(is.na(survDatat[,1])==FALSE)),surv_id=which(is.na(survDatat$Year)==FALSE),
                   nYears=nYears+t.pred,calves=round(data3t$Count),CountAntlerless=round(data4t$Count))

      sp.params <- c("S.annual.KM" ,"R", "Rfemale", "pop.growth","fpop.size","var.R.real","l.R","l.Saf","beta.Rec.anthro","beta.Rec.fire","beta.Saf")

      rr.surv <- jags(data=sp.data, parameters.to.save=sp.params,
                      model.file="INF_saf_INF_r_JAGS_expSurvival.txt",
                      n.chains=inp$Nchains, n.iter=inp$Niter, n.burnin=inp$Nburn, n.thin=inp$Nthin)

    }else{
      sp.data=list(surv=survDatat[,1:13],survYr = survDatat$Year-inp$startYear+1,anthro=disturbance$Anthro,fire=disturbance$fire_excl_anthro,
                   beta.Saf.Prior1=betaPriors$beta.Saf.Prior1,beta.Saf.Prior2=betaPriors$beta.Saf.Prior2,
                   beta.Rec.anthro.Prior1=betaPriors$beta.Rec.anthro.Prior1,beta.Rec.anthro.Prior2=betaPriors$beta.Rec.anthro.Prior2,
                   beta.Rec.fire.Prior1=betaPriors$beta.Rec.fire.Prior1,beta.Rec.fire.Prior2=betaPriors$beta.Rec.fire.Prior2,
                   l.Saf.Prior1 = betaPriors$l.Saf.Prior1,l.Saf.Prior2 = betaPriors$l.Saf.Prior2,
                   l.R.Prior1 = betaPriors$l.R.Prior1,l.R.Prior2 = betaPriors$l.R.Prior2,
                   sig.Saf.Prior1=betaPriors$sig.Saf.Prior1,sig.Saf.Prior2 = betaPriors$sig.Saf.Prior2,
                   sig.R.Prior1 = betaPriors$sig.R.Prior1,sig.R.Prior2=betaPriors$sig.R.Prior2,
                   Ninit=inp$N0,
                   nCounts=length(which(is.na(data3t$Count)==FALSE)),count_id=which(is.na(data3t$Count)==FALSE),
                   nSurvs=length(which(is.na(survDatat[,1])==FALSE)),surv_id=which(is.na(survDatat$Year)==FALSE),
                   nYears=nYears+t.pred,calves=round(data3t$Count),CountAntlerless=round(data4t$Count))

      sp.params <- c("S.annual.KM" ,"R", "Rfemale", "pop.growth","fpop.size","var.R.real","l.R","l.Saf","beta.Rec.anthro","beta.Rec.fire","beta.Saf")

      rr.surv <- jags(data=sp.data, parameters.to.save=sp.params,
                      model.file="INF_saf_INF_r_JAGS_expSurvival2.txt",
                      n.chains=inp$Nchains, n.iter=inp$Niter, n.burnin=inp$Nburn, n.thin=inp$Nthin)

    }
  }

  return(list(result=rr.surv,survInput=survDatat))
}

getKMSurvivalEstimates<-function(dSubset){

  sModel = survfit(Surv(enter, exit, event)~as.factor(Year), conf.type = "log-log",data=dSubset)
  reg.out<-summary(sModel)
  if(0){
    check= data.frame(strata=reg.out$strata,survival= reg.out$surv,time=reg.out$time)
    check$type="est"
    check$Year = as.numeric(gsub("as.factor(Year)=","",as.character(check$strata),fixed=T))
    check$strata=NULL
    tt = subset(oo$exData,select=c(Year,survival))
    tt$time=12;tt$type="truth"
    check=rbind(check,tt)
    base=ggplot(check,aes(x=time,y=survival,colour=type,shape=type))+geom_point()+facet_wrap(~Year)
    print(base)
  }

  if(is.null(reg.out$strata)){reg.out$strata=as.character(dSubset$Year[1])}
  data5<-data.frame(reg.out$strata, reg.out$surv)
  data.se<-data.frame(reg.out$strata, reg.out$std.err)
  data.l<-data.frame(reg.out$strata, reg.out$lower)
  data.u<-data.frame(reg.out$strata, reg.out$upper)
  data6<-data.frame(table(data5$reg.out.strata))
  num<-data6$Freq
  if(class(data5$reg.out.strata)=="character"){data5$reg.out.strata=as.factor(data5$reg.out.strata)}
  levs<-data.frame(levels(data5$reg.out.strata))
  names(levs)<-c("reg.out.strata")
  data7<-subset(data6, Freq>0)
  survs<-numeric(length(data7$Freq))
  se<-numeric(length(data7$Freq))
  lower<-numeric(length(data7$Freq))
  upper<-numeric(length(data7$Freq))
  index<-cumsum(data7$Freq)

  for(i in 1:length(survs)){
    survs[i]<-data5$reg.out.surv[index[i]]
    se[i] <- data.se$reg.out.std.err[index[i]]
    lower[i]<-data.l$reg.out.lower[index[i]]
    upper[i]<-data.u$reg.out.upper[index[i]]
  }

  indexS<-which(data6$Freq!=0)
  data6$surv<-numeric(length(data6$Freq))
  data6$se<-numeric(length(data6$Freq))
  data6$lower<-numeric(length(data6$Freq))
  data6$upper<-numeric(length(data6$Freq))
  data6$surv<-ifelse(data6$Freq<1, 1, NA )
  data6$lower<-ifelse(data6$Freq<1, NA, 0)
  data6$upper<-ifelse(data6$Freq<1, NA, 0)

  for(i in 1:length(data6$Freq)){
    data6$surv[indexS[i]]<-survs[i]
    data6$se[indexS[i]]<-se[i]
    data6$lower[indexS[i]]<-lower[i]
    data6$upper[indexS[i]]<-upper[i]
  }

  survData <- data6

  return(survData)
}

#install and load required packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

expandSurvivalRecord<-function(crow){
  #crow=subset(dSubset,exit==12)[1,]
  crow=as.numeric(crow)
  mnths= c(rep(1,crow[5]-crow[4]),rep(!crow[3],12-crow[5]+1))
  mnths = data.frame(t(mnths))
  mnths$Year=crow[2]
  mnths$enter=crow[4]
  mnths$exit=crow[5]
  mnths$event=crow[3]
  return(mnths)
}

getPriors<-function(modifiers=NULL,
                    expectMods = list(survivalModelNumber="M1",
                                      recruitmentModelNumber="M4",
                                      bre=3,
                                      bse=3,
                                      lse = 5,
                                      sse = 0.08696*0.4,
                                      ssv = 0.01,
                                      lre=9,
                                      sre = 0.46,
                                      srv = 0.22),
                    popGrowthTable = popGrowthTableJohnsonECCC,
                    modVer = "Johnson",returnValues=T){
  #modifiers=cs


  if(is.null(modifiers)){modifiers=expectMods}
  addMods<-setdiff(names(expectMods),names(modifiers))
  for(i in addMods){
    modifiers[[i]]=expectMods[[i]]
  }

  popGrowthPars <- demographicCoefficients(
    2,
    modelVersion = modVer,
    survivalModelNumber = modifiers$survivalModelNumber,
    recruitmentModelNumber = modifiers$recruitmentModelNumber,
    populationGrowthTable = popGrowthTable
  )

  rPriorCoefs <- popGrowthPars$coefSamples_Recruitment$coefValues
  rPriorStdErrs <- popGrowthPars$coefSamples_Recruitment$coefStdErrs
  sPriorCoefs <- popGrowthPars$coefSamples_Survival$coefValues
  sPriorStdErrs <- popGrowthPars$coefSamples_Survival$coefStdErrs

  # check for variables not programmed in jags
  diff_pars <- setdiff(names(rPriorCoefs),
                       c("Intercept", "Precision", "Anthro", "fire_excl_anthro"))
  if(length(diff_pars) > 0){
    stop("The recruitment model contains unrecognized coefficients: ",
         paste0(diff_pars, collapse = ", "))
  }

  diff_pars <- setdiff(names(sPriorCoefs),
                       c("Intercept", "Precision", "Anthro", "fire_excl_anthro"))
  if(length(diff_pars) > 0){
    stop("The survival model contains unrecognized coefficients: ",
         paste0(diff_pars, collapse = ", "))
  }

  if(returnValues){
    betaPriors = list(l.R.Prior1 = rPriorCoefs$Intercept,
                      l.R.Prior2 = rPriorStdErrs$Intercept*modifiers$lre,
                      beta.Rec.anthro.Prior1= rPriorCoefs$Anthro,
                      beta.Rec.anthro.Prior2 = rPriorStdErrs$Anthro*modifiers$bre,
                      beta.Rec.fire.Prior1 = rPriorCoefs$fire_excl_anthro,
                      beta.Rec.fire.Prior2 = rPriorStdErrs$fire_excl_anthro,
                      sig.R.Prior1=modifiers$sre,
                      sig.R.Prior2 =modifiers$srv,
                      l.Saf.Prior1 = sPriorCoefs$Intercept,
                      l.Saf.Prior2 = sPriorStdErrs$Intercept*modifiers$lse,
                      beta.Saf.Prior1=sPriorCoefs$Anthro,
                      beta.Saf.Prior2=sPriorStdErrs$Anthro*modifiers$bse,
                      sig.Saf.Prior1 = modifiers$sse,
                      sig.Saf.Prior2=modifiers$ssv
    )

    # replace NULL values with 0 for when anthro or fire is not included
    betaPriors <- lapply(betaPriors, function(x){
      if(is.null(x) || length(x) == 0){
        1e-10
      } else {
        x
      }
    })
  }else{
    betaPriors = list(l.R.Prior1 = rPriorCoefs$Intercept,
                      l.R.Prior2 = paste0(round(rPriorStdErrs$Intercept,4),"*",modifiers$lre),
                      beta.Rec.anthro.Prior1= rPriorCoefs$Anthro,
                      beta.Rec.anthro.Prior2 = paste0(round(rPriorStdErrs$Anthro,4),"*",modifiers$bre),
                      beta.Rec.fire.Prior1 = rPriorCoefs$fire_excl_anthro,
                      beta.Rec.fire.Prior2 = rPriorStdErrs$fire_excl_anthro,
                      sig.R.Prior1=modifiers$sre,
                      sig.R.Prior2 =modifiers$srv,
                      l.Saf.Prior1 = sPriorCoefs$Intercept,
                      l.Saf.Prior2 = paste0(round(sPriorStdErrs$Intercept,4),"*",modifiers$lse),
                      beta.Saf.Prior1=sPriorCoefs$Anthro,
                      beta.Saf.Prior2=paste0(round(sPriorStdErrs$Anthro,4),"*",modifiers$bse),
                      sig.Saf.Prior1 = modifiers$sse,
                      sig.Saf.Prior2=modifiers$ssv
    )

  }
  return(betaPriors)

}

simCovariates<-function(initAnthro,initFire,numYears,anthroSlope,anthroSlopeFuture,futureStep,fireSlope=0){
  covInit= data.frame(Anthro=initAnthro,fire_excl_anthro=initFire)
  for (t in 1:numYears) {
    #t=1s
    if(t==1){
      cov <- covInit
    }else{
      cov<- covPrev
    }
    if(t>1){
      if(t>=futureStep){
        cov$Anthro <- pmin(100,covPrev$Anthro + anthroSlopeFuture)
      }else{
        cov$Anthro <- pmin(100,covPrev$Anthro + anthroSlope)
      }
    }
    covPrev=cov
    if(fireSlope==0){
      cov$fire_excl_anthro<-pmax(0,rnorm(1,cov$fire_excl_anthro,0.0001))
    }else{
      cov$fire_excl_anthro<-pmin(100,cov$fire_excl_anthro+fireSlope*(t-1))
    }
    cov$Total_dist = cov$Anthro + cov$fire_excl_anthro
    cov$time=t
    if(t==1){
      covariates<-cov
    }else{
      covariates<-rbind(covariates,cov)
    }

  }
  return(covariates)

}

simTrajectory<-function(numYears,covariates,survivalModelNumber = "M1",recruitmentModelNumber = "M4",
                        popGrowthTable = popGrowthTableJohnsonECCC,
                        recSlopeMultiplier=1,sefSlopeMultiplier=1,recQuantile=0.5,sefQuantile=0.5,
                        stepLength=1,N0=1000){
  # survivalModelNumber = "M1";recruitmentModelNumber = "M4";
  #recSlopeMultiplier=1;sefSlopeMultiplier=1;recQuantile=0.5;sefQuantile=0.5
  #stepLength=1;N0=1000

  #alter coefficients
  growthTab<-popGrowthTable
  growthTab$Value[(growthTab$Coefficient=="Anthro")&(growthTab$responseVariable=="recruitment")]=recSlopeMultiplier*growthTab$Value[(growthTab$Coefficient=="Anthro")&(growthTab$responseVariable=="recruitment")]
  growthTab$Value[(growthTab$Coefficient=="Anthro")&(growthTab$responseVariable=="femaleSurvival")]=sefSlopeMultiplier*growthTab$Value[(growthTab$Coefficient=="Anthro")&(growthTab$responseVariable=="femaleSurvival")]

  popGrowthParsSmall <- demographicCoefficients(
    2,
    modelVersion = "Johnson",
    survivalModelNumber = survivalModelNumber,
    recruitmentModelNumber = recruitmentModelNumber,
    populationGrowthTable = growthTab,
    useQuantiles = c(recQuantile,recQuantile)
  )
  #manually set quantiles for example population
  popGrowthParsSmall$coefSamples_Survival$quantiles=sefQuantile
  #TO DO: across full covariate range, note proportion of this scenario that is outside distribution of observations across the country

  # Only use precision if included in the table for this model number for both rec and surv
  usePrec <- "Precision" %in% names(popGrowthParsSmall$coefSamples_Survival$coefValues) &
    "Precision" %in% names(popGrowthParsSmall$coefSamples_Recruitment$coefValues)
  # at each time,  sample demographic rates and project, save results
  pars <- data.frame(N0 = N0)
  for (t in 1:numYears) {
    #t=1
    covs<-subset(covariates,time==t)

    rateSamples <- demographicRates(
      covTable = covs,
      popGrowthPars = popGrowthParsSmall,
      ignorePrecision = !usePrec,
      returnSample = TRUE
    )
    if (is.element("N", names(pars))) {
      pars <- subset(pars, select = c(replicate, N))
      names(pars)[names(pars) == "N"] <- "N0"
    }
    pars <- merge(pars, rateSamples)
    pars <- cbind(pars,
                  popGrowthJohnson(pars$N0,
                                   R_bar = pars$R_bar, S_bar = pars$S_bar,
                                   numSteps = stepLength, K=F, l_R = 1e-06,adjustR=T))

    # add results to output set
    fds <- subset(pars, select = c(replicate, Anthro,fire_excl_anthro, S_t, R_t, N, lambda,n_recruits,surviving_adFemales))
    fds$replicate <- as.numeric(gsub("V", "", fds$replicate))
    names(fds) <- c("Replicate", "Anthro","fire_excl_anthro", "survival", "recruitment", "N", "lambda","n_recruits","n_cows")
    fds$n_recruits=fds$n_recruits*2 #all calves, not just females
    fds <- pivot_longer(fds, !Replicate, names_to = "MetricTypeID", values_to = "Amount")
    fds$Timestep <- t * stepLength
    if (t == 1) {
      popMetrics <- fds
    } else {
      popMetrics <- rbind(popMetrics, fds)
    }
  }

  popMetrics$MetricTypeID <- as.character(popMetrics$MetricTypeID)
  popMetrics$Replicate <- paste0("x", popMetrics$Replicate)
  return(subset(popMetrics,Replicate=="x1"))

}

simCalfCowRatios<-function(cowCounts,minYr,exData){
  ageRatioOut=subset(cowCounts,(Year>=minYr),select=c(Year,Class,Count))#assume info from only one herd
  ageRatioOut=pivot_wider(ageRatioOut,id_cols=c("Year"),names_from="Class",values_from="Count")
  ageRatioOut=merge(ageRatioOut,subset(exData,select=c("Year","n_recruits","n_cows")))
  # rbinom needs n_recruits to be <= n_cows and n_cows not 0
  n_recs <- pmin(ageRatioOut$n_cows, ageRatioOut$n_recruits)
  ageRatioOut$calf= ifelse(ageRatioOut$n_cows == 0, 0,
                           rbinom(n=nrow(ageRatioOut),size=ageRatioOut$cow,
                                  prob=n_recs/ageRatioOut$n_cows))
  ageRatioOut=subset(ageRatioOut,select=c(Year,calf,cow))
  ageRatioOut=pivot_longer(ageRatioOut,cols=c(calf,cow),names_to="Class",values_to="Count")
  return(ageRatioOut)
}

simSurvivalData<-function(freqStartsByYear,exData,collarNumYears,collarOffTime,collarOnTime){
  #for simplicity, ignore variation in survival probability among months
  initYear = min(exData$Year)
  freqStartsByYear=subset(freqStartsByYear,Year>=initYear)
  survivalSeries = subset(exData,select=c(survival,Year))

  #freqStartsByYear$numStarts=10000;collarNumYears=2

  #if(initYear<min(survivalSeries$Year)){
  #  missingYrs = seq(initYear,min(survivalSeries$Year)-1)
  #  addBit = subset(survivalSeries,Year==min(survivalSeries$Year))
  #  addBit$Year=NULL;addBit = merge(addBit, data.frame(Year=missingYrs))
  #  survivalSeries=rbind(survivalSeries,addBit)
  #}

  animalID = 1
  started=F
  for(k in 1:nrow(freqStartsByYear)){
    #k =1
    if(is.na(freqStartsByYear$numStarts[k])|(freqStartsByYear$numStarts[k]<=0)){
      next
    }
    startYear = freqStartsByYear$Year[k]
    for(n in 1:freqStartsByYear$numStarts[k]){
      #n=1
      addS <- simSurvivalObs(animalID,startYear,collarNumYears,collarOffTime,collarOnTime,survivalSeries)
      animalID= animalID+1
      if(!started){
        simSurvObs<-addS
        started=T
      }else{
        simSurvObs<-rbind(simSurvObs,addS)
      }
    }
  }



  #1-sum(simSurvObs$event)/nrow(simSurvObs)
  #exData

  simSurvObs<-subset(simSurvObs,is.element(Year,exData$Year))
  simSurvObs<-simSurvObs[order(simSurvObs$Year),]

  addBit <- unique(subset(freqStartsByYear,select=setdiff(names(freqStartsByYear),c("numStarts",names(simSurvObs)))))
  if(nrow(addBit)>1){
    stop()
  } else if(nrow(addBit) > 0){
    simSurvObs<-merge(simSurvObs,addBit)
  }

  return(simSurvObs)
}

simSurvivalObs<-function(animalID,startYear,collarNumYears,collarOffTime,collarOnTime,survivalSeries){
  #animalID =  1; startYear = 2016
  for(i in startYear:min((startYear+collarNumYears-1),max(survivalSeries$Year))){
    #i = startYear
    surv=survivalSeries$survival[survivalSeries$Year==i]^(1/12)#monthly survival

    if(i==startYear){
      enter=collarOnTime-1
    }else{enter=0}

    if(i==(startYear+collarNumYears-1)){
      exit=collarOffTime
    }else{exit=12}

    event=0
    for(j in (enter+1):exit){
      die = rbinom(1,1,1-surv)
      if(die){
        event=1;exit=j
        break
      }
    }

    addBit = data.frame(id=animalID,Year=i,event=die,enter=enter,exit=exit)

    if(i==startYear){
      survObs = addBit
    }else{
      survObs = rbind(survObs,addBit)
    }
    if(die){break}
  }

  return(survObs)
}


# Tables ------------------------------------------------------------------

getSumStats <- function(param, rrSurvMod, startYear, endYear,doSummary=T){
  #param = "S.annual.KM"
  paramNames <- data.frame(parameter = c("S.annual.KM", "R", "Rfemale", "pop.growth","fpop.size",
                                     "meanAFsurv", "meanR", "meanRfemale",
                                     "medianLambda", "meanLambda"),
                           name = c("Adult female survival", "Recruitment",
                                    "Female-only recruitment", "Population growth rate","Female population size",
                                    "Mean adult female survival",
                                    "Mean recruitment", "Mean female recruitment",
                                    "Median population growth rate",
                                    "Mean population growth rate"))

  paramNames<-subset(paramNames, is.element(parameter,rrSurvMod$parameters.to.save))

  if(grepl("mean|median", param)){
    yr <- NA_integer_
  } else {
    yr <- startYear:endYear
  }

  if(!param %in% paramNames$parameter){
    stop("param ", param, "is not recognized\n",
         "accepted params are: ", paramNames$parameter)
  }

  if(doSummary){
    lower.cri <- apply(rrSurvMod$BUGSoutput$sims.list[[param]], 2,
                       function(x){quantile(x, 0.025)})
    upper.cri <- apply(rrSurvMod$BUGSoutput$sims.list[[param]], 2,
                       function(x){quantile(x, 0.975)})

    results <- data.frame(Year = yr,
                          Parameter = subset(paramNames, paramNames$parameter == param,
                                             select = name, drop = T),
                          Mean = round(rrSurvMod$BUGSoutput$mean[[param]],digits=3),
                          SD = round(rrSurvMod$BUGSoutput$sd[[param]],digits=3),
                          `Lower 95% CRI` = round(lower.cri, digits=3),
                          `Upper 95% CRI` = round(upper.cri,digits=3),
                          check.names = FALSE)
  }else{
    #rrSurvMod= result
    wideRes <- data.frame(rrSurvMod$BUGSoutput$sims.list[[param]])
    names(wideRes)<-yr

    results<-wideRes %>%pivot_longer(cols=names(wideRes),names_to="Year",values_to="Value")
    results$Year=as.numeric(results$Year)
    results$Parameter = subset(paramNames, paramNames$parameter == param,
                               select = name, drop = T)
  }
  return(results)
}

tabAllRes <- function(rrSurvMod, startYear, endYear,doSummary=T){
  #rrSurvMod=rr.surv;startYear= minYr;endYear= maxYr

  allParams <- c("S.annual.KM", "R", "Rfemale", "pop.growth","fpop.size",
                 "meanAFsurv", "meanR", "meanRfemale",
                 "medianLambda", "meanLambda")
  allParams<-allParams[is.element(allParams,rrSurvMod$parameters.to.save)]

  allResults <- lapply(allParams, getSumStats, rrSurvMod, startYear, endYear,doSummary=doSummary)

  allResults <- do.call(rbind, allResults)

  allResults <- allResults[order(allResults$Year),]
  allResults <- allResults[order(allResults$Parameter),]
  row.names(allResults) <- 1:length(allResults$Year)
  allResults
}

getParamsFromEacker<-function(path){
  ###############
  #Use Eacker example data for sample sizes in each year.
  survData=paste0(path,"/tte_caribouFEMALES.csv")
  ageRatio.herd=paste0(path,"/ageRatio.herd.csv")
  ageRatio.herd2 <- read.csv(ageRatio.herd,header=T)
  tte_caribou2 <- read.csv(survData, header=T)

  #need table of observed number of cows each year as input for simulating calf:cow ratios. Use Eaker as example.
  cowCounts = subset(ageRatio.herd2,Class=="cow")
  write.csv(cowCounts,"tabs/cowCounts.csv")
  yrRange = max(tte_caribou2$Year)-min(tte_caribou2$Year)

  #get survival sampling parameters from example data
  animalStarts<-subset(tte_caribou2,select=c(id,Year))%>%group_by(id)%>%summarise(startYear=min(Year))
  freqStartsByYear<- as.data.frame(table(animalStarts$startYear));names(freqStartsByYear)=c("Year","numStarts")
  freqStartsByYear$Year = as.numeric(as.character(freqStartsByYear$Year))
  sm = tte_caribou2%>%group_by(id)%>%summarize(startYear=min(Year),endYear = max(Year),numYears=max(Year)-min(Year),died=sum(event))
  collarNumYears = median(subset(sm,!died)$numYears) #for simplicity, pick a single number of years that collars remain on

  addInfo=unique(subset(tte_caribou2,select=c(HerdDescription,HerdCode,Range_ID,RangeDescription,RangeCode)))
  #TO DO: remove requirement for additional herd ID and range ID columns in UI code
  freqStartsByYear = merge(freqStartsByYear,addInfo)
    write.csv(freqStartsByYear,"tabs/freqStartsByYear.csv")

  offSet = subset(sm,!died,select=c(id,endYear));names(offSet)=c("id","Year")
  offSet = merge(offSet,tte_caribou2,all.x=T)
  collarOffTime = median(offSet$exit) # for simplicity, pick a single month that collars fall off

  onSet = subset(sm,select=c(id,startYear));names(onSet)=c("id","Year")
  onSet = merge(onSet,tte_caribou2,all.x=T)
  collarOnTime = median(onSet$enter) # for simplicity, pick a single month that collars are put on
  #TO DO: allow users to set freqStartsByYear, collarNumYears, collarOffTime, and collarOnTime as parameters OR
  #provide a file formatted as in Eacker from which these parameters can be derived.

  #freqStartsByYear$numStarts=30
  return(list(cowCounts=cowCounts,freqStartsByYear=freqStartsByYear,collarOnTime=collarOnTime,collarOffTime=collarOffTime,collarNumYears=collarNumYears))
}

fillDefaults <- function(scns = NULL,
                         defList = list(
                           iF = 0, iA = 0, aS = 0, aSf = 4,
                           rS = 1, sS = 1,
                           rQ = 0.5, sQ = 0.5, J = 20, P = 2, st = 25, N0 = 1000
                         ), curYear = 2018) {
  if (is.null(scns)) {
    scns <- as.data.frame(defList)
  } else {
    fillSet <- setdiff(names(defList), names(scns))

    for (i in fillSet) {
      scns[[i]] <- defList[[i]]
    }
  }
  scns$ID <- seq(1:nrow(scns))
  scns$label <- ""
  for (n in names(scns)[(length(names(scns)) - 1):1]) {
    scns$label <- paste0(scns$label, n, scns[[n]], "_")
  }

  if (!is.element("iYr", names(scns))) {
    scns$iYr <- curYear - scns$P + 1
  }
  return(scns)
}

getOutputTables<-function(result,startYear,endYear,survInput,oo,simBig,getKSDists){
  #result=out$result;startYear=minYr;endYear=maxYr;survInput=out$survInput;oo=oo;simBig=simLow

  #get summary info for plots
  rr.summary<-tabAllRes(result, startYear, endYear)

  if(!is.element("surv",names(survInput))){
    if(sum(survInput$event,na.rm=T)>0){
      obsSurv = getKMSurvivalEstimates(survInput)
    }else{
      obsSurv =unique(subset(survInput,!is.na(enter),select=c(Year)))
      obsSurv$surv=NA
      obsSurv$Var1=obsSurv$Year
    }
  }else{
    obsSurv<-survInput
  }

  obsSurv$Mean=obsSurv$surv
  obsSurv$Year = as.numeric(gsub("as.factor(Year)=","",obsSurv$Var1,fixed=T))
  obsSurv<- subset(obsSurv,Year>1000)

  obsSurv$parameter="Adult female survival"
  obsSurv$type="observed"

  trueSurv = subset(oo$exData,select=c(Year,survival))
  names(trueSurv)=c("Year","Mean")
  trueSurv$parameter="Adult female survival"
  trueSurv$type="true"

  obsRec=subset(oo$ageRatioOut,select=c(Year,Count,Class))
  obsRec <-pivot_wider(obsRec,id_cols=c("Year"),names_from="Class",values_from="Count")
  obsRec$Mean=obsRec$calf/obsRec$cow
  obsRec$parameter="Recruitment"
  obsRec$type="observed"

  trueRec = subset(oo$exData,select=c(Year,recruitment))
  names(trueRec)=c("Year","Mean")
  trueRec$parameter="Recruitment"
  trueRec$type="true"

  obsLam = subset(oo$exData,select=c(Year,lambda))
  names(obsLam)=c("Year","Mean")
  obsLam$parameter="Population growth rate"
  obsLam$type="true"

  obsSize = subset(oo$exData,select=c(Year,N))
  names(obsSize)=c("Year","Mean")
  obsSize$parameter="Female population size"
  obsSize$type = "true"

  obsAll=rbind(obsLam,obsSize,subset(obsRec,select=names(obsLam)),trueRec,subset(obsSurv,select=names(obsLam)),trueSurv)

  simBigO<-subset(simBig$summary,select=c(Anthro,Mean,lower,upper,parameter))
  names(simBigO)<-c("Anthro","Mean","Lower 95% CRI","Upper 95% CRI","parameter")

  # combine cs and simDisturbance and add to all output tables, nest params in a list
  dist_params <- merge(oo$simDisturbance, oo$cs)

  rr.summary=merge(rr.summary,dist_params)
  simBigO=merge(simBigO,dist_params)
  obsAll=merge(obsAll,dist_params)
  rr.summary.all = rr.summary
  sim.all = simBigO
  obs.all=obsAll

  if(getKSDists){
    #get Kolmogorov smirnov distance between samples at each point

    variables = unique(simBig$summary$parameter)
    anthroPts = unique(subset(rr.summary,select=c(Year,Anthro)))
    #TO DO: make this step faster
    bmSamples = tabAllRes(result, startYear, endYear,doSummary=F)
    bmSamples$type="local"

    simSamples=merge(anthroPts,simBig$samples);simSamples$Anthro=NULL
    simSamples$type = "national"

    allSamples=rbind(subset(bmSamples,is.element(Parameter,unique(simSamples$Parameter))),simSamples)

    ksDists <- allSamples %>% group_by(Year,Parameter) %>% group_modify(~ {getKSDist(.x$Value,.x$type)})
  }else{
    ksDists <-unique(subset(rr.summary,select=c(Year,Parameter)))
    ksDists$KSDistance=NA;ksDists$KSpvalue=NA

  }
  return(list(rr.summary.all=rr.summary.all,sim.all=sim.all,obs.all=obs.all,ksDists=ksDists))
}

getKSDist<-function(Value,type){
  #sampleBit=subset(allSamples,(Year==2017)&(Parameter==allSamples$Parameter[1]))
  #Value=sampleBit$Value;type=sampleBit$type

  if(length(Value[type=="national"])==0){
    out=data.frame(KSDistance=NA,KSpvalue=NA)
    return(out)
  }
  res = ks.test(Value[type=="local"],Value[type=="national"])

  out=data.frame(KSDistance=res$statistic,KSpvalue=res$p.value)
  return(out)
}

makeInterceptPlots<-function(scResults,addBit="",facetVars=c("P","sQ"),loopVars = NULL,
                             whichPlots=c("Adult female survival","Population growth rate","Recruitment","Female population size"),
                             survLow=0.6,type="png"){
  #facetVars=c("lse","sse");loopVars="ssv"

  if(!is.null(loopVars)){

    loopSet = unique(subset(scResults$rr.summary.all,select=loopVars))
    loopSet$dummy = 1

  }else{
    loopSet=data.frame(dummy=1)
  }


  for(l in 1:nrow(loopSet)){
    #l = 1
    crow=loopSet[l,]

    aa = ""
    for(n in names(crow)){
      if(n=="dummy"){next}
      aa=paste0(aa,n,crow[[n]])
    }

    addBitO = paste0(addBit,aa)

    if(is.element("Adult female survival",whichPlots)){
      if(type=="png"){
        png(here::here(paste0("figs/Surv",addBitO,".png")),
            height = 6, width = 7.48, units = "in",res=600)
      }else{pdf(paste0("figs/Surv",addBitO,".pdf"),width=10,height=7)}
      print(plotRes(merge(scResults$rr.summary.all,crow), "Adult female survival",obs=merge(scResults$obs.all,crow),
                    lowBound=survLow,simRange=merge(scResults$sim.all,crow),facetVars=facetVars))
      dev.off()
    }

    if(is.element("Population growth rate",whichPlots)){
      if(type=="png"){
        png(here::here(paste0("figs/Lambda",addBitO,".png")),
            height = 6, width = 7.48, units = "in",res=600)
      }else{pdf(paste0("figs/Lambda",addBitO,".pdf"),width=10,height=7)}
      print(plotRes(scResults$rr.summary.all, "Population growth rate",obs=scResults$obs.all,
                    lowBound=0,simRange=scResults$sim.all,facetVars=facetVars))
      dev.off()
    }

    if(is.element("Recruitment",whichPlots)){
      if(type=="png"){
        png(here::here(paste0("figs/Rec",addBitO,".png")),
            height = 6, width = 7.48, units = "in",res=600)
      }else{pdf(paste0("figs/Rec",addBitO,".pdf"),width=10,height=7)}
      print(plotRes(scResults$rr.summary.all, "Recruitment",obs=scResults$obs.all,
                    lowBound=0,simRange=scResults$sim.all,facetVars=facetVars))
      dev.off()
    }

    if(is.element("Female population size",whichPlots)){
      if(type=="png"){
        png(here::here(paste0("figs/FPOP",addBitO,".png")),
            height = 6, width = 7.48, units = "in",res=600)
      }else{pdf(paste0("figs/FPOP",addBitO,".pdf"),width=10,height=7)}
      print(plotRes(scResults$rr.summary.all, "Female population size",obs=scResults$obs.all,
                    lowBound=0,facetVars=facetVars))
      dev.off()
    }
  }
}

getSimsNational<-function(reps=500,N0=1000,Anthro=seq(0,100,by=1),fire_excl_anthro=0,quants=NULL,wdir=NULL){
  #reps=500;N0=500;Anthro=seq(0,100,by=2);fire_excl_anthro=0;quants=c(0.025,0.025)

  if(is.null(wdir)){wdir=getwd()}
  doSave <- FALSE

  check <- as.list(match.call());check$wdir=NULL

  if(length(check) == 1){

    if(file.exists(paste0(wdir,"/simsNational.rds"))){
      message("Using saved object")
      return(readRDS(paste0(wdir,"/simsNational.rds")))
    } else {
      message("Object will be saved for future use")
      doSave <- TRUE
    }

  }
  covTableObs <- expand.grid(Anthro = Anthro,
                             fire_excl_anthro = fire_excl_anthro)
  covTableObs$Total_dist = covTableObs$Anthro + covTableObs$fire_excl_anthro
  if(is.null(quants)){
    popGrowthPars <- demographicCoefficients(reps)
    rateSamplesAll <- demographicRates(covTable = covTableObs,popGrowthPars = popGrowthPars,returnSample=T,useQuantiles = F)
  }else{
    popGrowthPars <- demographicCoefficients(reps,useQuantiles = quants)
    rateSamplesAll <- demographicRates(covTable = covTableObs,popGrowthPars = popGrowthPars,returnSample=T)

  }
  pars <- merge(data.frame(N0 = N0), rateSamplesAll)
  pars <- cbind(pars,popGrowthJohnson(pars$N0,R_bar = pars$R_bar, S_bar = pars$S_bar,numSteps = 1,K=F,adjustR=T))
  simSurvBig<-pars %>% select(Anthro,S_t) %>% group_by(Anthro) %>% summarize(Mean=mean(S_t),lower=quantile(S_t, 0.025),upper=quantile(S_t,0.975))
  simSurvBig$parameter="Adult female survival"
  simRecBig<-pars %>% select(Anthro,R_t) %>% group_by(Anthro) %>% summarize(Mean=mean(R_t),lower=quantile(R_t, 0.025),upper=quantile(R_t,0.975))
  simRecBig$parameter="Recruitment"
  simLamBig<-pars %>% select(Anthro,lambda) %>% group_by(Anthro) %>% summarize(Mean=mean(lambda),lower=quantile(lambda, 0.025),upper=quantile(lambda,0.975))
  simLamBig$parameter="Population growth rate"
  simFpopBig<-pars %>% select(Anthro,N) %>% group_by(Anthro) %>% summarize(Mean=mean(N),lower=quantile(N, 0.025),upper=quantile(N,0.975))
  simFpopBig$parameter="Female population size"
  simBig =rbind(simSurvBig,simRecBig,simLamBig,simFpopBig)

  parsSelect = subset(pars,select=c(Anthro,S_t,R_t,lambda,N))
  names(parsSelect)=c("Anthro","Adult female survival","Recruitment","Population growth rate","Female population size")
  parsSelect = parsSelect %>% pivot_longer(!Anthro,names_to="Parameter",values_to="Value")

  simBig=list(summary=simBig,samples=parsSelect)

  if(doSave){
    saveRDS(simBig, paste0(wdir,"/simsNational.rds"))
  }
  return(simBig)
}

# test a popGrowthTable has the right format
testPopGrowthTable <- function(df) {
  # required columns
  missed_nms <- setdiff(c("responseVariable", "Coefficient", "Value"),
                        colnames(df))

  # need stdErr or CI
  if(! "StdErr" %in% colnames(df)){
    missed_nms <- c(missed_nms,
                    setdiff(c("lowerCI", "upperCI"), colnames(df)))
    df <- mutate(df, StdErr = NA_real_)
  } else if(!all(c("lowerCI", "upperCI") %in% colnames(df))){
    df <- mutate(df, lowerCI = NA_real_, upperCI = NA_real_)
  }

  if(length(missed_nms) > 0){
    stop("The model coefficient file loaded is missing the columns ",
         paste(missed_nms, collapse = ", "), call. = FALSE)
  }

  # should only give one model number in table because no method to select a mod num
  if(!is.null(df$ModelNumber)){
    nmod <- df %>% group_by(responseVariable) %>%
      summarise(nmod = n_distinct(ModelNumber)) %>%
      pull(nmod)

    if(any(nmod > 1)){
      stop("The model coefficient file loaded contains more than one model per response variable",
           call. = FALSE)
    }
  }

  # add expected columns that should never change
  df <- mutate(df, modelVersion = "Johnson",
               ModelNumber = ifelse(responseVariable == "recruitment", "M4", "M1"),
               Type = "National")

  # expected values
  diff_res <- setdiff(unique(df$responseVariable), unique(popGrowthTableJohnsonECCC$responseVariable))

  if(!setequal(unique(popGrowthTableJohnsonECCC$responseVariable), unique(df$responseVariable))){
    stop("The model coefficient file loaded contains unrecognized responseVariable: ",
         paste(diff_res, collapse = ", "), call. = FALSE)
  }

  diff_coef <- setdiff(unique(df$Coefficient),  c("Intercept", "Precision",
                                                  "Anthro", "fire_excl_anthro"))

  if(length(diff_coef) > 0){
    stop("The model coefficient file loaded contains unrecognized Coefficient: ",
         paste(diff_coef, collapse = ", "), call. = FALSE)
  }

  testStdCI <- df %>% mutate(stdOrCI = !is.na(StdErr)|(!is.na(lowerCI) & !is.na(upperCI)))

  if(!all(testStdCI$stdOrCI)){
    stop("The model coefficient file loaded is missing StdErr or lowerCI and upperCI for:\n",
         "femaleSurvival: ",
         testStdCI %>% filter(!stdOrCI, responseVariable == "femaleSurvival") %>%
           pull(Coefficient) %>% paste0(collapse = ", "), "\n",
         "recruitment: ",
         testStdCI %>% filter(!stdOrCI, responseVariable == "recruitment") %>%
           pull(Coefficient) %>% paste0(collapse = ", "),
         call. = FALSE)
  }

  return(df)
}

# Plots -------------------------------------------------------------------

plotRes <- function(allRes, parameter,obs=NULL,lowBound=0,highBound=1,simRange=NULL,facetVars=NULL){
  #allRes=scResults$ksDists; parameter="Recruitment";obs=scResults$obs.all;lowBound=0; highBound=1;simRange=scResults$sim.all;facetVars=c("P","sQ")

  if(is.null(facetVars)){
    titleFontSize = 16
    labFontSize=14
    breakInterval=1
  }else{
    titleFontSize = 11
    labFontSize=10
    breakInterval=2
  }
  if(is.element("KSDistance",names(allRes))){
    #plot Kolmogorov Smirnov distances
    KS=T
    allRes$Mean=allRes$KSDistance
  }else{
    KS=F
  }

  df <- subset(allRes, allRes$Parameter == parameter)

  if(nrow(df) < 1){
    stop()
  }

  if(!is.null(obs)){
    pr = parameter
    obs=subset(obs,parameter==pr)
  }

  if(!KS&!is.null(simRange)){
    pr = parameter
    simRange=subset(simRange,parameter==pr)

    df$Type="local";simRange$Type="national"
    nameSel<-c(c("Year","Mean","Lower 95% CRI", "Upper 95% CRI","Type"),facetVars)
    df=rbind(subset(df,select=nameSel),subset(simRange,select=nameSel))
    df$grp = df$Type
    if(!is.null(facetVars)){
      for(i in facetVars)
        df$grp=paste0(df$grp,df[[i]])
    }
    x1 <- ggplot(df, aes(x = Year, y=Mean,fill=Type,col=Type))

  }else{
    x1 <- ggplot(df, aes(x = Year, y=Mean))

  }
  x2 <- x1 + theme_classic() + xlab("Year")+ ylab(parameter) +
    geom_line(aes(x = Year, y=Mean), size = 1.75) +
    theme(axis.text.y = element_text(size=labFontSize),
          axis.text.x = element_text(angle = 90, vjust = 0.5, size=labFontSize),
          axis.title.x=element_text(size=titleFontSize,face="bold"),
          axis.title.y=element_text(size=titleFontSize,face="bold")) +
    scale_x_continuous(breaks=seq(min(df$Year, na.rm = TRUE),
                                  max(df$Year, na.rm = TRUE),breakInterval))

  if(!KS){
    x2 <- x2+    geom_ribbon(aes(ymin = `Lower 95% CRI`, ymax = `Upper 95% CRI`),
                             show.legend = FALSE, alpha = 0.25,colour=NA)+
      scale_y_continuous(limits=c(ifelse(any(df$`Lower 95% CRI` < lowBound), NA, lowBound),
                                  ifelse(any(df$`Upper 95% CRI` > 1), NA, highBound)))
  }

  if(!KS&!is.null(obs)){
    obs$Type="local"
    obs$obsError=F
    obs$obsError[obs$type=="observed"]=T
    x2<-x2+geom_point(data=obs,aes(x=Year,y=Mean, shape=obsError),col="black",show.legend=T)+
      scale_shape_manual(values = c(16,2))
  }

  if(!is.null(facetVars)){
    if(length(facetVars)==2){
      x2<-x2+facet_grid(as.formula(paste(facetVars[1],"~", facetVars[2])),labeller = "label_both")
    }
  }
  if(!KS & (parameter=="Population growth rate")){
    x2<-x2+geom_hline(yintercept=1, color = "black")
  }

  x2
}

simulateObservations<-function(cs,printPlot=F,cowCounts="cowCounts.csv",
                               freqStartsByYear="freqStartsByYear.csv",
                               collarNumYears=4,collarOffTime=5,
                               collarOnTime=8,distScen = NULL,
                               popGrowthTable = popGrowthTableJohnsonECCC,
                               survivalModelNumber = "M1",
                               recruitmentModelNumber = "M4"){
  #printPlot=T;cowCounts=ePars$cowCounts;freqStartsByYear=ePars$freqStartsByYear;collarNumYears=ePars$collarNumYears;collarOffTime=ePars$collarOffTime;collarOnTime=ePars$collarOnTime
  #distScen = NULL;popGrowthTable = popGrowthTableJohnsonECCC;survivalModelNumber = "M1";recruitmentModelNumber = "M4"
  if(is.character(cowCounts)){
    cowCounts =read.csv2(paste0("tabs/",cowCounts))
  }
  if(is.character(freqStartsByYear)){
    freqStartsByYear=read.csv2(paste0("tabs/",freqStartsByYear))
  }

  #Simulate covariate table
  #TO DO: in UI, option for user to load a table rather than simulate one
  if(is.null(distScen)){
    covariates<-simCovariates(cs$iA,cs$iF,cs$P+cs$J,cs$aS,cs$aSf,cs$P+1)
    simDisturbance=covariates
    simDisturbance$Year = cs$iYr+simDisturbance$time-1

    write.csv(simDisturbance,"tabs/simDisturbance.csv") #note default file for UI is always the last scenario run
    write.csv(simDisturbance,paste0("tabs/simDisturbance",cs$label,".csv"))
  } else {
    simDisturbance <- distScen
    simDisturbance$time = simDisturbance$Year - cs$iYr + 1
    simDisturbance <- filter(simDisturbance, Year <= (cs$iYr + cs$P - 1 + cs$J) &
                               Year >= cs$iYr)
  }

  #simulate true population trajectory
  popMetrics<-simTrajectory(numYears=cs$P+cs$J,covariates=simDisturbance,
                            popGrowthTable = popGrowthTable,
                            survivalModelNumber = survivalModelNumber,
                            recruitmentModelNumber = recruitmentModelNumber,
                            recSlopeMultiplier=cs$rS,
                            sefSlopeMultiplier=cs$sS,recQuantile=cs$rQ,sefQuantile=cs$sQ,N0=cs$N0)

  simDisturbance$time=NULL
  if(printPlot){
    #TO DO: save info on true population dynamics, add to projection plots for comparison
    base1 <- ggplot(data = popMetrics, aes(x = Timestep, y = Amount, colour = Replicate,
                                           group = Replicate)) +
      geom_line() +
      facet_wrap(~MetricTypeID, scales = "free") +
      xlab("Time") +
      theme(legend.position = "none")
    print(base1)
  }

  popMetricsWide<-pivot_wider(popMetrics,id_cols=c(Replicate,Timestep),names_from=MetricTypeID,values_from=Amount)
  exData<-subset(popMetricsWide,(Timestep<=cs$P))
  exData$Year = cs$iYr+exData$Timestep-1

  #Now apply observation process model to get simulated calf:cow and survival data.
  #Use sample sizes in example input data e.g. Eaker

  #reduce sim data tables to length of observations prior to max year
  minYr = cs$iYr
  maxYr = cs$iYr+cs$P+cs$J-1

  #given observed total animals & proportion calfs/cows from simulation - get calf/cow ratio
  ageRatioOut<-simCalfCowRatios(cowCounts,minYr,exData)
  ageRatioOut$HerdCode="ALAP" #TO DO: remove option for more than one herd in input files, UI, and all associated code...
  ageRatioOut=subset(ageRatioOut,select=c(names(cowCounts)))
  #TO DO: ensure UI code uses column names rather than column positions, and is not sensitive to rearrangement of columns
  write.csv(ageRatioOut,"tabs/simAgeRatio.csv")
  write.csv(ageRatioOut,paste0("tabs/simAgeRatio",cs$label,".csv"))

  #simulate survival data from survival propability.
  #TO DO: allow user to enter table of starts by year, collarOffTime, and collarNumYears
  #for each animal, construct simulated observations
  simSurvObs<-simSurvivalData(freqStartsByYear,exData,collarNumYears,collarOffTime,collarOnTime)

  #TO DO: ensure UI code uses column names rather than column positions, and is not sensitive to rearrangement of columns
  write.csv(simSurvObs,"tabs/simSurvData.csv")
  write.csv(simSurvObs,paste0("tabs/simSurvData",cs$label,".csv"))
  #TO DO: UI option to easily select among available scenarios.

  return(list(minYr=minYr,maxYr=maxYr,simDisturbance=simDisturbance,simSurvObs=simSurvObs,ageRatioOut=ageRatioOut,exData=exData,cs=cs))
}


runScnSet<-function(scns,ePars,simBig,survAnalysisMethod = "KaplanMeier",getKSDists=T){
  #ePars=eParsIn;survAnalysisMethod="KaplanMeier";getKSDists=T
  scns<-fillDefaults(scns)
  for(p in 1:nrow(scns)){
    #p=1
    cs =scns[p,]

    if(is.element("st",names(cs))){
      ePars$freqStartsByYear$numStarts=cs$st
    }
    if(is.element("cw",names(cs))){
      ePars$cowCounts$Count=cs$cw
    }
    oo = simulateObservations(cs,cowCounts=ePars$cowCounts,freqStartsByYear=ePars$freqStartsByYear,collarNumYears=ePars$collarNumYears,collarOffTime=ePars$collarOffTime,collarOnTime=ePars$collarOnTime)
    betaPriors<-getPriors(cs)
    minYr=min(oo$exData$Year)
    maxYr=max(oo$simDisturbance$Year)
    out<-runRMModel(survData=oo$simSurvObs,ageRatio.herd=oo$ageRatioOut,disturbance=oo$simDisturbance,
                    betaPriors=betaPriors,startYear = minYr,endYear=maxYr,N0=cs$N0,survAnalysisMethod = survAnalysisMethod)
    outTabs<-getOutputTables(result=out$result,startYear=minYr,endYear=maxYr,survInput=out$survInput,oo=oo,simBig=simBig,getKSDists = getKSDists)

    if(p==1){
      rr.summary.all = outTabs$rr.summary.all
      sim.all = outTabs$sim.all
      obs.all=outTabs$obs.all
      ksDists = merge(outTabs$ksDists,cs)

    }else{
      rr.summary.all = rbind(rr.summary.all,outTabs$rr.summary.all)
      sim.all = rbind(sim.all,outTabs$sim.all)
      obs.all=rbind(obs.all,outTabs$obs.all)
      ksDists = rbind(ksDists,merge(outTabs$ksDists,cs))
    }

  }
  return(list(rr.summary.all=rr.summary.all,sim.all=sim.all,obs.all=obs.all,ksDists=ksDists))
}


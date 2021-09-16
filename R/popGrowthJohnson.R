#' Implementation of the Johnson 2020 population model
#'
#' @param N
#' @param numSteps
#' @param Rec_bar
#' @param S_bar
#' @param P_0
#' @param P_K
#' @param alpha
#' @param beta
#' @param Kmultiplier
#' @param r_max
#' @param sexRatio
#' @param interannualVar
#' @param probOption Choices are "binomial","continuous" or "matchJohnson2020".
#' 
#' @export
popGrowthJohnson <- function(N,
                             numSteps,
                             Rec_bar,
                             S_bar,
                             P_0 = 1,
                             P_K = 0.6,
                             alpha = 1,
                             beta = 4,
                             Kmultiplier = 100,
                             r_max = 1.3,
                             sexRatio=0.5,
                             minRec=0,
                             maxRec=0.41,
                             minSadF=0.61,
                             maxSadF=1,
                             interannualVar = list(Rec_CV=0.46,S_CV=0.08696),
                             probOption="binomial"){
  rr=data.frame(N=N)
  Rec_bar[Rec_bar<0]=0
  S_bar[S_bar<0]=0

  if(!is.element("Rec_phi",names(interannualVar))){
    Rec_bar=sexRatio*Rec_bar 
  }else{
    #Phi is precision of calf cow ratio, not recruitment.
    minRec=minRec/sexRatio
    maxRec=maxRec/sexRatio
  }
  if(probOption=="matchJohnson2020"){
    roundDigits=0
    doBinomial=F
  }else{
    if(probOption=="continuous"){
      roundDigits=100
      doBinomial=F
    }else{
      roundDigits=0
      doBinomial=T
    }
    rK <- Kmultiplier * N 
  }
  
  for(t in 1:numSteps){
    print(paste("projecting step ",t))
    if(is.null(interannualVar)||is.na(interannualVar)||((length(interannualVar)==1)&&!interannualVar)){
      Rec_t= Rec_bar
      S_t = S_bar
    }else{
      Rec_t = addInterannualVar(Rec_bar,interannualVar,type="Rec",minV =minRec,maxV=maxRec)  
      S_t = addInterannualVar(S_bar,interannualVar,type="S",minV =minSadF,maxV=maxSadF)      
    }

    if(is.element("Rec_phi",names(interannualVar))){
      Rec_t=sexRatio*Rec_t 
    }  
    
    Ntm1=N
    
    #Note: following SpaDES code from ECCC_CaribouPopnProjections.Rmd
    if(doBinomial){
      n_deaths <- rbinom(length(N),N,(1 - S_t))
    }else{
      n_deaths <- round(N * (1 - S_t),roundDigits)
    }

    surviving_adFemales <- N - n_deaths

    if(probOption=="matchJohnson2020"){
      rK <- Kmultiplier * N 
    }

    n_recruitsUnadjDD <- surviving_adFemales * Rec_t 

    if(Kmultiplier){
      adjDDRtProportion <- (P_0 -
                              ((P_0 - P_K) *
                                 (surviving_adFemales/rK)^beta)) * 
        surviving_adFemales/(surviving_adFemales+1e-6 + alpha)
      
      adjDDRtProportion[adjDDRtProportion<0] <- 0
      adjDDRtProportion[adjDDRtProportion>1] <- 1
    }else{
      adjDDRtProportion=1
    }
    if(doBinomial){
      n_recruits <- rbinom(length(N),surviving_adFemales,Rec_t*adjDDRtProportion)
    }else{
      n_recruits <- round(n_recruitsUnadjDD * adjDDRtProportion,roundDigits) 
    }  
    N <- surviving_adFemales + n_recruits
    
    if(probOption=="matchJohnson2020"){
      ad = adjustN(N,Ntm1,r_max,denominatorAdjust=1e-06,roundDigits=roundDigits)
    }else{
      ad = adjustN(N,Ntm1,r_max,roundDigits=roundDigits)
    }
    if(sum(is.na(ad$N))>0){stop()}
    
    N=ad$N
    rr[paste0("lam",t)]=ad$Lambda    
  }
  
  lamBits = names(rr)[grepl("lam",names(rr))]
  rr$lambda=matrixStats::rowMeans2(as.matrix(subset(rr,select=lamBits)),na.rm=T)
  rr=subset(rr,select=setdiff(names(rr),lamBits))
  rr$N=N 
  return(rr)
}

adjustN<-function(N,Ntm1,r_max=NA,denominatorAdjust=0,roundDigits=50){
  N[N<0]=0
  N[Ntm1==0]=0
  Lambda = N/(Ntm1+ denominatorAdjust)
  
  if(!is.na(r_max)){
    Lambda[Ntm1==0]=0
    N[Lambda>r_max]=round(Ntm1[Lambda>r_max] * r_max,roundDigits)
    Lambda = N/(Ntm1+ denominatorAdjust)
  }
  
  if(denominatorAdjust==0){
    Lambda[Ntm1==0]=0
  }
  return(list(N=N,Lambda=Lambda))
}

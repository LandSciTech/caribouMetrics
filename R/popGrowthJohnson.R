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
#' 

popGrowthJohnson <- function(N,
                             numSteps,
                             Rec_bar,
                             S_bar,
                             P_0 = 0.95,
                             P_K = 0.6,
                             alpha = 1,
                             beta = 4,
                             Kmultiplier = 100,
                             r_max = 1.3,
                             sexRatio=0.5,
                             interannualVar = list(Rec_shape1=1.62,Rec_shape2=3.44,S_shape1=13.98,S_shape2=2.51)){
  #N=pars$N0;numSteps=numSteps;Rec_bar=pars$R_bar;S_bar=pars$S_bar;interannualVar=list(Rec_shape1=1.62,Rec_shape2=3.44,S_shape1=13.98,S_shape2=2.51);P_0=0.95;P_K=0.6;alpha=1;beta=4;Kmultiplier=100;r_max=1.3
  #N=48;numSteps=20;Rec_bar=0.2436859;S_bar=0.8933278;sexRatio=0.5;interannualVar=NA;P_0=0.95;P_K=0.6;alpha=1;beta=4;Kmultiplier=100;r_max=1.3
  
  rr=data.frame(N=N)
  rK <- Kmultiplier * N
  Rec_bar[Rec_bar<0]=0
  S_bar[S_bar<0]=0
  
  Rec_bar=sexRatio*Rec_bar 
  rK <- Kmultiplier * N
  
  for(t in 1:numSteps){
    print(paste("projecting step ",t))
    #t=1

    #TO DO: interannual variability should be in generatePopGrowthPredictions function, not here.
    #Note interannual variability slows things down a lot, and doesn't make much difference unless initial population size is very low.
    #If worrying about small populations sizes need to switch from continuous to stochastic model with discrete number of births/deaths. 
    if(is.null(interannualVar)||is.na(interannualVar)||((length(interannualVar)==1)&!interannualVar)){
      Rec_t= Rec_bar
      S_t = S_bar
    }else{
      Rec_phi = unique(interannualVar$Rec_shape1+interannualVar$Rec_shape2)
      if(length(Rec_phi)>1){
        stop("handle vector of recruitment precision parameters")
      }
      S_phi = unique(interannualVar$S_shape1+interannualVar$S_shape2)
      if(length(S_phi)>1){
        stop("handle vector of survival precision parameters")
      }
      
      Rec_t = betaSample(Rec_bar,Rec_phi)
      S_t = betaSample(S_bar,S_phi)
    }
    
    Ntm1=N
    
    #Note: following SpaDES code from ECCC_CaribouPopnProjections.Rmd. Description in the paper is not adequate for reproducing the model.  
    n_deaths <- (N * (1 - S_t)) #note this is not the right way to do stochastic modelling
    #bernoulli sampling would be the way to go, but slower to vectorize
    #Sutherland implementation uses rounding, which in turn yields wierd results for low population sizes with no interannual variability
    
    surviving_adFemales <- N - n_deaths
    
    n_recruitsUnadjDD <- (surviving_adFemales * Rec_t) # projected rec (after DD effects):cow at yr end
    
    adjDDRtProportion <- (P_0 -
                            ((P_0 - P_K) *
                               (surviving_adFemales/rK)^beta)) * 
      (surviving_adFemales/((surviving_adFemales+1e-6) + alpha))
    
    adjDDRtProportion[adjDDRtProportion<0] <- 0
    adjDDRtProportion[adjDDRtProportion>1] <- 1
    n_recruits <- (n_recruitsUnadjDD * adjDDRtProportion) # projected rec (after DD effects):cow at yr end
    
    N <- surviving_adFemales + n_recruits;
    
    ad = adjustN(N,Ntm1,r_max)
    N=ad$N
    rr[paste0("lam",t)]=ad$Lambda    
  }
  
  lamBits = names(rr)[grepl("lam",names(rr))]
  rr$lambda=matrixStats::rowMeans2(as.matrix(subset(rr,select=lamBits)))
  rr=subset(rr,select=setdiff(names(rr),lamBits))
  rr$N=N #report final population size
  return(rr)
}

#' Implementation of the ECCC 2011 population model
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
#' 

popGrowthECCC2011 <- function(N,
                             numSteps,
                             Rec_bar,
                             S_bar,
                             Kmultiplier = 20,
                             r_max = 1.3,
                             sexRatio=0.5,
                             interannualVar = list(Rec_shape1=1.62,Rec_shape2=3.44,S_shape1=13.98,S_shape2=2.51)){
  #N=pars$N0;numSteps=numSteps;Rec_bar=pars$R_bar;S_bar=pars$S_bar;interannualVar=NA;P_0=0.95;P_K=0.6;alpha=1;beta=4;Kmultiplier=20;r_max=1.3
  #list(Rec_shape1=1.62,Rec_shape2=3.44,S_shape1=13.98,S_shape2=2.51)
  #N=48;numSteps=20;Rec_bar=0;S_bar=0.8933278;sexRatio=0.5;interannualVar=NA;P_0=0.95;P_K=0.6;alpha=1;beta=4;Kmultiplier=100;r_max=1.3
  
  rr=data.frame(N=N)
  rK <- Kmultiplier * N
  Rec_bar[Rec_bar<0]=0
  S_bar[S_bar<0]=0
  
  Rec_bar=sexRatio*Rec_bar #only 50% of calves are females.???
  
  for(t in 1:numSteps){
    print(paste("projecting step ",t))
    #t=1
    
    #TO DO: interannual variability should be in generatePopGrowthPredictions function, not here.
    #Note interannual variability slows things down a lot, and doesn't make much difference unless initial population size is very low.
    #If worrying about small populations sizes need to switch from continuous to stochastic model with discrete number of births/deaths. 
    if(is.null(interannualVar)||is.na(interannualVar)||((length(interannualVar)==1)&!interannualVar)){
      Rec_t= Rec_bar
      S_t = S_bar
    }else{
      if(is.element("R_var",names(interannualVar))){
        #interannual variability as in ECCC2011
        Rec_t=rnorm(length(Rec_bar),Rec_bar,interannualVar$R_var)
        S_t=rnorm(length(S_bar),S_bar,interannualVar$S_var)
      }else{
        
        Rec_phi = unique(interannualVar$Rec_shape1+interannualVar$Rec_shape2)
        if(length(Rec_phi)>1){
          stop("handle vector of recruitment precision parameters")
        }
        S_phi = unique(interannualVar$S_shape1+interannualVar$S_shape2)
        if(length(S_phi)>1){
          stop("handle vector of survival precision parameters")
        }
        Rec_t = betaSample(Rec_bar,Rec_phi)
        S_t = betaSample(S_bar,S_phi)
      }
    }
    Ntm1=N
    
    N = N -N*(1-S_t)+N*S_t*Rec_t
    
    N[N>rK]=N[N>rK]+r_max*N[N>rK]*(1-N[N>rK]/rK[N>rK])
    
    ad = adjustN(N,Ntm1)
    N=ad$N
    rr[paste0("lam",t)]=ad$Lambda    
  }
  
  lamBits = names(rr)[grepl("lam",names(rr))]
  rr$lambda=matrixStats::rowMeans2(as.matrix(subset(rr,select=lamBits)),na.rm=T)
  rr=subset(rr,select=setdiff(names(rr),lamBits))
  rr$N=N #report final population size
  return(rr)
}

adjustN<-function(N,Ntm1,r_max=NA){
  N[N<0]=0
  N[Ntm1==0]=0
  Lambda = N/Ntm1
  
  if(!is.na(r_max)){
    Lambda[Ntm1==0]=0
    N[Lambda>r_max]=(Ntm1[Lambda>r_max] * r_max)
    Lambda = N/Ntm1
  }
  Lambda[Ntm1==0]=NA
  return(list(N=N,Lambda=Lambda))
}
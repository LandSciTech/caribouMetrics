#' Take a sample from a beta distribution
#'
#' @param x
#' @param phi
#' @param useQuantiles
#' 
betaSample<-function(x,phi,useQuantiles=F){
  #x=predictedTableSD[1,]
  bShapes = simstudy::betaGetShapes(x,phi) 
  
  if((length(useQuantiles)==1)&&!useQuantiles){
    return(rbeta(length(x),bShapes$shape1,bShapes$shape2))
  }else{
    if(length(useQuantiles)!=length(x)){
      q=getQuantiles(x)
    }else{
      q=useQuantiles
    }
    
    qq = qbeta(q,bShapes$shape1,bShapes$shape2)
    return(qq)
  }
}

#' Take a sample from a normal distribution
#'
#' @param x
#' @param sd
#' @param useQuantiles
#' 

normalSample<-function(x,sd,useQuantiles=F){
  #x=predictedTableSD[1,]
  
  if((length(useQuantiles)==1)&&!useQuantiles){
    return(rnorm(length(x),sd))
  }else{
    if(length(useQuantiles)!=length(x)){
      q=getQuantiles(x)
    }else{
      q=useQuantiles
    }
    qq = qnorm(q,mean=x,sd=sd)
    return(qq)
  }
}

#' Take a sample from a lognormal distribution
#'
#' @param x
#' @param sd
#' @param useQuantiles
#' 

lnormSample<-function(x,sd,useQuantiles=F){
  #x=predictedTableSD[1,]
  
  if((length(useQuantiles)==1)&&!useQuantiles){
    return(rnorm(length(x),sd))
  }else{
    if(length(useQuantiles)!=length(x)){
      q=getQuantiles(x)
    }else{
      q=useQuantiles
    }
    
    qq = qlnorm(q,sd)
    return(qq)
  }
}

#' Take a sample from a beta distribution
#'
#' @param x
#' @param phi
#' @param useQuantiles
#' 
betaSample<-function(x,phi,useQuantiles=F){
  #x=predictedTableSD[1,]
  bShapes = simstudy::betaGetShapes(x,phi) 
  
  if((length(useQuantiles)==1)&&!useQuantiles){
    return(rbeta(length(x),bShapes$shape1,bShapes$shape2))
  }else{
    if(length(useQuantiles)!=length(x)){
      q=getQuantiles(x)
    }else{
      q=useQuantiles
    }
    
    qq = qbeta(q,bShapes$shape1,bShapes$shape2)
    return(qq)
  }
}


estBetaParams <- function(mu, sigma){
  
  if(any(mu<0)){
    print("ERROR (estBetaParam): mu must not be less than 0. Returning NULL.")
    return()
  } 
  
  if(any(mu>1)){
    print("ERROR (estBetaParam): mu must not be greater than 1. Returning NULL.")
    return()
  } 
  
  if(any(sigma<0)){
    print("ERROR (estBetaParam): sigma must not be negative. Returning NULL.")
    return()
  } 
  
  alpha <- ((1-mu)/sigma^2 - 1/mu) * mu^2
  beta <- alpha * (1/mu - 1)
  return(list(alpha=alpha, beta=beta))
  
}
fillNAsWithMean <- function(vector) {
  vector[is.na(vector)] <- mean(vector, na.rm = TRUE)
  vector
}

addInterannualVar<-function(bar,interannualVar,type,minV,maxV){
  #bar=pars$expectedRec;type="Rec";minV=pars$minRec; maxV=pars$maxRec; interannualVar = list(Rec_CV=pars$procEARCV,S_CV=pars$procESadFCV)
  #bar=pars$expectedSadF;type="S";minV=pars$minSadF; maxV=pars$maxSadF; interannualVar = list(Rec_CV=pars$procEARCV,S_CV=pars$procESadFCV)
  
  if(is.element(paste0(type,"_CV"),names(interannualVar))){
    #reproducing ECCC_CaribouPopnProjection - see line 143 etc of functions.R
    ProcVar <- (bar * interannualVar[[paste0(type,"_CV")]])^2
    ProcVar=fillNAsWithMean(ProcVar)
    BetaPars  <- estBetaParams(bar, ProcVar)
    BetaPars$alpha[BetaPars$alpha < 0] <- 0.01
    BetaPars$beta[BetaPars$beta < 0] <- 0.01
    #median(BetaPars$alpha+BetaPars$beta)
    interannualVar[[paste0(type,"_alpha")]]=BetaPars$alpha
    interannualVar[[paste0(type,"_beta")]]=BetaPars$beta
  }
  if(is.element(paste0(type,"_phi"),names(interannualVar))){
    bShapes = simstudy::betaGetShapes(bar,interannualVar[[paste0(type,"_phi")]]) 
    interannualVar[[paste0(type,"_alpha")]]=bShapes$shape1
    interannualVar[[paste0(type,"_beta")]]=bShapes$shape2
  }
  
  bar_t = truncdist::rtrunc(length(bar), 
                 spec="beta", 
                 shape1=interannualVar[[paste0(type,"_alpha")]], 
                 shape2= interannualVar[[paste0(type,"_beta")]],
                 a=minV,
                 b=maxV)
  
  return(bar_t)
}

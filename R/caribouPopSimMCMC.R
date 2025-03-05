#' Get a set of simulation results from fitted demographic models in raw form
#'  
#' Assumes that rec_pred and surv_pred each include the same years and populations.TO DO: check this.
#' @param popInfo 
#' @param rec_pred mcmcarray returned by predict_calf_cow function of bboutools R package. 
#' @param surv_pred mcmcarray returned by predict_survival function of bboutools R package.
#' @param initYear numeric. Initial year.
#' @inheritParams caribouPopGrowth
#'
#' @return a data frame with results from [caribouPopGrowth()] for each set of survival/recruitment predictions.
#' 
#' @family demography
#' @export
#'
caribouPopSimMCMC <- function(popInfo, rec_pred, surv_pred, initYear=NULL,...) {
  #assumes rec_pred and surv_pred are mcmcarrays returned by predict_survival and predict_calf_cow
  #functions of bboutools R package, and that each includes the same years and populations.
  #TO DO: checks to ensure these conditions are met.
  #initYear = 2024
  data_sur = surv_pred$data
  data_rec = rec_pred$data
  sur=collapse_chains(surv_pred$samples)
  rec=collapse_chains(rec_pred$samples)
  years = sort(unique(surv_pred$data$Annual))

  if(is.null(initYear)){initYear = as.numeric(as.character(years))}
  years = years[as.numeric(as.character(years))>=initYear]

  if(length(popInfo)>1){
    if(!is.element("id",names(popInfo))){popInfo=merge(popInfo,data.frame(id=1:dim(rec)[2]))}
    popInfo = popInfo[order(popInfo$id,popInfo$pop_name),]
    for(nn in setdiff(names(popInfo),c("pop_name","id"))){
      txt  = paste0(nn," = popInfo[['",nn,"']]")
      eval(parse(text=txt))
    }
  }else{
    N0=popInfo
  }
  for (ts in 1:length(years)) {
    #ts=1
    #print(ts)
    yr = years[ts]
    S_samp = sur[,,data_sur$Annual %in% yr]
    R_samp = rec[,,data_rec$Annual %in% yr]

    if(min(dim(R_samp)[2],dim(S_samp)[2])==0){break}

    if(class(S_samp)[1]=="numeric"){S_samp = matrix(S_samp)}
    if(class(R_samp)[1]=="numeric"){R_samp = matrix(S_samp)}

    rownames(S_samp)=seq(1,nrow(S_samp))
    colnames(S_samp)=levels(data_sur$PopulationID)
    S_samp_long = matrix(S_samp,dimnames=list(t(outer(colnames(S_samp), rownames(S_samp), FUN=paste)), NULL))

    rownames(R_samp)=seq(1:nrow(R_samp))
    colnames(R_samp)=levels(data_rec$PopulationID)
    R_samp_long = matrix(R_samp,dimnames=list(t(outer(colnames(R_samp), rownames(R_samp), FUN=paste)), NULL))

    if(length(N0)==1){N0=rep(N0,length(S_samp_long))}
    if (ts == 1) {
      out <- caribouPopGrowth(N0=N0,
                              numSteps = 1,
                              interannualVar = F,
                              R_bar = R_samp_long, S_bar = S_samp_long, l_S = 0, h_R = 1,...
      )
      out$lab <- matrix(rownames(R_samp_long),nrow(out),1)
      out$year <- yr
      out$time <- ts
      outBit <- out
    } else {
      outBit <- caribouPopGrowth(outBit$N,
                                 numSteps = 1, interannualVar = F,
                                 R_bar = R_samp_long, S_bar = S_samp_long, l_S = 0, h_R = 1,...)
      outBit$lab <- matrix(rownames(R_samp_long),nrow(outBit),1)
      outBit$year <- yr
      outBit$time <- ts
      out <- rbind(out, outBit)
    }
  }
  out[c("PopulationName", "id")] <- do.call(rbind, strsplit(out$lab, " "))
  out$id=as.numeric(out$id)
  return(out)
}

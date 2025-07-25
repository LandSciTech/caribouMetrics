#' Get a set of simulation results from fitted demographic models in raw form
#'  
#' Assumes that rec_pred and surv_pred each include the same years and populations.TO DO: check this.
#' @param popInfo If NA (default) predictions are made without populations size, density dependence, or demographic stochasticity. See [caribouPopGrowth()] for details.
#' @param rec_pred results returned by [bb_fit_recruitment()] or [bb_predict_calf_cow_ratio()] functions of bboutools R package, or recruit_fit returned by [bbouMakeSummaryTable()].
#' @param surv_pred bboufit object return by [bb_fit_survival()] or [bb_predict_survival()] functions of bboutools R package, or surv_fit returned by [bbouMakeSummaryTable()].
#' @param initYear numeric. Initial year.
#' @param correlateRates logical. Set TRUE to force correlation between recruitment and survival. Ignored 
#' @param returnExpected logical. Default FALSE. Set TRUE to return expected values of R, S, and lambda (without interannual variation). Ignored if rec_pred/surv_pred are [bb_predict_calf_cow_ratio()]/[bb_predict_survival()] results.
#' @inheritParams caribouPopGrowth
#'
#' @return a data frame with results from [caribouPopGrowth()] for each set of survival/recruitment predictions.
#' 
#' @family demography
#' @export
#'
caribouPopSimMCMC <- function(popInfo=NA, rec_pred, surv_pred, initYear=NULL,correlateRates=FALSE,returnExpected=FALSE,c=formals(caribouPopGrowth)$c,...) {
  #TO DO: checks to ensure assumptions about form of rec_pred and surv_pred are correct

  inyears<-unique(c(rec_pred$data$Year,surv_pred$data$Year))
  if(!is.null(initYear)){
    inyears<-inyears[inyears>=initYear]
  }
  
  if(returnExpected&&(is.element("bboufit",class(rec_pred))!=is.element("bboufit",class(surv_pred)))){
    stop("To get expected values rec_pred and surv_pred should be [bb_fit_recruitment()]/[bb_fit_survival()] results")
  }
  if(is.element("bboufit",class(rec_pred))){
    if(returnExpected){
      rec_pred <- bboutools::bb_predict_calf_cow_ratio(rec_pred,year=F,conf_level=F)
    }else{
      rec_pred <- bboutools::bb_predict_calf_cow_ratio(rec_pred,year=T,conf_level=F)
    }
  }else{returnExpected=F}
  
  if(is.element("bboufit",class(surv_pred))){
    if(returnExpected){
      surv_pred <- bboutools::bb_predict_survival(surv_pred,year=F,conf_level=F)
    }else{
      surv_pred <- bboutools::bb_predict_survival(surv_pred,year=T,conf_level=F)
    }
  }
  
  if(is.element(NA,rec_pred$data$Annual)){
    rec_pred$data$Annual=0
  }
  if(is.element(NA,surv_pred$data$Annual)){
    surv_pred$data$Annual=0
  }
  
  surv_pred$data$Annual = as.numeric(as.character(surv_pred$data$Annual))
  rec_pred$data$Annual = as.numeric(as.character(rec_pred$data$Annual))  
  data_sur = surv_pred$data
  data_rec = rec_pred$data
  
  S_lookup = unique(subset(data_sur,select=c(Annual,PopulationID)))
  R_lookup =  subset(data_rec,select=c(Annual,PopulationID))
  
  if(class(surv_pred$samples)=="mcmc.list"){
    nns <- colnames(surv_pred$samples[[1]])
    nni <- seq(1,length(nns))
    if(returnExpected){nni <- nni[grepl("bar",nns,fixed=T)]}else{nni <- nni[!grepl("bar",nns,fixed=T)]}
    for(i in 1:length(surv_pred$samples)){
      surv_pred$samples[[i]]<-surv_pred$samples[[i]][,nni]
    }
  }
  if(class(rec_pred$samples)=="mcmc.list"){
    nns <- colnames(rec_pred$samples[[1]])
    nni <- seq(1,length(nns))
    if(returnExpected){nni <- nni[grepl("bar",nns,fixed=T)]}else{nni <- nni[!grepl("bar",nns,fixed=T)]}
    for(i in 1:length(rec_pred$samples)){
      rec_pred$samples[[i]]<-rec_pred$samples[[i]][,nni]
    }
  }
  
  if(class(surv_pred$samples)=="mcmcarray"){
    sur = mcmcr::collapse_chains(surv_pred$samples)
  }else{
    sur2 <- as.matrix(data.table::rbindlist(lapply(surv_pred$samples, as.data.frame)))
    sur <- array(0,dim=c(1,dim(sur2)))
    sur[1,1:nrow(sur2),1:ncol(sur2)]<-sur2[1:nrow(sur2),1:ncol(sur2)]
    S_lookup = S_lookup[order(S_lookup$PopulationID,S_lookup$Annual),]
  }
  
  if(class(rec_pred$samples)=="mcmcarray"){
    rec = mcmcr::collapse_chains(rec_pred$samples)
  }else{
    rec2 <- as.matrix(data.table::rbindlist(lapply(rec_pred$samples, as.data.frame)))
    rec <- array(0,dim=c(1,dim(rec2)))
    rec[1,1:nrow(rec2),1:ncol(rec2)]<-rec2[1:nrow(rec2),1:ncol(rec2)]
    R_lookup = R_lookup[order(R_lookup$PopulationID,R_lookup$Annual),]
  }

  years = sort(unique(rec_pred$data$Annual))
  
  #force correlation between mean recruitment and mean survival
  if(correlateRates){
    recMeans = matrix(0,nrow=dim(rec)[2],ncol=length(unique(R_lookup$PopulationID)))
    surMeans = matrix(0,nrow=dim(sur)[2],ncol=length(unique(S_lookup$PopulationID)))
    for(j in 1:ncol(recMeans)){
      if(length(years)>1){
        recMeans[,j] = recMeans[,j] + rowSums(rec[,,as.numeric(R_lookup$PopulationID) %in% j])
        surMeans[,j] = surMeans[,j] + rowSums(sur[,,as.numeric(S_lookup$PopulationID) %in% j])
      }else{
        recMeans[,j] = recMeans[,j] + rec[,,as.numeric(R_lookup$PopulationID) %in% j]
        surMeans[,j] = surMeans[,j] + sur[,,as.numeric(S_lookup$PopulationID) %in% j]
      }
    }
  }
  
  if(is.null(initYear)){initYear = as.numeric(as.character(years))}
  years = years[as.numeric(as.character(years))>=initYear]
  

  if(length(popInfo)>1){
    popInfo=merge(popInfo,data.frame(id=1:dim(rec)[2]))
    popInfo = popInfo[order(popInfo$id,popInfo$pop_name),]
    for(nn in setdiff(names(popInfo),c("pop_name","id"))){
      txt  = paste0(nn," = popInfo[['",nn,"']]")
      eval(parse(text=txt))
    }
  }else{
    N0=popInfo
  }
  
  if(is.null(c)){c<-1}
  first=T
  
  for (ts in 1:length(years)) {
    #ts=1
    yr <- years[ts]
    
    S_samp <- sur[,,S_lookup$Annual %in% yr]
    R_samp <- rec[,,R_lookup$Annual %in% yr]
    
    if(class(S_samp)=="numeric"){S_samp<-as.matrix(S_samp,ncol=1)}
    if(class(R_samp)=="numeric"){R_samp<-as.matrix(R_samp,ncol=1)}
    
    if((ncol(S_samp)==0)|(ncol(R_samp)==0)){next}
    
    if(correlateRates){
      for(j in 1:ncol(recMeans)){
        R_samp[,j] =  R_samp[,j][order(recMeans[,j])]   
        S_samp[,j] =  S_samp[,j][order(surMeans[,j])]   
      }
    }
    
    if(is.null(dim(S_samp))|is.null(dim(R_samp))){
      if(!is.null(dim(R_samp))&is.null(dim(S_samp))){stop("Handle this case")}
      minDim <- min(length(S_samp),length(R_samp))
      if(minDim==0){next}
      S_samp_long <- S_samp[1:minDim]
      R_samp_long <- R_samp[1:minDim]
      #lab <- paste(levels(data_sur$PopulationID),seq(1:length(S_samp_long)))
      labs <- paste(levels(data_rec$PopulationID),seq(1:length(R_samp_long)))
    }else{
      minDim = min(dim(R_samp)[1],dim(S_samp)[1])
      if(minDim==0){next}

      S_samp <- as.matrix(S_samp[1:minDim,],ncol=ncol(S_samp))
      rownames(S_samp)=seq(1,nrow(S_samp))
      colnames(S_samp)=levels(data_sur$PopulationID)
      S_samp_long = matrix(S_samp,dimnames=list(t(outer(colnames(S_samp), rownames(S_samp), FUN=paste)), NULL))

      R_samp <- as.matrix(R_samp[1:minDim,],ncol=ncol(R_samp))
      rownames(R_samp)=seq(1:nrow(R_samp))
      colnames(R_samp)=levels(data_rec$PopulationID)
      R_samp_long = matrix(R_samp,dimnames=list(t(outer(colnames(R_samp), rownames(R_samp), FUN=paste)), NULL))
      
      labs <- matrix(rownames(R_samp_long),ncol=1)
    }
    
    if(length(N0)==1){N0=rep(N0,length(S_samp_long))}
    

    if (first) {
      out <- caribouPopGrowth(N0=N0,
                              numSteps = 1,
                              interannualVar = F,
                              R_bar = R_samp_long, S_bar = S_samp_long,c=c, l_S = 0, h_R = 1,K=FALSE,...
      )
      out$lab <- labs 
      out$Year <- yr
      out$time <- ts
      outBit <- out
      first <- F
    } else {
      outBit <- caribouPopGrowth(outBit$N,
                                 numSteps = 1, interannualVar = F,
                                 R_bar = R_samp_long, S_bar = S_samp_long,c=c, l_S = 0, h_R = 1,K=FALSE,...)
      
      outBit$lab <- labs
      outBit$Year <- yr
      outBit$time <- ts
      out <- rbind(out, outBit)
    }
  }
  
  if(min(out$Year)==0){
    out$Year=NULL
    out=merge(out,data.frame(Year=inyears))
  }
  
  for (cc in names(out)){
    #cc = "lambda"
    if(is.element("matrix",class(out[[cc]]))){
      if(dim(out[[cc]])[2]>1){
        stop("Unexpected dimensions from caribouPopSimMCMC. Handle this error.")
      }
      out[[cc]] <- out[[cc]][,1]
    }
  }
  out[c("PopulationName", "id")] <- do.call(rbind, strsplit(out$lab, " "))
  out$id=as.numeric(out$id)
  out$time=NULL
  
  if(returnExpected){
    changeNames <- names(out)
    changeNames[!is.element(changeNames,c("Year","lab","id","PopulationName"))] <-
      paste0(changeNames[!is.element(changeNames,c("Year","lab","id","PopulationName"))],"_bar")
    changeNames <- gsub("_t","",changeNames,fixed=T)
    names(out) <- changeNames  
  }    
  return(out)
}

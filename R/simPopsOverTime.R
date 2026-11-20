#' Simulate caribou population over time
#'
#' If `dynamicRates = FALSE` then `R_samp` and `S_samp` are constant rates over
#' time and their length is the number of populations. If `dynamicRates = TRUE`
#' and `R_samp` and `S_samp` are vectors they represent the rate at each time
#' step and there length should be equal to `numSteps`. If `dynamicRates = TRUE`
#' and `R_samp` and `S_samp` are matrices then rows represent populations
#' and columns represent timesteps, so cell `[i,j]` of the matrix is the rate
#' for population i at timestep j
#' 
#' `numSteps` is the number of timesteps in the simulation while `stepLength` is
#' the `numSteps` in each call to `caribouPopGrowth`
#'
#' @noRd
simPopsOverTime <- function(N0, numSteps, R_samp, S_samp, interannualVar, dynamicRates = FALSE, stepLength = 1, ...) {
  if(dynamicRates){
    if(!is.null(dim(R_samp))){
      stopifnot(ncol(R_samp) == numSteps)
      onePop <- FALSE
    }else {
      stopifnot(length(R_samp) == numSteps)
      onePop <- TRUE
    }
  }
  for (ts in 1:numSteps) {
    if(dynamicRates){
      if(onePop){
        R_use <- data.frame(value = R_samp[ts]) %>% setNames(ts)
        S_use <- data.frame(value = S_samp[ts]) %>% setNames(ts)
      } else {
        R_use <- R_samp[,ts, drop = FALSE]
        S_use <- S_samp[,ts, drop = FALSE]
      } 
    } else {
      R_use <- data.frame(value = R_samp)
      S_use <- data.frame(value = S_samp)
    }
    R_use <- na.omit(R_use)
    S_use <- na.omit(S_use)
    if (ts == 1) {
      
      if(length(unique(N0)) == 1){
        N0 <- rep(unique(N0), nrow(R_use))
      } else if(inherits(N0, "data.frame")){
        N0 <- addN0Variation(N0)$N0
      }
      if (length(N0) != nrow(R_use)) {
        stop("The length of N0 and R_samp do not match.",
             " To add N0 variation specify a data frame that addN0Variation() can use.")
        
        #N0 <- seq(from = N0[1], to = N0[2], by  = 1) %>% round() %>% 
        #  sample(size = nrow(R_use), replace = TRUE)
      }
      
      out <- caribouPopGrowth(N0,
                              numSteps = stepLength,
                              interannualVar = interannualVar,
                              R_bar = R_use[,1], S_bar = S_use[,1],warn=F, ...
      )
      
      if(is.null(rownames(R_use))){
        out$id <- seq(1, nrow(out))
      } else {
        out$id <- rownames(R_use)
      }
      
      out$time <- ifelse(!is.null(colnames(R_samp)), colnames(R_use), ts)
      
      outBit <- out
    } else {
      if(length(outBit$N) != nrow(R_use)){
        if(length(unique(N0)) > 1){
          stop("Range of N0 only supported for static rates when there are NAs in R_samp")
        }
        # Add rows missing in previous rounds and set N to N0
        outBit <- left_join(R_use %>% tibble::rownames_to_column(),
                            outBit %>% select(-N0) %>% mutate(id = as.character(id)),
                            by = join_by(rowname == id)) %>% 
          mutate(N = ifelse(is.na(N), unique(N0), N)) %>% 
          select(any_of(colnames(outBit))) 
      }
      outBit <- caribouPopGrowth(outBit$N,
                                 numSteps = stepLength, interannualVar = interannualVar,
                                 R_bar = R_use[,1], S_bar = S_use[,1],warn=F, ...
      )
      if(is.null(rownames(R_use))){
        outBit$id <- seq(1, nrow(outBit))
      } else {
        outBit$id <- rownames(R_use)
      }
      
      outBit$time <- ifelse(!is.null(colnames(R_samp)), colnames(R_use), ts)
      
      out <- rbind(out, outBit)
    }
  }
  return(out)
}


prepareTrajectories <- function(trajs, returnSamples){
  #trajs = mod_samps
  #get the lambda percentile for each id - to allow users to select extreme examples
  simSum <- trajs  %>%
    group_by(id) %>%
    summarize(MeanLam = mean(lambda,na.rm=T))
  simSum <-simSum[order(simSum$MeanLam),]
  simSum$lamPercentile <- round(100*seq(1:nrow(simSum))/nrow(simSum))
  simSum$MeanLam=NULL
  trajs <- merge(trajs,simSum)
  
  trajs <- convertTrajectories(trajs)
  
  simSum <- summarizeTrajectories(trajs, returnSamples = returnSamples)
  return(simSum)
}

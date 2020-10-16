#' Calculate the moving window average of raster
#'
#' Calculates a moving window average for each raster cell using a circular
#' window. If \code{veloxpkg = TRUE} it uses the velox package.
#'
#' @param rast a Raster* object 
#' @param radius the radius of the circular window
#' @param nms the names of each raster layer
#' @param veloxpkg should the velox package be used? 

movingWindowAvg <- function(rast, radius, nms, veloxpkg = FALSE, offset = TRUE, 
                            na.rm = TRUE, pad = FALSE, padValue = NA){
  
  if(veloxpkg){
    if(!requireNamespace("velox", quietly = TRUE)){
      message("velox package not available, using raster package instead")
      out <-  movingWindowAvg(rast, radius, nms, veloxpkg = FALSE)
      return(out)
    }
    nl <- nlayers(rast)
    cf2 <- focalWeight(rast, radius, "circle")
    # cf2[which(cf2 > 0)] <- 1
    rast <- velox::velox(rast)
    rast$sumFocal(cf2, bands = 1:nl)
    
    out <- rast$as.RasterBrick() %>% `names<-`(nms)
    
    # out <- raster::as.data.frame(out) %>% 
    #   mutate(PID = 1:n()) %>% 
    #   set_names(nm = c(nms, "PID"))
    return(out)
  } else {

    nl <- nlayers(rast)
    cf2 <- focalWeight(rast, radius, "circle")
    
    if(offset > 0){
      # add rows to make only far edge points overlap for 4 horizontal/vertical
      # offsets
      nToAdd2 <- nrow(cf2)-1
      cf2_off <- matrix(0, nrow = nrow(cf2) + nToAdd2, ncol = nrow(cf2) + nToAdd2)
      
      cf2_off5 <- cf2_off
      cf2_off5[(nToAdd2/2+1):(nToAdd2/2+nrow(cf2)), (nToAdd2+1):(nToAdd2+nrow(cf2))] <- cf2
      cf2_off6 <- cf2_off
      cf2_off6[(nToAdd2+1):(nToAdd2+nrow(cf2)), (nToAdd2/2+1):(nToAdd2/2+nrow(cf2))] <- cf2
      cf2_off7 <- cf2_off
      cf2_off7[0:(nToAdd2+1), (nToAdd2/2+1):(nToAdd2/2+nrow(cf2))] <- cf2
      cf2_off8 <- cf2_off
      cf2_off8[(nToAdd2/2+1):(nToAdd2/2+nrow(cf2)), 0:(nToAdd2+1)] <- cf2
      
      # The centre "offset"
      cf2_off9 <- cf2_off
      cf2_off9[(nToAdd2/2+1):(nToAdd2/2+nrow(cf2)), (nToAdd2/2+1):(nToAdd2/2+nrow(cf2))] <- cf2
      
      # find top left edge of circle
      test <- FALSE
      for (i in 1:nrow(cf2)) {
        x <- which(cf2[i, ] > 0)
        y <- which(cf2[, i] > 0)
        test <- any(x == i & x == y)
        if(test){
          break
        }
      }
      
      # make cell at i,i overlap point for 4 diagonal offsets
      nToAdd <- nToAdd2-i-2
      
      cf2_off1 <- cf2_off
      cf2_off1[(nToAdd+1):(nToAdd+nrow(cf2)),(nToAdd+1):(nToAdd+nrow(cf2))] <- cf2 
      cf2_off2 <- cf2_off
      cf2_off2[(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1),(nToAdd+1):(nToAdd+nrow(cf2))] <- cf2
      cf2_off3 <- cf2_off
      cf2_off3[(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1),(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1)] <- cf2
      cf2_off4 <- cf2_off
      cf2_off4[(nToAdd+1):(nToAdd+nrow(cf2)),(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1)] <- cf2
      
      if(offset == 16){
        # add 8 more diagonal offsets
        
        nToAdd3 <- nToAdd2 - (i*2 - 1)
        nToAdd4 <- nToAdd + (nToAdd - nToAdd3)/2 
        
        cf2_off10 <- cf2_off
        cf2_off10[(nToAdd3+1):(nToAdd3+nrow(cf2)),(nToAdd4+1):(nToAdd4+nrow(cf2))] <- cf2 
        cf2_off11 <- cf2_off
        cf2_off11[(nToAdd4+1):(nToAdd4+nrow(cf2)),(nToAdd3+1):(nToAdd3+nrow(cf2))] <- cf2 
        
        cf2_off14 <- cf2_off
        cf2_off14[(nrow(cf2)-nToAdd3):(nrow(cf2)*2-nToAdd3-1),(nToAdd4+1):(nToAdd4+nrow(cf2))] <- cf2 
        cf2_off15 <- cf2_off
        cf2_off15[(nrow(cf2)-nToAdd4):(nrow(cf2)*2-nToAdd4-1),(nToAdd3+1):(nToAdd3+nrow(cf2))] <- cf2 
        
        cf2_off12 <- cf2_off
        cf2_off12[(nrow(cf2)-nToAdd4):(nrow(cf2)*2-nToAdd4-1),(nrow(cf2)-nToAdd3):(nrow(cf2)*2-nToAdd3-1)] <- cf2
        cf2_off13 <- cf2_off
        cf2_off13[(nToAdd3+1):(nToAdd3+nrow(cf2)),(nrow(cf2)-nToAdd4):(nrow(cf2)*2-nToAdd4-1)] <- cf2
        
        cf2_off16 <- cf2_off
        cf2_off16[(nrow(cf2)-nToAdd3):(nrow(cf2)*2-nToAdd3-1),(nrow(cf2)-nToAdd4):(nrow(cf2)*2-nToAdd4-1)] <- cf2
        cf2_off17 <- cf2_off
        cf2_off17[(nToAdd4+1):(nToAdd4+nrow(cf2)),(nrow(cf2)-nToAdd3):(nrow(cf2)*2-nToAdd3-1)] <- cf2
        
        cf2 <- cf2_off1 + cf2_off2 + cf2_off3 + cf2_off4 + cf2_off5 + cf2_off6 + 
          cf2_off7 + cf2_off8 + cf2_off10 + cf2_off11+ cf2_off12 + cf2_off13 +
          cf2_off14 + cf2_off15 + cf2_off16 + cf2_off17
        cf2 <- cf2/16
      } else {
        cf2 <- cf2_off1 + cf2_off2 + cf2_off3 + cf2_off4 + cf2_off5 + cf2_off6 + 
          cf2_off7 + cf2_off8 + cf2_off9
        cf2 <- cf2/9
      }
      
    }
    
    
    if(nl == 1){
      rast <- raster::focal(rast, w = cf2, na.rm = na.rm, pad = pad,
                            padValue = padValue)
    } else {
      for(i in 1:nl){
        rast[[i]] <- raster::focal(rast[[i]], w = cf2, na.rm = na.rm, pad = pad,
                                   padValue = padValue)
      }
    }
    
    rast <- rast %>% `names<-`(nms)
    
    return(rast)
  }

   
}
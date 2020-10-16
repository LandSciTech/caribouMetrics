mat <- matrix(c(rep(1,5000), rep(2, 5000)), nrow = 100, ncol = 100)
exps <- raster(mat, xmn = 0, xmx = 40000, ymn = 0, ymx = 40000) 
lyrs <- layerize(exps)
pt <- st_sfc(st_point(c(20000, 20000))) %>% st_sf(PID = 1) %>%
  set_names(c("PID", "geometry")) %>% st_set_geometry("geometry")

# Internal of function helpful for interactive testing
cf2 <- focalWeight(lyrs, 5600, "circle")

# add rows to make only far edge points overlap
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

# make cell at i,i centre add rows and cols all 4 sides
nToAdd <- nToAdd2-i-2

cf2_off1 <- cf2_off
cf2_off1[(nToAdd+1):(nToAdd+nrow(cf2)),(nToAdd+1):(nToAdd+nrow(cf2))] <- cf2 
cf2_off2 <- cf2_off
cf2_off2[(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1),(nToAdd+1):(nToAdd+nrow(cf2))] <- cf2
cf2_off3 <- cf2_off
cf2_off3[(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1),(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1)] <- cf2
cf2_off4 <- cf2_off
cf2_off4[(nToAdd+1):(nToAdd+nrow(cf2)),(nrow(cf2)-nToAdd):(nrow(cf2)*2-nToAdd-1)] <- cf2

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


cfGrid1 <- cf2_off1 %>% as.data.frame() %>% mutate(X = 1:nrow(cf2_off1)) %>% 
  gather(key, value, -X) %>% 
  mutate(key = stringr::str_remove(key, "V") %>% as.numeric()) %>% 
  filter(value > 0) %>% 
  select(-value)

cfGrid10 <- cf2_off10 %>% as.data.frame() %>% mutate(X = 1:nrow(cf2_off1)) %>% 
  gather(key, value, -X) %>% 
  mutate(key = stringr::str_remove(key, "V") %>% as.numeric()) %>% 
  filter(value > 0) %>% 
  select(-value)

cfGrid9 <- cf2_off9 %>% as.data.frame() %>% mutate(X = 1:nrow(cf2_off9)) %>% 
  gather(key, value, -X) %>% 
  mutate(key = stringr::str_remove(key, "V") %>% as.numeric()) %>% 
  filter(value > 0) %>% 
  select(-value)

plot(cfGrid9 %>% rbind(c(40,40)))
points(cfGrid1, col = "red")
points(data.frame(X = 29, key = 29), col = "yellow")

cf2rasoff1 <- raster(cf2_off1)
cf2rasoff2 <- raster(cf2_off2)
cf2rasoff3 <- raster(cf2_off3)
cf2rasoff4 <- raster(cf2_off4)
cf2rasoff5 <- raster(cf2_off5)
cf2rasoff6 <- raster(cf2_off6)
cf2rasoff7 <- raster(cf2_off7)
cf2rasoff8 <- raster(cf2_off8)
cf2rasoff9 <- raster(cf2_off9)
cf2rasoff10 <- raster(cf2_off10)
cf2rasoff11 <- raster(cf2_off11)
cf2rasoff12 <- raster(cf2_off12)
cf2rasoff13 <- raster(cf2_off13)
cf2rasoff14 <- raster(cf2_off14)
cf2rasoff15 <- raster(cf2_off15)
cf2rasoff16 <- raster(cf2_off16)
cf2rasoff17 <- raster(cf2_off17)


cf2rasoffAll16 <- cf2rasoff4+cf2rasoff1+cf2rasoff2+cf2rasoff3+cf2rasoff5+cf2rasoff6+
  cf2rasoff7+cf2rasoff8+cf2rasoff10+cf2rasoff11+cf2rasoff12+cf2rasoff13+cf2rasoff14+
  cf2rasoff15+cf2rasoff16+cf2rasoff17
cf2rasoffAll16 <- cf2rasoffAll16 / 16

cf2rasoffAll9 <- cf2rasoff4+cf2rasoff1+cf2rasoff2+cf2rasoff3+cf2rasoff5+cf2rasoff6+
  cf2rasoff7+cf2rasoff8+cf2rasoff9
cf2rasoffAll9 <- cf2rasoffAll9 / 9

plot(cf2rasoffAll9)
plot(cf2rasoffAll16)

cf2ras <- raster(cf2, xmn = 15800, xmx = 24200, ymn = 15800, ymx = 24200)

# Try with gaussian window
cfG <- focalWeight(lyrs, c(3960, 11200), "Gauss")

cfGras <- raster(cfG)

plot(1:nrow(cf2_off), cf2rasoffAll[ceiling(nrow(cf2_off)/2),])
lines(1:nrow(cfG), cfGras[ceiling(nrow(cfG)/2),])

plot(cfGras)

# how to generalize to any winRad?
# winRad 5600 then sig: 3960 and size: 11200, actual radius: 5800
# winRad 4000 then sig: 2850 and size: 8000, actual radius: 4200
# size seems to be (nrow(cf2)-1)*res(rast)[1]
# sig

nl <- nlayers(rast)

if(nl == 1){
  rast <- raster::focal(rast, w = cf2)
} else {
  for(i in 1:nl){
    rast[[i]] <- raster::focal(rast[[i]], w = cf2)
  }
}

rast <- rast %>% `names<-`(nms)

res <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"))
expect_equal(raster::extract(res$X1, pt), 13/21)


# Use cluster to pick farthest apart cells
cf2 <- focalWeight(lyrs, 5600, "circle")

cfGrid <- cf2 %>% as.data.frame() %>% mutate(X = 1:nrow(cf2)) %>% 
  gather(key, value, -X) %>% 
  mutate(key = stringr::str_remove(key, "V") %>% as.numeric()) %>% 
  filter(value > 0) %>% 
  select(-value)

f_new <- function(dat, n) {
  dmat <- as.matrix(stats::dist(dat))
  r <- sample.int(nrow(dmat), n)
  #r <- which(dat[[1]] == 15)[seq(from = 1, by = floor(max(dat[[1]])/n),
    #                             length.out = n)]
  # pick outer points and center
  nrw <- max(dat[[1]])
  mid <- (nrw/2) %>% ceiling()
  r5 <- c(which(dat[[1]] == mid & dat[[2]] == mid), 
          which(dat[[1]] == 1 & dat[[2]] == mid),
          which(dat[[1]] == nrw & dat[[2]] == mid),
          which(dat[[1]] == mid & dat[[2]] == 1),
          which(dat[[1]] == mid & dat[[2]] == nrw))
  r <- c(r5, sample.int(nrow(dmat), n-5))
  repeat {
    r_old <- r
    for (i in 1:n) {
      # drop the one in the group with the min dist to another in group
      dgroup <- dmat[r, r, drop = FALSE]
      mindgroup <- min(dgroup[lower.tri(dgroup)])
      toDrop <- apply(dgroup, 2, function(x) which(x == mindgroup))
      toDrop<- toDrop[vapply(toDrop, function(x) length(x) > 0, 
                             FUN.VALUE = c(TRUE))] %>% 
        names()  %>% as.numeric() %>% 
        .[!. %in% r5] %>% .[ceiling(length(.)/2)]
      
      tD <- which(r == toDrop)
   
      mm <- dmat[r[-tD], -r[-tD], drop = FALSE]
      k <- which.max(mm[(1:ncol(mm) - 1) * nrow(mm) + max.col(t(-mm))])
      r[tD] <- as.numeric(dimnames(mm)[[2]][k])

    }
    if (identical(r_old, r)) {
      plot(dat)
      points(dat[r,], col = "red")
      return(r)
    }
  }
}

f_new(cfGrid, 16)

plot(cfGrid)
points(cluster::pam(cfGrid, k = 16)$medoids, col = "red")
points(cfGrid[f_new(cfGrid, 9),], col = "red")

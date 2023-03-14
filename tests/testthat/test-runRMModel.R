library(caribouMetrics)
library(dplyr)
library(R2jags)
rjags::load.module("glm")

mod <- runRMModel()

#https://raw.githubusercontent.com/petrkeil/Blog/master/2015_05_13_Survival_analysis/survival_analysis.pdf

cancer <- read.table("http://goo.gl/3cnoam", header=TRUE)
summary(cancer)

censored <- cancer$status==0
is.censored <- censored*1
t.to.death <- cancer$death
t.to.death[censored] <- NA
t.to.death
t.cen <- rep(0, times=length(censored))
t.cen[censored] <- cancer$death[censored]
t.cen
cancer.data <- list(t.to.death = t.to.death,
                    t.cen = t.cen,
                    N = nrow(cancer),
                    group = rep(1:4, each=30))

cat("
model
{
# priors
for(j in 1:4)
{
# prior lambda for each group
lambda[j] ~ dgamma(0.001, 0.001)
mu[j] <- 1/lambda[j] # mean time to death
}
# likelihood
for(i in 1:N)
{
is.censored[i] ~ dinterval(t.to.death[i], t.cen[i])
t.to.death[i] ~ dexp(lambda[group[i]])
}
}
", file="survival_cancer.txt")

library(runjags)
library(coda)
cancer.fit <- run.jags(data = cancer.data,
                       model = "survival_cancer.txt",
                       monitor = c("mu"),
                       sample = 1000, burnin = 1000, n.chains = 3)

# Simulation of monitoring scenarios informed by local data

## 0.1 Examples of initial Bayesian models informed by local data and national-demographic disturbance relationships

We simulate monitoring scenarios informed by example data from the
[`bboudata::bbousurv_a`](https://poissonconsulting.github.io/bboudata/reference/bbousurv_a.html)
and
[`bboudata::bbourecruit_a`](https://poissonconsulting.github.io/bboudata/reference/bbourecruit_a.html)
data sets:

- Informative local data: Observations from 2010 to 2016.
- Limited local data: Observations from 2014 to 2016.

To examine and compare predictions for unobserved years (2018 to 2021)
we add missing data.

Our starting points for analysis of monitoring scenarios are Bayesian
Beta models informed by local data and prior knowledge of national
demographic-disturbance relationships. See compare vignette for details.
We examine the influence of prior information from national-demographic
disturbance relationships in low disturbance (5% anthropogenic, 5% fire)
and high disturbance (80% anthropogenic, 10% fire) scenarios.

``` r
library(caribouMetrics)
# use local version on local and installed on GH
if(requireNamespace("devtools", quietly = TRUE)) devtools::load_all()
library(bboudata)
library(dplyr)
library(ggplot2)

figWidth <- 8
figHeight <- 10

#Note - set niters to 100 to run quickly when testing. Set to 1000 for complete results.
niters <- 1000

useSaved <- T #option to skip slow step of fitting bboutools model
surv_data = bboudata::bbousurv_a %>% filter(Year > 2010)
surv_data_add = expand.grid(Year=seq(2017,2025),Month=seq(1:12),PopulationName=unique(surv_data$PopulationName))
surv_data=merge(surv_data,surv_data_add,all.x=T,all.y=T)
surv_data$StartTotal[is.na(surv_data$StartTotal)]=1

recruit_data=bboudata::bbourecruit_a %>% filter(Year > 2010)
recruit_data_add = expand.grid(Year=seq(2017,2025),PopulationName=unique(recruit_data$PopulationName))
recruit_data=merge(recruit_data,recruit_data_add,all.x=T,all.y=T)
recruit_data$Month[is.na(recruit_data$Month)]=3;recruit_data$Day[is.na(recruit_data$Day)]=15

surv_dataNone <- surv_data %>% filter(Year > 2017)
recruit_dataNone <- recruit_data %>% filter(Year > 2017)

surv_dataLimited <- surv_data %>% filter(Year > 2014)
recruit_dataLimited <- recruit_data %>% filter(Year > 2014)

disturbance = data.frame(Year = unique(surv_data$Year), Anthro = 5, Fire_excl_anthro = 5)
disturbanceHigh <- data.frame(Year = disturbance$Year,Anthro = 80,Fire_excl_anthro = 10)

betaLimited <- estimateBayesianRates(surv_dataLimited, recruit_dataLimited, disturbance=disturbance,niters=niters)
betaInformative <- estimateBayesianRates(surv_data, recruit_data, disturbance=disturbance,niters=niters)
betaLimitedHigh <- estimateBayesianRates(surv_dataLimited, recruit_dataLimited, disturbance=disturbanceHigh,niters=niters)

simInformative <- trajectoriesFromBayesian(betaInformative)
simLimited <- trajectoriesFromBayesian(betaLimited)
simLimitedHigh <- trajectoriesFromBayesian(betaLimitedHigh)
```

## 0.2 Simulating additional monitoring of plausible example trajectories from fitted Bayesian models

We begin by specifying monitoring scenario parameters, and selecting an
example trajectory with a relatively low expected growth rate (i.e. in
the 1st percentile of the distribution of mcmc samples from the initial
model). In this example, we specify 5 additional years of monitoring,
with a target of 30 collars per year, and 6 cows per collared cow in the
composition surveys. We then simulate additional years of monitoring,
combine the simulated and observed data, and reanalyze the combined
data. Methods for simulating monitoring of example trajectories are
described in Hughes et al. ([2025](#ref-hughes_integration_2025)).

``` r
scns=list()
scns$lQuantile=0.01 #select example trajectories with low growth rate
correlateRates = T #Force correlation among demographic rates to examine extreme cases
scns$N0 = NA
scns$projYears = 3
scns$curYear = max(simInformative$recruit_data$Year)- scns$projYears
scns$obsYears = scns$curYear - min(simInformative$recruit_data$Year)
scns$cowMult = 6
scns$collarCount = 30

trajectories <- subset(simInformative$samples,LambdaPercentile == round(scns$lQuantile*100))
trajectories <- subset(trajectories,Replicate==sample(unique(trajectories$Replicate),1))

oo <- simulateObservations(getScenarioDefaults(scns), trajectories,
                           surv_data = simInformative$surv_data, recruit_data=simInformative$recruit_data)
informativeMoreMonitoring <- bayesianTrajectoryWorkflow(surv_data = oo$simSurvObs, recruit_data = oo$simRecruitObs, disturbance=disturbance,niters=niters)
out_tbls <- compareTrajectories(informativeMoreMonitoring, simInitial = simInformative)
typeLabsI <- c("More", "Informative")

#Use bayesianScenariosWorkflow to select trajectory, simulate observations, fit new model and summarize results for the case when initial data is limited.
limitedMoreMonitoring = bayesianScenariosWorkflow(scns,simLimited,niters=niters) 
typeLabsL <- c("More", "Limited")

limitedMoreHigh = bayesianScenariosWorkflow(scns,simLimitedHigh,niters=niters) 
#TO DO: observation points in fig 2?
```

``` r
plotSurvivalSeries(oo$simSurvObs)
```

![Observed (2010-2017) and simulated (2018-2022) survival data for an
example trajectory with a relatively low expected growth rate (i.e. in
the 1st percentile of the distribution of mcmc samples from the initial
model).](combine-observed-simulated_files/figure-html/fig-plot1-1.png)

Figure 1: Observed (2010-2017) and simulated (2018-2022) survival data
for an example trajectory with a relatively low expected growth rate
(i.e. in the 1st percentile of the distribution of mcmc samples from the
initial model).

``` r
#TO DO: fix data gap
```

``` r
rec <- plotCompareTrajectories(out_tbls, "Recruitment", typeLabels = typeLabsI)
surv <- plotCompareTrajectories(out_tbls, "Adult female survival", typeLabels = typeLabsI)
lam <- plotCompareTrajectories(out_tbls, "Population growth rate", typeLabels = typeLabsI,
               lowBound = 0, highBound = 1.5)
rec / surv / lam
```

![95% prior predictive intervals from Beta models with and without
additional monitoring (5 years) of an example trajectory with
informative initial data (5 years), low disturbance, and a relatively
low expected growth rate (Figure
1).](combine-observed-simulated_files/figure-html/fig-plot2-1.png)

Figure 2: 95% prior predictive intervals from Beta models with and
without additional monitoring (5 years) of an example trajectory with
informative initial data (5 years), low disturbance, and a relatively
low expected growth rate (Figure 1).

``` r
rec <- plotCompareTrajectories(limitedMoreMonitoring, "Recruitment", typeLabels = typeLabsL)
surv <- plotCompareTrajectories(limitedMoreMonitoring, "Adult female survival", typeLabels = typeLabsL)
lam <- plotCompareTrajectories(limitedMoreMonitoring, "Population growth rate", typeLabels = typeLabsL,
               lowBound = 0, highBound = 1.5)
rec / surv / lam
```

![95% prior predictive intervals from Beta models with and without
additional monitoring (5 years) of an example trajectory with limited
initial data (2 years), low disturbance, and a relatively low expected
growth rate (Figure
1).](combine-observed-simulated_files/figure-html/fig-plot3-1.png)

Figure 3: 95% prior predictive intervals from Beta models with and
without additional monitoring (5 years) of an example trajectory with
limited initial data (2 years), low disturbance, and a relatively low
expected growth rate (Figure 1).

``` r
rec <- plotCompareTrajectories(limitedMoreHigh, "Recruitment", typeLabels = typeLabsL)
surv <- plotCompareTrajectories(limitedMoreHigh, "Adult female survival", typeLabels = typeLabsL)
lam <- plotCompareTrajectories(limitedMoreHigh, "Population growth rate", typeLabels = typeLabsL,
               lowBound = 0, highBound = 1.5)
rec / surv / lam
```

![95% prior predictive intervals from Beta models with and without
additional monitoring (5 years) of an example trajectory with limited
initial data (2 years), high disturbance, and a relatively low expected
growth rate (Figure
1).](combine-observed-simulated_files/figure-html/fig-plot4-1.png)

Figure 4: 95% prior predictive intervals from Beta models with and
without additional monitoring (5 years) of an example trajectory with
limited initial data (2 years), high disturbance, and a relatively low
expected growth rate (Figure 1).

In these examples, an additional 5 years of monitoring of a plausible
low growth rate trajectory with 5 intial years of data does not provide
additional clarity about whether the population is likely to decline or
not in future (Figure [2](#fig:fig-plot2)). In contrast, an additional 5
years of monitoring does substantially reduce uncertainty when there is
less initial data and disturbance is low (Figure [3](#fig:fig-plot3)).
When disturbance is high, additional monitoring does less to reduce
uncertainty because we have clearer prior expectations from the national
model even when local data is limited (Figure [4](#fig:fig-plot4)).

Examining multiple example trajectories clarifies that there are a range
of possible outcomes from additional monitoring.In the case where true
growth rate is plausibly high (99th percentile of the distribution of
expected growth rate) 5 years of additional monitoring clarifies that
the population is likely viable. In the case where true growth rate is
plausibly low (1st percentile of the distribution of expected growth
rate) 5 years of additional monitoring does not clarify whether the
population is stable or declining. Additional collars do help
distinguish the cases from one another, but additional benefits of that
information should be weighed against additional costs. More formal and
thorough analyses of the value of additional monitoring can be done by
examining outcomes from a larger representative sample of plausible
trajectories from the initial model (as in [Hughes et al.
2025](#ref-hughes_integration_2025)).

``` r
#Add rows to the scenario table to examine multiple monitoring scenarios and example trajectories
scnsMulti=scns;scnsMulti$lQuantile=NULL;scnsMulti$collarCount=NULL
scnsMulti = merge(scnsMulti,expand.grid(lQuantile=c(0.01,0.5,0.99),collarCount=c(15,30)))

limitedMoreMulti = bayesianScenariosWorkflow(scnsMulti,simLimited,niters=niters) 

lmmView <- subset(limitedMoreMulti$rr.summary.all,(Parameter=="Expected growth rate")&(Year==2021),select=c(lQuantile,collarCount,Mean,lower,upper,probViable))
base <- ggplot(lmmView,aes(x=collarCount,y=Mean,ymax=upper,ymin=lower,group=lQuantile,colour=probViable))+geom_point()+geom_errorbar()+ylab("Expected growth rate")
```

``` r
base
```

![Outcomes of additional monitoring vary among example trajectories
(points), and with monitoring effort (collarCount). In these examples
there is 5 years additional monitoring, limited initial data (2 years),
and low disturbance. Example trajectories are from the 99th, 50th, and
1st percentiles of the distribution of expected growth rate in the
initial model informed by 2 years of local
data.](combine-observed-simulated_files/figure-html/fig-plot5-1.png)

Figure 5: Outcomes of additional monitoring vary among example
trajectories (points), and with monitoring effort (collarCount). In
these examples there is 5 years additional monitoring, limited initial
data (2 years), and low disturbance. Example trajectories are from the
99th, 50th, and 1st percentiles of the distribution of expected growth
rate in the initial model informed by 2 years of local data.

TO DO: As in Hughes et al we can also explore monitoring scenarios using
simulated trajectories from the national model instead of trajectories
from a fitted Bayesian model.

## References

Hughes, Josie, Sarah Endicott, Anna M. Calvert, and Cheryl A. Johnson.
2025. “Integration of National Demographic-Disturbance Relationships and
Local Data Can Improve Caribou Population Viability Projections and
Inform Monitoring Decisions.” *Ecological Informatics* 87 (July):
103095. <https://doi.org/10.1016/j.ecoinf.2025.103095>.

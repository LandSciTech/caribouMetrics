# Simulation of monitoring scenarios informed by local data

## 1 Examples of initial Bayesian models informed by local data and national-demographic disturbance relationships

We simulate monitoring scenarios informed by example data from the
[`bboudata::bbousurv_a`](https://LandSciTech.github.io/bboudata/reference/bbousurv_a.html)
and
[`bboudata::bbourecruit_a`](https://LandSciTech.github.io/bboudata/reference/bbourecruit_a.html)
data sets:

- Informative local data: Observations from 2010 to 2015.
- Limited local data: Observations from 2014 and 2015.

To work with model predictions for unobserved years (2016 to 2023) we
add missing years to the input data.

Our starting points for analysis of monitoring scenarios are Bayesian
Beta models informed by local data and prior knowledge of national
demographic-disturbance relationships. See compare vignette for details.
We show a how to work with multiple populations and examine the
influence of prior information from national-demographic disturbance
relationships by constructing a data set that includes two populations
with limited data; population A experiences low disturbance (3%
anthropogenic, 5% fire), and population B experiences high disturbance
(80% anthropogenic, 10% fire).

``` r

library(caribouMetrics)
library(bboutools)
# use local version on local and installed on GH
if(requireNamespace("devtools", quietly = TRUE)) devtools::load_all()
if(requireNamespace("devtools", quietly = TRUE)) devtools::load_all(path = "C:/Users/HughesJo/Documents/gitprojects/bboutools/")
library(bboudata)
library(dplyr)
library(ggplot2)
library(patchwork)

figWidth <- 8
figHeight <- 10

#Note - set niters to 100 to run quickly when testing. Set to 1000 for complete results.
niters <- 1000

surv_data <- addMissingYears(bboudata::bbousurv_a %>% getCaribouYear() %>% 
                               filter(between(CaribouYear, 2010, 2015)),
                             seq(2016,2023)) 
recruit_data <- addMissingYears(bboudata::bbourecruit_a %>% getCaribouYear() %>% 
                                  filter(between(CaribouYear, 2010, 2015)),
                                seq(2016,2023)) 
surv_dataLimited <- surv_data %>% filter(CaribouYear > 2013)
recruit_dataLimited <- recruit_data %>% filter(CaribouYear > 2013)

#Construct data set that includes two example populations with the same initial observation data to show workflow for multiple populations
surv_dataLimited$PopulationName = NULL; recruit_dataLimited$PopulationName = NULL
surv_dataLimited <- merge(surv_dataLimited,data.frame(PopulationName=c("A","B")))
recruit_dataLimited <- merge(recruit_dataLimited,data.frame(PopulationName=c("A","B")))
                          
disturbance = data.frame(Year = unique(surv_data$CaribouYear), Anthro = 3, Fire_excl_anthro = 5)
disturbance$PopulationName = "A"
disturbanceHigh <- data.frame(Year = disturbance$Year,Anthro = 80,Fire_excl_anthro = 10)
disturbanceHigh$PopulationName = "B"
#Disturbance differs between population A and B
disturbanceBoth <- rbind(disturbance,disturbanceHigh)

betaLimited <- estimateBayesianRates(surv_dataLimited, recruit_dataLimited, disturbance=disturbanceBoth,niters=niters)

betaInformative <- estimateBayesianRates(surv_data, recruit_data, disturbance=disturbance,niters=niters)

simInformative <- trajectoriesFromBayesian(betaInformative)
simLimited <- trajectoriesFromBayesian(betaLimited)
```

## 2 Simulating additional monitoring of plausible example trajectories from fitted Bayesian models

We begin by specifying monitoring scenario parameters, and selecting an
example trajectory with a relatively low expected growth rate (i.e. in
the 1st percentile of the distribution of mcmc samples from the initial
model). In this example, we specify 5 additional years of monitoring,
with a target of 30 collars per year, and 6 cows per collared cow in the
composition surveys. Default settings specify that deployed collars last
4 years, and lost collars are replaced in April of each year. We then
simulate additional years of monitoring, combine the simulated and
observed data (Figure [1](#fig:fig-plot1)), and reanalyze the combined
data. Methods for simulating monitoring of example trajectories are
described in Hughes et al. ([2025](#ref-hughes_integration_2025)).

In these examples, an additional 5 years of monitoring of a plausible
low growth rate trajectory with 5 initial years of data does not provide
additional clarity about whether the population is likely to decline or
not in future (Figure [2](#fig:fig-plot2)). An additional 5 years of
monitoring is more likely to reduce uncertainty when there is limited
initial data, particularly when disturbance is low (e.g. Population A in
Figure [3](#fig:fig-plot3)), but outcomes of single trajectory analysis
are stochastic. Thorough analysis of the value of additional monitoring
requires examining outcomes from a larger representative sample of
plausible trajectories from the initial model (as in [Hughes et al.
2025](#ref-hughes_integration_2025)).

Thorough analysis of representative samples of plausible trajectories is
time consuming, and beyond the scope of this vignette. Here we show an
example workflow with three example trajectories for each case (Figure
[4](#fig:fig-plot5)). This example clarifies that there are a range of
possible outcomes from additional monitoring. In the case where true
growth rate is plausibly high (99th percentile of the distribution of
expected growth rate) 5 years of additional monitoring clarifies that
this population is likely viable. In the case where true growth rate is
plausibly low (1st percentile of the distribution of expected growth
rate) 5 years of additional monitoring does not clarify whether the
population is stable or declining. Additional collars do help
distinguish the cases from one another, but additional benefits of that
information should be weighed against additional costs.

``` r

scns=list()
scns$lQuantile = 0.01 #select example trajectories with low growth rate
scns$correlateRates = T #Force correlation among demographic rates to examine extreme cases
scns$N0 = NA
scns$projYears = 3
scns$curYear = max(simInformative$recruit_data$Year)- scns$projYears
scns$obsYears = scns$curYear - min(simInformative$recruit_data$Year)
scns$cowMult = 6
scns$collarCount = 30

trajectories <- subset(simInformative$samples,LambdaPercentile == round(scns$lQuantile*100))
trajectories <- subset(trajectories,Replicate==sample(unique(trajectories$Replicate),1))

oo <- simulateObservations(getScenarioDefaults(scns), trajectories,
                           surv_data = simInformative$surv_data %>%
                             filter(!is.na(Mortalities)),
                           recruit_data = simInformative$recruit_data %>% 
                             filter(!is.na(Cows)))

informativeMoreMonitoring <- bayesianTrajectoryWorkflow(
  surv_data = oo$simSurvObs, recruit_data = oo$simRecruitObs,
  disturbance=disturbance, niters=niters
)

out_tbls <- compareTrajectories(informativeMoreMonitoring, simInitial = simInformative)
typeLabsI <- c("More", "Informative")

#Use bayesianScenariosWorkflow to select trajectory, simulate observations, fit new model and summarize results for the case of two population with limited initial data. In this case, we can see both the "true" population trajectory and the observations in result plots.
#TO DO: include additional projection years
limitedMoreMonitoring = bayesianScenariosWorkflow(scns,simLimited,niters=niters) 
typeLabsL <- c("More", "Limited")
```

``` r

plotSurvivalSeries(oo$simSurvObs)
```

![Observed (2011-2015) and simulated (2016-2020) survival data for an
example trajectory with a relatively low expected growth rate (i.e. in
the 1st percentile of the distribution of mcmc samples from the initial
model). X's indicate mortalities, and lines indicate the number of
deployed collars. Each simulated collar lasts for 4 years. In April of
each year, collars lost throughout the year are
replaced.](combine-observed-simulated_files/figure-html/fig-plot1-1.png)

Figure 1: Observed (2011-2015) and simulated (2016-2020) survival data
for an example trajectory with a relatively low expected growth rate
(i.e. in the 1st percentile of the distribution of mcmc samples from the
initial model). X’s indicate mortalities, and lines indicate the number
of deployed collars. Each simulated collar lasts for 4 years. In April
of each year, collars lost throughout the year are replaced.

``` r

rec <- plotCompareTrajectories(out_tbls, "Recruitment", typeLabels = typeLabsI)
surv <- plotCompareTrajectories(out_tbls, "Adult female survival", typeLabels = typeLabsI)
lam <- plotCompareTrajectories(out_tbls, "Expected growth rate", typeLabels = typeLabsI,
               lowBound = 0.7, highBound = 1.2)
rec / surv / lam
```

![95% predictive intervals from Beta models with and without additional
monitoring (5 years) of an example trajectory with informative initial
data (5 years), low disturbance, and a relatively low expected growth
rate (Figure 1). In this example, additional monitoring does not
substantially reduce uncertainty about population growth
rate.](combine-observed-simulated_files/figure-html/fig-plot2-1.png)

Figure 2: 95% predictive intervals from Beta models with and without
additional monitoring (5 years) of an example trajectory with
informative initial data (5 years), low disturbance, and a relatively
low expected growth rate (Figure 1). In this example, additional
monitoring does not substantially reduce uncertainty about population
growth rate.

``` r

rec <- plotCompareTrajectories(limitedMoreMonitoring, "Recruitment", typeLabels = typeLabsL)
surv <- plotCompareTrajectories(limitedMoreMonitoring, "Adult female survival", typeLabels = typeLabsL)
lam <- plotCompareTrajectories(limitedMoreMonitoring, "Expected growth rate", typeLabels = typeLabsL,
               lowBound = 0.7, highBound = 1.2)
rec / surv / lam
```

![95% predictive intervals from Beta models with and without additional
monitoring (5 years) of example trajectories with limited initial data
(2 years) and a relatively low expected growth rate (Figure 1).
Population A experiences low disturbance, and Population B experiences
high
disturbance.](combine-observed-simulated_files/figure-html/fig-plot3-1.png)

Figure 3: 95% predictive intervals from Beta models with and without
additional monitoring (5 years) of example trajectories with limited
initial data (2 years) and a relatively low expected growth rate (Figure
1). Population A experiences low disturbance, and Population B
experiences high disturbance.

``` r

#Add rows to the scenario table to examine multiple monitoring scenarios and example trajectories
scnsMulti=scns;scnsMulti$lQuantile=NULL;scnsMulti$collarCount=NULL
scnsMulti = merge(scnsMulti,expand.grid(lQuantile=c(0.01,0.5,0.99),collarCount=c(15,30)))

limitedMoreMulti = bayesianScenariosWorkflow(scnsMulti,simLimited,niters=niters) 

lmmView <- subset(limitedMoreMulti$rr.summary.all,(Metric=="Expected growth rate")&(Year==2020),select=c(PopulationName, lQuantile,collarCount,Mean,lower,upper,probViable))
base <- ggplot(lmmView,aes(x=PopulationName,y=Mean,ymax=upper,ymin=lower,group=lQuantile,colour=probViable))+geom_point()+geom_errorbar()+facet_wrap(~collarCount,labeller=label_both)+ylab("Expected growth rate")
```

``` r

base
```

![Outcomes of additional monitoring vary among example trajectories
(points), with monitoring effort (collarCount), and among populations (A
and B). In these examples there is 5 years additional monitoring, and
limited initial data (2 years). Population A experiences low
disturbance, and population B experiences high disturbance. Example
trajectories are from the 99th, 50th, and 1st percentiles of the
distribution of expected growth rate in the initial model informed by 2
years of local data. probViable is the posterior probability of
viability (i.e. the proportion of the estimated posterior predictive
distribution of expected population growth rate that is greater than
0.99).](combine-observed-simulated_files/figure-html/fig-plot5-1.png)

Figure 4: Outcomes of additional monitoring vary among example
trajectories (points), with monitoring effort (collarCount), and among
populations (A and B). In these examples there is 5 years additional
monitoring, and limited initial data (2 years). Population A experiences
low disturbance, and population B experiences high disturbance. Example
trajectories are from the 99th, 50th, and 1st percentiles of the
distribution of expected growth rate in the initial model informed by 2
years of local data. probViable is the posterior probability of
viability (i.e. the proportion of the estimated posterior predictive
distribution of expected population growth rate that is greater than
0.99).

## 3 Specifying a more complex staggered monitoring scenario in tables

``` r

freqStartsByYrA <- data.frame(Year = seq(2017,2021), PopulationName = "A", numTarget = c(30,0,30,0,0)) 
freqStartsByYrB <- data.frame(Year = seq(2017,2021), PopulationName = "B", numTarget = c(0,30,0,30,0)) 
freqStartsByYr <- rbind(freqStartsByYrA,freqStartsByYrB)
#Note these are calendar years. Collars deployed in March before start of 2017 caribou year in April

cowCounts <- data.frame(PopulationName = c("A","B"),cowMult=c(3,4))
cowCounts <- merge(cowCounts,data.frame(Year = unique(freqStartsByYr$Year)))
#Note these are caribou years. 2017 survey happens in spring 2018 and so on. 2021 survey won't happen because no collared animals to survey after March 2021.

eParsIn <- eval(formals(caribouMetrics::bayesianScenariosWorkflow)$ePars)
eParsIn$collarOnTime <- 3
eParsIn$collarOffTime <- 2
eParsIn$cowCounts <- cowCounts
eParsIn$freqStartsByYear <- freqStartsByYr
eParsIn$collarNumYears <- 6

scnsB <- data.frame(obsYears = max(freqStartsByYr$Year)-min(simLimited$recruit_data$CaribouYear)+1,
                   startYear = min(simLimited$summary$Year))
scnsB$projYears <- max(simLimited$summary$Year)-scnsB$obsYears-scnsB$startYear
scnsB$N0 <- scns$N0
scnsB$lQuantile <- scns$lQuantile
scnsB$correlateRates <- T

posteriorResult = bayesianScenariosWorkflow(scnsB,simLimited,eParsIn,
                                            niters=niters,returnSamples = T)

surv_data_sim <- subset(posteriorResult$out$result$surv_data,!is.na(StartTotal))
rec_data_sim <- posteriorResult$out$result$recruit_data
```

``` r

plotSurvivalSeries(surv_data_sim)
```

![Observed (2014-2016) and simulated (2017-2021) survival data for an
example trajectory with a relatively low expected growth rate (i.e. in
the 1st percentile of the distribution of mcmc samples from the initial
model). X's indicate mortalities, and lines indicate the number of
deployed collars. Each simulated collar lasts for 6 years. In March of
each year, collars lost throughout the year are replaced. Top up to 30
collars happens in March 2017/2019 for population A, and March 2018/2020
for population B, as
specified.](combine-observed-simulated_files/figure-html/fig-plot6-1.png)

Figure 5: Observed (2014-2016) and simulated (2017-2021) survival data
for an example trajectory with a relatively low expected growth rate
(i.e. in the 1st percentile of the distribution of mcmc samples from the
initial model). X’s indicate mortalities, and lines indicate the number
of deployed collars. Each simulated collar lasts for 6 years. In March
of each year, collars lost throughout the year are replaced. Top up to
30 collars happens in March 2017/2019 for population A, and March
2018/2020 for population B, as specified.

``` r

base <- ggplot(rec_data_sim,aes(x=Year,y=Cows))+geom_point()+facet_wrap(~PopulationName)+
  ylab("Cows in rec survey")
base
```

![Observed (2014-2016) and simulated (2017-2020) cows in recruitment
survey for an example trajectory with a relatively low expected growth
rate (i.e. in the 1st percentile of the distribution of mcmc samples
from the initial model). The simulated recruitment survey occurs in
caribou years 2017 to 2020 but not in caribou year 2021 because there
are no collared animals in spring 2022 at the end of caribou year
2021.](combine-observed-simulated_files/figure-html/fig-plot7-1.png)

Figure 6: Observed (2014-2016) and simulated (2017-2020) cows in
recruitment survey for an example trajectory with a relatively low
expected growth rate (i.e. in the 1st percentile of the distribution of
mcmc samples from the initial model). The simulated recruitment survey
occurs in caribou years 2017 to 2020 but not in caribou year 2021
because there are no collared animals in spring 2022 at the end of
caribou year 2021.

## References

Hughes, Josie, Sarah Endicott, Anna M. Calvert, and Cheryl A. Johnson.
2025. “Integration of National Demographic-Disturbance Relationships and
Local Data Can Improve Caribou Population Viability Projections and
Inform Monitoring Decisions.” *Ecological Informatics* 87 (July):
103095. <https://doi.org/10.1016/j.ecoinf.2025.103095>.

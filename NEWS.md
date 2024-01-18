# caribouMetrics 0.3.2
**Breaking changes**
* Changed function name:
  - caribouBayesianIPM -> caribouBayesianPM

* Added composition bias
* allow simulated trajectory to start before monitoring to allow fair comparisons. 
* Fix some edge case errors
* Fix saving in caribouHabitat using terra

# caribouMetrics 0.3.1
**Breaking changes**
Changed argument names in `getPriors()`:
rInterannualVar -> rSigmaMean
rInterannualVarSE -> rSigmaSD
sInterannualVar -> sSigmaMean
sInterannualVarSE -> sSigmaSD

# caribouMetrics 0.3.0
* Changed to use terra for all spatial processing. Still accepts RasterLayer inputs but outputs SpatRasters.

# caribouMetrics 0.2.0

* Added a `NEWS.md` file to track changes to the package.
* Added functions for Bayesian population model

**Breaking Changes**
* Changed functions names:
    - fillDefaults -> getScenarioDefaults
    - runRMmodel -> caribouBayesianIPM
    - popGrowthJohnson -> caribouPopGrowth
* Changed argument names throughout Bayesian model functions and N in popGrowthJohnson (now caribouPopGrowth) changed to N0
    
  

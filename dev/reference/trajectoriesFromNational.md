# Get a set of simulation results from the national demographic model

Simulate demograhic rates based on the National demographic -
disturbance model If a disturbance scenario containing Years is supplied
trajectories will show growth of a population over time based on the
National demographic - disturbance model

## Usage

``` r
trajectoriesFromNational(
  replicates = 1000,
  N0 = 1000,
  useQuantiles = NULL,
  populationGrowthTable = NULL,
  cPars = subset(getScenarioDefaults(), select = -iAnthro),
  interannualVar = eval(formals(caribouPopGrowth)$interannualVar),
  disturbance = NULL,
  skipSave = FALSE,
  forceUpdate = FALSE,
  doSummary = TRUE,
  returnSamples = "default",
  numSteps = 1
)
```

## Arguments

- replicates:

  integer. Number of replicate populations.

- N0:

  initial population size

- useQuantiles:

  logical or numeric. If it is a numeric vector it must be length 2 and
  give the low and high limits of the quantiles to use. If
  `useQuantiles != FALSE`, each replicate population is assigned to a
  quantile of the distribution of variation around the expected values,
  and remains in that quantile as covariates change. If
  `useQuantiles = TRUE`, replicate populations will be assigned to
  quantiles in the default range of 0.025 and 0.975.

- populationGrowthTable:

  data.frame.[popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md)
  is included in the package and should be used in most cases. A custom
  table of model coefficients and standard errors or confidence
  intervals can be provided but it must match the column names of
  [popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md).
  If the table does not contain the standard error it is calculated from
  the confidence interval.

- cPars:

  optional. Parameters for calculating composition survey bias term.

- interannualVar:

  list or logical. List containing interannual variability parameters.
  These can be either coefficients of variation (R_CV, S_CV), beta
  precision parameters (R_phi, S_phi), or random effects parameters from
  a logistic glmm (R_annual, S_annual). Set to `FALSE` to ignore
  interannual variability.

- disturbance:

  data frame with Anthro, fire_excl_anthro and Year numeric columns.
  Anthro and fire_excl_anthro are vectors of numbers between 0 and 100
  representing the percentage of the landscape covered by anthropogenic
  disturbance buffered by 500 m, and the percentage covered by fire that
  does not overlap anthropogenic disturbance.

- doSummary:

  logical. Default TRUE. If FALSE returns unprocessed outcomes from
  caribouPopGrowth. If TRUE returns summaries and (if returnSamples = T)
  sample trajectories from prepareTrajectories.

- returnSamples:

  logical. If FALSE returns only summaries. If TRUE returns example
  trajectories as well. By default summaries are not returned unless the
  disturbance data provided contains a column named "Year".

- numSteps:

  numeric. Number of steps to run
  [`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md)
  at each disturbance level.

## Value

Output from caribouPopGrowth function.

## See also

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
[`bbouNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/bbouNationalPriors.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateTrajectoriesFromPosterior.md),
[`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/dev/reference/dataFromSheets.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographicProjectionApp.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateNationalRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/dev/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md)

## Examples

``` r
trajectoriesFromNational()
#> Using saved object
#> $summary
#>     MetricTypeID PopulationName AnthroID fire_excl_anthroID       Mean
#> 102         Rbar       National       15                  0 0.28085341
#> 103         Rbar       National       20                  0 0.25636957
#> 104         Rbar       National       19                  0 0.26020270
#> 105         Rbar       National       16                  0 0.27936749
#> 106         Rbar       National       18                  0 0.27075662
#> 107         Rbar       National       17                  0 0.27223632
#> 108         Rbar       National       57                  0 0.14041246
#> 109         Rbar       National        5                  0 0.32836650
#> 110         Rbar       National       56                  0 0.14000992
#> 111         Rbar       National       27                  0 0.22629788
#> 112         Rbar       National       21                  0 0.25034579
#> 113         Rbar       National       32                  0 0.20642098
#> 114         Rbar       National        3                  0 0.34029490
#> 115         Rbar       National       31                  0 0.20952003
#> 116         Rbar       National        2                  0 0.34479999
#> 117         Rbar       National       30                  0 0.21556853
#> 118         Rbar       National        7                  0 0.32335938
#> 119         Rbar       National       58                  0 0.13988129
#> 120         Rbar       National        6                  0 0.32874617
#> 121         Rbar       National        0                  0 0.35866240
#> 122         Rbar       National       28                  0 0.22373509
#> 123         Rbar       National       22                  0 0.24750586
#> 124         Rbar       National       33                  0 0.20863450
#> 125         Rbar       National        4                  0 0.34098310
#> 126         Rbar       National       55                  0 0.14258881
#> 127         Rbar       National       26                  0 0.22853237
#> 128         Rbar       National       37                  0 0.19378088
#> 129         Rbar       National        8                  0 0.31171082
#> 130         Rbar       National       59                  0 0.13284802
#> 131         Rbar       National       70                  0 0.10911689
#> 132         Rbar       National        1                  0 0.35078628
#> 133         Rbar       National       29                  0 0.22287688
#> 134         Rbar       National       23                  0 0.23795164
#> 135         Rbar       National       34                  0 0.19734160
#> 136         Rbar       National       45                  0 0.16902686
#> 137         Rbar       National      100                  0 0.07103899
#> 138         Rbar       National       67                  0 0.11603610
#> 139         Rbar       National       38                  0 0.18782572
#> 140         Rbar       National        9                  0 0.29985152
#> 141         Rbar       National       60                  0 0.13322832
#> 142         Rbar       National       71                  0 0.10793328
#> 143         Rbar       National       42                  0 0.17574488
#> 144         Rbar       National       13                  0 0.28602676
#> 145         Rbar       National       24                  0 0.24053592
#> 146         Rbar       National       35                  0 0.19965526
#> 147         Rbar       National       46                  0 0.16312032
#> 148         Rbar       National       97                  0 0.07158611
#> 149         Rbar       National       68                  0 0.10971025
#> 150         Rbar       National       39                  0 0.19146742
#> 151         Rbar       National       10                  0 0.30472171
#> 152         Rbar       National       61                  0 0.12713210
#> 153         Rbar       National       72                  0 0.10603053
#> 154         Rbar       National       43                  0 0.17716502
#> 155         Rbar       National       14                  0 0.28518688
#> 156         Rbar       National       25                  0 0.23784257
#> 157         Rbar       National       36                  0 0.19401162
#> 158         Rbar       National       47                  0 0.16164348
#> 159         Rbar       National       98                  0 0.07037824
#> 160         Rbar       National       69                  0 0.10997982
#> 161         Rbar       National       40                  0 0.18552120
#> 162         Rbar       National       11                  0 0.29357478
#> 163         Rbar       National       62                  0 0.12876061
#> 164         Rbar       National       73                  0 0.10055450
#> 165         Rbar       National       44                  0 0.16964768
#> 166         Rbar       National       95                  0 0.07441214
#> 167         Rbar       National       66                  0 0.12207423
#> 168         Rbar       National       77                  0 0.09849421
#> 169         Rbar       National       48                  0 0.15646769
#> 170         Rbar       National       99                  0 0.06681503
#> 171         Rbar       National       53                  0 0.15102282
#> 172         Rbar       National       41                  0 0.18099907
#> 173         Rbar       National       12                  0 0.29360293
#> 174         Rbar       National       63                  0 0.12356999
#> 175         Rbar       National       74                  0 0.10452321
#> 176         Rbar       National       85                  0 0.08276527
#> 177         Rbar       National       96                  0 0.07390205
#> 178         Rbar       National       50                  0 0.15591464
#> 179         Rbar       National       78                  0 0.09715106
#> 180         Rbar       National       49                  0 0.16075387
#> 181         Rbar       National       83                  0 0.09011596
#> 182         Rbar       National       54                  0 0.14401095
#> 183         Rbar       National       82                  0 0.09373984
#> 184         Rbar       National       76                  0 0.10035451
#> 185         Rbar       National       64                  0 0.12284928
#> 186         Rbar       National       75                  0 0.10015130
#> 187         Rbar       National       86                  0 0.08314366
#> 188         Rbar       National       80                  0 0.09224989
#> 189         Rbar       National       51                  0 0.15281346
#> 190         Rbar       National       79                  0 0.09604928
#> 191         Rbar       National       81                  0 0.08911149
#> 192         Rbar       National       84                  0 0.09032287
#> 193         Rbar       National       92                  0 0.07847262
#> 194         Rbar       National       89                  0 0.07891795
#> 195         Rbar       National       91                  0 0.07745230
#> 196         Rbar       National       65                  0 0.11670691
#> 197         Rbar       National       90                  0 0.08174093
#> 198         Rbar       National       87                  0 0.08511965
#> 199         Rbar       National       93                  0 0.07766967
#> 200         Rbar       National       52                  0 0.15161212
#> 201         Rbar       National       88                  0 0.08278466
#> 202         Rbar       National       94                  0 0.07393087
#> 203         Sbar       National       16                  0 0.86759709
#> 204         Sbar       National       22                  0 0.86079347
#> 205         Sbar       National       21                  0 0.86208118
#> 206         Sbar       National       17                  0 0.86409482
#> 207         Sbar       National       20                  0 0.86000560
#> 208         Sbar       National       26                  0 0.85820748
#> 209         Sbar       National       19                  0 0.86395815
#> 210         Sbar       National        4                  0 0.87477966
#> 211         Sbar       National        9                  0 0.86953515
#> 212         Sbar       National        7                  0 0.87055409
#> 213         Sbar       National       18                  0 0.86156209
#> 214         Sbar       National        8                  0 0.86755109
#> 215         Sbar       National       23                  0 0.85778345
#> 216         Sbar       National       11                  0 0.86959050
#> 217         Sbar       National        5                  0 0.87432653
#> 218         Sbar       National        6                  0 0.87088300
#> 219         Sbar       National       27                  0 0.85765368
#> 220         Sbar       National       15                  0 0.86585586
#> 221         Sbar       National       49                  0 0.84303699
#> 222         Sbar       National       10                  0 0.87003128
#> 223         Sbar       National       31                  0 0.85540883
#> 224         Sbar       National       59                  0 0.83323005
#> 225         Sbar       National       13                  0 0.86700313
#> 226         Sbar       National       24                  0 0.86103530
#> 227         Sbar       National       12                  0 0.86625989
#> 228         Sbar       National       46                  0 0.84367824
#> 229         Sbar       National       34                  0 0.85390598
#> 230         Sbar       National       28                  0 0.85783720
#> 231         Sbar       National       56                  0 0.83740781
#> 232         Sbar       National       50                  0 0.83869169
#> 233         Sbar       National       61                  0 0.83283341
#> 234         Sbar       National       32                  0 0.85276300
#> 235         Sbar       National       60                  0 0.83219367
#> 236         Sbar       National       14                  0 0.86626947
#> 237         Sbar       National       25                  0 0.85984110
#> 238         Sbar       National       36                  0 0.84855934
#> 239         Sbar       National       47                  0 0.84275932
#> 240         Sbar       National       35                  0 0.85125858
#> 241         Sbar       National       29                  0 0.85525304
#> 242         Sbar       National       57                  0 0.83686054
#> 243         Sbar       National       51                  0 0.84063221
#> 244         Sbar       National       62                  0 0.83403026
#> 245         Sbar       National       33                  0 0.85045515
#> 246         Sbar       National        0                  0 0.87658590
#> 247         Sbar       National       55                  0 0.83721225
#> 248         Sbar       National       66                  0 0.83009757
#> 249         Sbar       National       37                  0 0.85237549
#> 250         Sbar       National       48                  0 0.84197178
#> 251         Sbar       National       99                  0 0.80910655
#> 252         Sbar       National       30                  0 0.85379748
#> 253         Sbar       National       58                  0 0.83614429
#> 254         Sbar       National       52                  0 0.83963310
#> 255         Sbar       National       63                  0 0.83277772
#> 256         Sbar       National       74                  0 0.82464468
#> 257         Sbar       National        1                  0 0.87714243
#> 258         Sbar       National       96                  0 0.81101667
#> 259         Sbar       National       67                  0 0.82747494
#> 260         Sbar       National       38                  0 0.84845942
#> 261         Sbar       National       89                  0 0.81272133
#> 262         Sbar       National      100                  0 0.80815547
#> 263         Sbar       National       71                  0 0.83048276
#> 264         Sbar       National        2                  0 0.87496426
#> 265         Sbar       National       53                  0 0.84033137
#> 266         Sbar       National       64                  0 0.83338241
#> 267         Sbar       National       75                  0 0.82510165
#> 268         Sbar       National       86                  0 0.81542270
#> 269         Sbar       National       97                  0 0.81194226
#> 270         Sbar       National       68                  0 0.82840438
#> 271         Sbar       National       39                  0 0.84801008
#> 272         Sbar       National       90                  0 0.81340314
#> 273         Sbar       National       44                  0 0.84654338
#> 274         Sbar       National       72                  0 0.82654137
#> 275         Sbar       National        3                  0 0.87514044
#> 276         Sbar       National       54                  0 0.83765879
#> 277         Sbar       National       65                  0 0.82974559
#> 278         Sbar       National       76                  0 0.82193925
#> 279         Sbar       National       87                  0 0.81407106
#> 280         Sbar       National       98                  0 0.80937627
#> 281         Sbar       National       69                  0 0.82552595
#> 282         Sbar       National       40                  0 0.84661539
#> 283         Sbar       National       91                  0 0.81410107
#> 284         Sbar       National       45                  0 0.84466563
#> 285         Sbar       National       73                  0 0.82582168
#> 286         Sbar       National       88                  0 0.81555876
#> 287         Sbar       National       95                  0 0.81128975
#> 288         Sbar       National       70                  0 0.82670522
#> 289         Sbar       National       77                  0 0.82552725
#> 290         Sbar       National       92                  0 0.81309650
#> 291         Sbar       National       42                  0 0.84520040
#> 292         Sbar       National       80                  0 0.82012525
#> 293         Sbar       National       41                  0 0.84433265
#> 294         Sbar       National       79                  0 0.82118763
#> 295         Sbar       National       94                  0 0.81167801
#> 296         Sbar       National       78                  0 0.82263194
#> 297         Sbar       National       93                  0 0.81151520
#> 298         Sbar       National       43                  0 0.84731266
#> 299         Sbar       National       81                  0 0.81856923
#> 300         Sbar       National       83                  0 0.82094789
#> 301         Sbar       National       84                  0 0.81927628
#> 302         Sbar       National       82                  0 0.82046917
#> 303         Sbar       National       85                  0 0.81664433
#> 304            X       National       11                  0 0.14675583
#> 305            X       National        7                  0 0.15612034
#> 306            X       National        6                  0 0.16011315
#> 307            X       National        8                  0 0.15017402
#> 308            X       National        4                  0 0.16329048
#> 309            X       National       38                  0 0.09323467
#> 310            X       National        9                  0 0.14720026
#> 311            X       National       12                  0 0.14708154
#> 312            X       National       35                  0 0.09640157
#> 313            X       National       19                  0 0.12817330
#> 314            X       National       13                  0 0.14212860
#> 315            X       National        5                  0 0.15669887
#> 316            X       National       39                  0 0.09575917
#> 317            X       National       10                  0 0.14587415
#> 318            X       National       17                  0 0.13519115
#> 319            X       National       36                  0 0.09398619
#> 320            X       National       20                  0 0.13155270
#> 321            X       National       14                  0 0.14230826
#> 322            X       National        2                  0 0.16934256
#> 323            X       National       40                  0 0.09202076
#> 324            X       National        1                  0 0.17146978
#> 325            X       National       18                  0 0.13334062
#> 326            X       National       50                  0 0.07923068
#> 327            X       National        0                  0 0.17284331
#> 328            X       National       55                  0 0.07020772
#> 329            X       National       16                  0 0.13679657
#> 330            X       National       37                  0 0.09443992
#> 331            X       National       21                  0 0.12035996
#> 332            X       National       15                  0 0.13533651
#> 333            X       National        3                  0 0.15776812
#> 334            X       National       41                  0 0.08883941
#> 335            X       National       52                  0 0.07562680
#> 336            X       National       63                  0 0.06075969
#> 337            X       National       51                  0 0.07468143
#> 338            X       National       45                  0 0.08503763
#> 339            X       National       56                  0 0.07010727
#> 340            X       National       23                  0 0.11738787
#> 341            X       National       34                  0 0.09782050
#> 342            X       National       22                  0 0.12019519
#> 343            X       National       60                  0 0.06667739
#> 344            X       National       48                  0 0.08012446
#> 345            X       National       42                  0 0.09135383
#> 346            X       National       53                  0 0.07452431
#> 347            X       National       64                  0 0.06011803
#> 348            X       National       31                  0 0.10457320
#> 349            X       National       46                  0 0.08429157
#> 350            X       National       57                  0 0.06943796
#> 351            X       National       24                  0 0.11739614
#> 352            X       National       79                  0 0.04814899
#> 353            X       National       90                  0 0.04081561
#> 354            X       National       61                  0 0.06214078
#> 355            X       National       49                  0 0.08122247
#> 356            X       National       43                  0 0.08720641
#> 357            X       National       54                  0 0.07173680
#> 358            X       National       65                  0 0.05834154
#> 359            X       National       32                  0 0.10697172
#> 360            X       National       47                  0 0.07949841
#> 361            X       National       58                  0 0.06914724
#> 362            X       National       25                  0 0.11910843
#> 363            X       National       80                  0 0.04565175
#> 364            X       National       91                  0 0.03848247
#> 365            X       National       62                  0 0.06470555
#> 366            X       National       29                  0 0.11062569
#> 367            X       National       44                  0 0.08565642
#> 368            X       National       95                  0 0.03641676
#> 369            X       National       66                  0 0.05953605
#> 370            X       National       33                  0 0.10306901
#> 371            X       National       88                  0 0.04100914
#> 372            X       National       59                  0 0.06750172
#> 373            X       National       26                  0 0.11343712
#> 374            X       National       81                  0 0.04470095
#> 375            X       National       92                  0 0.03970082
#> 376            X       National       86                  0 0.04155117
#> 377            X       National       30                  0 0.10710600
#> 378            X       National       85                  0 0.04022665
#> 379            X       National       96                  0 0.03716068
#> 380            X       National       67                  0 0.05617474
#> 381            X       National       78                  0 0.04870701
#> 382            X       National       89                  0 0.04073350
#> 383            X       National      100                  0 0.03618068
#> 384            X       National       27                  0 0.11216489
#> 385            X       National       82                  0 0.04598881
#> 386            X       National       93                  0 0.03865604
#> 387            X       National       87                  0 0.04217296
#> 388            X       National       75                  0 0.04848611
#> 389            X       National       69                  0 0.05517534
#> 390            X       National       97                  0 0.03591241
#> 391            X       National       68                  0 0.05350180
#> 392            X       National       83                  0 0.04601953
#> 393            X       National       73                  0 0.04992421
#> 394            X       National       84                  0 0.04527090
#> 395            X       National       28                  0 0.11132766
#> 396            X       National       70                  0 0.05446737
#> 397            X       National       94                  0 0.03617123
#> 398            X       National       71                  0 0.05369838
#> 399            X       National       76                  0 0.05036530
#> 400            X       National       74                  0 0.05111505
#> 401            X       National       98                  0 0.03439781
#> 402            X       National       77                  0 0.04877167
#> 403            X       National       99                  0 0.03355460
#> 404            X       National       72                  0 0.05367693
#> 405            c       National       25                  0 1.00000000
#> 406            c       National        8                  0 1.00000000
#> 407            c       National        6                  0 1.00000000
#> 408            c       National       24                  0 1.00000000
#> 409            c       National        7                  0 1.00000000
#> 410            c       National       35                  0 1.00000000
#> 411            c       National       10                  0 1.00000000
#> 412            c       National        0                  0 1.00000000
#> 413            c       National        9                  0 1.00000000
#> 414            c       National       26                  0 1.00000000
#> 415            c       National       37                  0 1.00000000
#> 416            c       National       48                  0 1.00000000
#> 417            c       National       36                  0 1.00000000
#> 418            c       National       30                  0 1.00000000
#> 419            c       National        1                  0 1.00000000
#> 420            c       National       12                  0 1.00000000
#> 421            c       National       23                  0 1.00000000
#> 422            c       National       11                  0 1.00000000
#> 423            c       National        5                  0 1.00000000
#> 424            c       National       33                  0 1.00000000
#> 425            c       National       27                  0 1.00000000
#> 426            c       National       38                  0 1.00000000
#> 427            c       National       49                  0 1.00000000
#> 428            c       National       20                  0 1.00000000
#> 429            c       National       31                  0 1.00000000
#> 430            c       National        2                  0 1.00000000
#> 431            c       National       13                  0 1.00000000
#> 432            c       National       68                  0 1.00000000
#> 433            c       National       79                  0 1.00000000
#> 434            c       National       46                  0 1.00000000
#> 435            c       National       34                  0 1.00000000
#> 436            c       National       28                  0 1.00000000
#> 437            c       National       39                  0 1.00000000
#> 438            c       National       50                  0 1.00000000
#> 439            c       National       21                  0 1.00000000
#> 440            c       National       32                  0 1.00000000
#> 441            c       National        3                  0 1.00000000
#> 442            c       National       14                  0 1.00000000
#> 443            c       National       69                  0 1.00000000
#> 444            c       National       80                  0 1.00000000
#> 445            c       National       47                  0 1.00000000
#> 446            c       National       18                  0 1.00000000
#> 447            c       National       29                  0 1.00000000
#> 448            c       National       40                  0 1.00000000
#> 449            c       National       51                  0 1.00000000
#> 450            c       National       22                  0 1.00000000
#> 451            c       National       77                  0 1.00000000
#> 452            c       National        4                  0 1.00000000
#> 453            c       National       15                  0 1.00000000
#> 454            c       National       70                  0 1.00000000
#> 455            c       National       81                  0 1.00000000
#> 456            c       National       92                  0 1.00000000
#> 457            c       National       19                  0 1.00000000
#> 458            c       National       74                  0 1.00000000
#> 459            c       National       41                  0 1.00000000
#> 460            c       National       52                  0 1.00000000
#> 461            c       National       63                  0 1.00000000
#> 462            c       National       78                  0 1.00000000
#> 463            c       National       45                  0 1.00000000
#> 464            c       National       16                  0 1.00000000
#> 465            c       National       71                  0 1.00000000
#> 466            c       National       82                  0 1.00000000
#> 467            c       National       93                  0 1.00000000
#> 468            c       National       60                  0 1.00000000
#> 469            c       National       75                  0 1.00000000
#> 470            c       National       42                  0 1.00000000
#> 471            c       National       53                  0 1.00000000
#> 472            c       National       64                  0 1.00000000
#> 473            c       National       58                  0 1.00000000
#> 474            c       National       90                  0 1.00000000
#> 475            c       National       17                  0 1.00000000
#> 476            c       National       72                  0 1.00000000
#> 477            c       National       83                  0 1.00000000
#> 478            c       National       94                  0 1.00000000
#> 479            c       National       61                  0 1.00000000
#> 480            c       National       76                  0 1.00000000
#> 481            c       National       43                  0 1.00000000
#> 482            c       National       54                  0 1.00000000
#> 483            c       National       65                  0 1.00000000
#> 484            c       National       59                  0 1.00000000
#> 485            c       National       91                  0 1.00000000
#> 486            c       National       85                  0 1.00000000
#> 487            c       National       73                  0 1.00000000
#> 488            c       National       84                  0 1.00000000
#> 489            c       National       95                  0 1.00000000
#> 490            c       National       62                  0 1.00000000
#> 491            c       National       56                  0 1.00000000
#> 492            c       National       44                  0 1.00000000
#> 493            c       National       55                  0 1.00000000
#> 494            c       National       66                  0 1.00000000
#> 495            c       National       87                  0 1.00000000
#> 496            c       National       98                  0 1.00000000
#> 497            c       National       86                  0 1.00000000
#> 498            c       National       97                  0 1.00000000
#> 499            c       National       89                  0 1.00000000
#> 500            c       National       96                  0 1.00000000
#> 501            c       National       67                  0 1.00000000
#> 502            c       National       57                  0 1.00000000
#> 503            c       National       99                  0 1.00000000
#> 504            c       National       88                  0 1.00000000
#> 505            c       National      100                  0 1.00000000
#> 506       lambda       National       15                  0 0.98566800
#> 507       lambda       National        5                  0 1.01474700
#> 508       lambda       National       17                  0 0.98471900
#> 509       lambda       National       19                  0 0.97625200
#> 510       lambda       National       30                  0 0.94712600
#> 511       lambda       National       12                  0 0.99198300
#> 512       lambda       National        6                  0 1.01417700
#> 513       lambda       National       57                  0 0.89813100
#> 514       lambda       National       18                  0 0.97642600
#> 515       lambda       National       16                  0 0.98812700
#> 516       lambda       National       10                  0 1.00040400
#> 517       lambda       National       21                  0 0.96950400
#> 518       lambda       National       28                  0 0.95486000
#> 519       lambda       National       20                  0 0.97467900
#> 520       lambda       National       31                  0 0.94355700
#> 521       lambda       National        2                  0 1.02333300
#> 522       lambda       National       13                  0 0.98947000
#> 523       lambda       National        7                  0 1.00987700
#> 524       lambda       National       58                  0 0.89346300
#> 525       lambda       National       29                  0 0.95420900
#> 526       lambda       National        0                  0 1.03114300
#> 527       lambda       National       11                  0 0.99637400
#> 528       lambda       National       22                  0 0.96472500
#> 529       lambda       National       33                  0 0.93667500
#> 530       lambda       National        4                  0 1.02201900
#> 531       lambda       National       32                  0 0.94287100
#> 532       lambda       National        3                  0 1.01707100
#> 533       lambda       National       14                  0 0.98952000
#> 534       lambda       National        8                  0 1.00246400
#> 535       lambda       National       59                  0 0.89265600
#> 536       lambda       National       70                  0 0.87462400
#> 537       lambda       National        1                  0 1.02811500
#> 538       lambda       National       52                  0 0.90400900
#> 539       lambda       National       23                  0 0.95876700
#> 540       lambda       National       34                  0 0.93994000
#> 541       lambda       National       45                  0 0.91480800
#> 542       lambda       National       56                  0 0.89175800
#> 543       lambda       National       27                  0 0.95774900
#> 544       lambda       National       55                  0 0.90079500
#> 545       lambda       National        9                  0 0.99736100
#> 546       lambda       National       60                  0 0.88976800
#> 547       lambda       National       71                  0 0.87816000
#> 548       lambda       National       42                  0 0.92292800
#> 549       lambda       National       53                  0 0.90570900
#> 550       lambda       National       24                  0 0.95622100
#> 551       lambda       National       35                  0 0.93610700
#> 552       lambda       National       46                  0 0.91501900
#> 553       lambda       National       97                  0 0.84335600
#> 554       lambda       National       68                  0 0.87855300
#> 555       lambda       National       39                  0 0.93010700
#> 556       lambda       National       50                  0 0.90638500
#> 557       lambda       National       61                  0 0.88445500
#> 558       lambda       National       72                  0 0.87591700
#> 559       lambda       National       43                  0 0.92086400
#> 560       lambda       National       54                  0 0.90235700
#> 561       lambda       National       25                  0 0.96199000
#> 562       lambda       National       36                  0 0.92488800
#> 563       lambda       National       47                  0 0.91262900
#> 564       lambda       National       98                  0 0.84204600
#> 565       lambda       National       69                  0 0.87521300
#> 566       lambda       National       40                  0 0.92466500
#> 567       lambda       National       51                  0 0.90679100
#> 568       lambda       National       62                  0 0.88953600
#> 569       lambda       National       73                  0 0.87143000
#> 570       lambda       National       44                  0 0.92624700
#> 571       lambda       National       95                  0 0.84259000
#> 572       lambda       National       26                  0 0.95730600
#> 573       lambda       National       37                  0 0.93302400
#> 574       lambda       National       48                  0 0.91908500
#> 575       lambda       National       99                  0 0.83830500
#> 576       lambda       National       93                  0 0.84592300
#> 577       lambda       National       41                  0 0.92410700
#> 578       lambda       National       92                  0 0.84653800
#> 579       lambda       National       63                  0 0.88801400
#> 580       lambda       National       74                  0 0.87164600
#> 581       lambda       National       85                  0 0.85071100
#> 582       lambda       National       96                  0 0.84275100
#> 583       lambda       National       67                  0 0.87689800
#> 584       lambda       National       38                  0 0.93080600
#> 585       lambda       National       49                  0 0.91153000
#> 586       lambda       National      100                  0 0.84045100
#> 587       lambda       National       94                  0 0.84310100
#> 588       lambda       National       82                  0 0.86046300
#> 589       lambda       National       76                  0 0.86792800
#> 590       lambda       National       64                  0 0.88635100
#> 591       lambda       National       75                  0 0.86310100
#> 592       lambda       National       86                  0 0.85020600
#> 593       lambda       National       80                  0 0.85739800
#> 594       lambda       National       91                  0 0.84870800
#> 595       lambda       National       79                  0 0.86305400
#> 596       lambda       National       90                  0 0.85162800
#> 597       lambda       National       84                  0 0.85568800
#> 598       lambda       National       78                  0 0.86561500
#> 599       lambda       National       83                  0 0.86014800
#> 600       lambda       National       77                  0 0.86617700
#> 601       lambda       National       65                  0 0.88407700
#> 602       lambda       National       88                  0 0.85660600
#> 603       lambda       National       87                  0 0.85307300
#> 604       lambda       National       81                  0 0.85609400
#> 605       lambda       National       89                  0 0.84167700
#> 606       lambda       National       66                  0 0.88198100
#> 607   lambda_bar       National       19                  0 0.97638805
#> 608   lambda_bar       National       17                  0 0.98166005
#> 609   lambda_bar       National       18                  0 0.97815844
#> 610   lambda_bar       National        6                  0 1.01420717
#> 611   lambda_bar       National       21                  0 0.97001039
#> 612   lambda_bar       National       11                  0 0.99721662
#> 613   lambda_bar       National       22                  0 0.96732892
#> 614   lambda_bar       National       10                  0 1.00267989
#> 615   lambda_bar       National        8                  0 1.00276600
#> 616   lambda_bar       National       59                  0 0.88857886
#> 617   lambda_bar       National       20                  0 0.97025176
#> 618   lambda_bar       National        1                  0 1.03111468
#> 619   lambda_bar       National       12                  0 0.99340349
#> 620   lambda_bar       National       23                  0 0.95979226
#> 621   lambda_bar       National        7                  0 1.01125840
#> 622   lambda_bar       National        5                  0 1.01790324
#> 623   lambda_bar       National       16                  0 0.98878162
#> 624   lambda_bar       National        4                  0 1.02386740
#> 625   lambda_bar       National       34                  0 0.93815680
#> 626   lambda_bar       National        9                  0 0.99990820
#> 627   lambda_bar       National       60                  0 0.88762334
#> 628   lambda_bar       National       31                  0 0.94510168
#> 629   lambda_bar       National        2                  0 1.02579797
#> 630   lambda_bar       National       13                  0 0.99090890
#> 631   lambda_bar       National       24                  0 0.96474167
#> 632   lambda_bar       National       35                  0 0.93627138
#> 633   lambda_bar       National       46                  0 0.91249068
#> 634   lambda_bar       National       57                  0 0.89560663
#> 635   lambda_bar       National       28                  0 0.95368240
#> 636   lambda_bar       National       39                  0 0.92914400
#> 637   lambda_bar       National       50                  0 0.90413837
#> 638   lambda_bar       National       61                  0 0.88576258
#> 639   lambda_bar       National       32                  0 0.94076122
#> 640   lambda_bar       National        3                  0 1.02398526
#> 641   lambda_bar       National       14                  0 0.98982721
#> 642   lambda_bar       National       25                  0 0.96200023
#> 643   lambda_bar       National       36                  0 0.93089138
#> 644   lambda_bar       National       47                  0 0.91096495
#> 645   lambda_bar       National       58                  0 0.89452870
#> 646   lambda_bar       National       29                  0 0.95053748
#> 647   lambda_bar       National       40                  0 0.92522539
#> 648   lambda_bar       National       51                  0 0.90493234
#> 649   lambda_bar       National       62                  0 0.88773513
#> 650   lambda_bar       National       33                  0 0.93923658
#> 651   lambda_bar       National        0                  0 1.03374786
#> 652   lambda_bar       National       15                  0 0.98750935
#> 653   lambda_bar       National       26                  0 0.95633402
#> 654   lambda_bar       National       37                  0 0.93505281
#> 655   lambda_bar       National       48                  0 0.90786086
#> 656   lambda_bar       National       99                  0 0.83613168
#> 657   lambda_bar       National       30                  0 0.94576586
#> 658   lambda_bar       National       41                  0 0.92064306
#> 659   lambda_bar       National       52                  0 0.90325891
#> 660   lambda_bar       National       63                  0 0.88425756
#> 661   lambda_bar       National       74                  0 0.86774165
#> 662   lambda_bar       National       45                  0 0.91606078
#> 663   lambda_bar       National       56                  0 0.89597117
#> 664   lambda_bar       National       27                  0 0.95454472
#> 665   lambda_bar       National       38                  0 0.92817678
#> 666   lambda_bar       National       49                  0 0.91077642
#> 667   lambda_bar       National      100                  0 0.83685754
#> 668   lambda_bar       National       71                  0 0.87543852
#> 669   lambda_bar       National       42                  0 0.91949563
#> 670   lambda_bar       National       53                  0 0.90377832
#> 671   lambda_bar       National       64                  0 0.88461343
#> 672   lambda_bar       National       75                  0 0.86643478
#> 673   lambda_bar       National       86                  0 0.84931139
#> 674   lambda_bar       National       97                  0 0.84101399
#> 675   lambda_bar       National       68                  0 0.87387271
#> 676   lambda_bar       National       79                  0 0.86062389
#> 677   lambda_bar       National       90                  0 0.84655801
#> 678   lambda_bar       National       44                  0 0.91827664
#> 679   lambda_bar       National       72                  0 0.87035717
#> 680   lambda_bar       National       43                  0 0.92235276
#> 681   lambda_bar       National       54                  0 0.89805290
#> 682   lambda_bar       National       65                  0 0.87814663
#> 683   lambda_bar       National       76                  0 0.86319663
#> 684   lambda_bar       National       87                  0 0.84862073
#> 685   lambda_bar       National       98                  0 0.83787154
#> 686   lambda_bar       National       69                  0 0.87106068
#> 687   lambda_bar       National       80                  0 0.85798807
#> 688   lambda_bar       National       91                  0 0.84565990
#> 689   lambda_bar       National       85                  0 0.85045406
#> 690   lambda_bar       National       73                  0 0.86744330
#> 691   lambda_bar       National       67                  0 0.87551963
#> 692   lambda_bar       National       55                  0 0.89680322
#> 693   lambda_bar       National       66                  0 0.88075072
#> 694   lambda_bar       National       77                  0 0.86617401
#> 695   lambda_bar       National       88                  0 0.84934771
#> 696   lambda_bar       National       82                  0 0.85895483
#> 697   lambda_bar       National       70                  0 0.87175651
#> 698   lambda_bar       National       81                  0 0.85505256
#> 699   lambda_bar       National       92                  0 0.84498527
#> 700   lambda_bar       National       94                  0 0.84169836
#> 701   lambda_bar       National       78                  0 0.86262527
#> 702   lambda_bar       National       89                  0 0.84475431
#> 703   lambda_bar       National       96                  0 0.84100726
#> 704   lambda_bar       National       93                  0 0.84307067
#> 705   lambda_bar       National       83                  0 0.85793628
#> 706   lambda_bar       National       84                  0 0.85620468
#> 707   lambda_bar       National       95                  0 0.84150433
#> 708  recruitment       National        8                  0 0.30034803
#> 709  recruitment       National       19                  0 0.25634661
#> 710  recruitment       National        6                  0 0.32022630
#> 711  recruitment       National        4                  0 0.32658096
#> 712  recruitment       National       20                  0 0.26310540
#> 713  recruitment       National       35                  0 0.19280313
#> 714  recruitment       National        7                  0 0.31224068
#> 715  recruitment       National        9                  0 0.29440052
#> 716  recruitment       National       24                  0 0.23479229
#> 717  recruitment       National       39                  0 0.19151835
#> 718  recruitment       National        2                  0 0.33868511
#> 719  recruitment       National       13                  0 0.28425721
#> 720  recruitment       National        5                  0 0.31339773
#> 721  recruitment       National       16                  0 0.27359314
#> 722  recruitment       National       10                  0 0.29174830
#> 723  recruitment       National       21                  0 0.24071992
#> 724  recruitment       National       36                  0 0.18797238
#> 725  recruitment       National        3                  0 0.31553624
#> 726  recruitment       National       14                  0 0.28461652
#> 727  recruitment       National       25                  0 0.23821686
#> 728  recruitment       National       40                  0 0.18404153
#> 729  recruitment       National       51                  0 0.14936286
#> 730  recruitment       National       18                  0 0.26668125
#> 731  recruitment       National       50                  0 0.15846137
#> 732  recruitment       National       17                  0 0.27038230
#> 733  recruitment       National       11                  0 0.29351167
#> 734  recruitment       National       22                  0 0.24039039
#> 735  recruitment       National       37                  0 0.18887983
#> 736  recruitment       National       48                  0 0.16024891
#> 737  recruitment       National       15                  0 0.27067303
#> 738  recruitment       National       26                  0 0.22687424
#> 739  recruitment       National       41                  0 0.17767882
#> 740  recruitment       National       52                  0 0.15125360
#> 741  recruitment       National       63                  0 0.12151937
#> 742  recruitment       National       34                  0 0.19564100
#> 743  recruitment       National        1                  0 0.34293957
#> 744  recruitment       National       12                  0 0.29416309
#> 745  recruitment       National       23                  0 0.23477573
#> 746  recruitment       National       38                  0 0.18646934
#> 747  recruitment       National       49                  0 0.16244494
#> 748  recruitment       National       60                  0 0.13335478
#> 749  recruitment       National       27                  0 0.22432977
#> 750  recruitment       National       42                  0 0.18270767
#> 751  recruitment       National       53                  0 0.14904862
#> 752  recruitment       National       64                  0 0.12023605
#> 753  recruitment       National       31                  0 0.20914639
#> 754  recruitment       National       46                  0 0.16858314
#> 755  recruitment       National       57                  0 0.13887593
#> 756  recruitment       National       68                  0 0.10700360
#> 757  recruitment       National       79                  0 0.09629799
#> 758  recruitment       National       90                  0 0.08163122
#> 759  recruitment       National       61                  0 0.12428156
#> 760  recruitment       National       28                  0 0.22265532
#> 761  recruitment       National       43                  0 0.17441283
#> 762  recruitment       National       54                  0 0.14347361
#> 763  recruitment       National       65                  0 0.11668308
#> 764  recruitment       National       32                  0 0.21394345
#> 765  recruitment       National       47                  0 0.15899683
#> 766  recruitment       National       58                  0 0.13829448
#> 767  recruitment       National       69                  0 0.11035068
#> 768  recruitment       National       80                  0 0.09130350
#> 769  recruitment       National       91                  0 0.07696494
#> 770  recruitment       National       62                  0 0.12941109
#> 771  recruitment       National       29                  0 0.22125138
#> 772  recruitment       National       44                  0 0.17131283
#> 773  recruitment       National       55                  0 0.14041543
#> 774  recruitment       National       66                  0 0.11907211
#> 775  recruitment       National       33                  0 0.20613802
#> 776  recruitment       National        0                  0 0.34568663
#> 777  recruitment       National       59                  0 0.13500345
#> 778  recruitment       National       70                  0 0.10893474
#> 779  recruitment       National       81                  0 0.08940191
#> 780  recruitment       National       92                  0 0.07940164
#> 781  recruitment       National       86                  0 0.08310234
#> 782  recruitment       National       30                  0 0.21421200
#> 783  recruitment       National       45                  0 0.17007525
#> 784  recruitment       National       56                  0 0.14021455
#> 785  recruitment       National       67                  0 0.11234947
#> 786  recruitment       National       78                  0 0.09741402
#> 787  recruitment       National       89                  0 0.08146700
#> 788  recruitment       National      100                  0 0.07236136
#> 789  recruitment       National       71                  0 0.10739676
#> 790  recruitment       National       82                  0 0.09197763
#> 791  recruitment       National       93                  0 0.07731209
#> 792  recruitment       National       87                  0 0.08434593
#> 793  recruitment       National       75                  0 0.09697221
#> 794  recruitment       National       73                  0 0.09984842
#> 795  recruitment       National       97                  0 0.07182483
#> 796  recruitment       National       74                  0 0.10223011
#> 797  recruitment       National       85                  0 0.08045331
#> 798  recruitment       National       94                  0 0.07234245
#> 799  recruitment       National       84                  0 0.09054180
#> 800  recruitment       National       72                  0 0.10735385
#> 801  recruitment       National       83                  0 0.09203906
#> 802  recruitment       National       98                  0 0.06879562
#> 803  recruitment       National       88                  0 0.08201827
#> 804  recruitment       National       76                  0 0.10073060
#> 805  recruitment       National       95                  0 0.07283352
#> 806  recruitment       National       96                  0 0.07432137
#> 807  recruitment       National       99                  0 0.06710919
#> 808  recruitment       National       77                  0 0.09754335
#> 809     survival       National        4                  0 0.87836006
#> 810     survival       National       69                  0 0.82890469
#> 811     survival       National       68                  0 0.83345873
#> 812     survival       National        2                  0 0.87563892
#> 813     survival       National        5                  0 0.87726255
#> 814     survival       National        7                  0 0.87289026
#> 815     survival       National       10                  0 0.87257538
#> 816     survival       National        6                  0 0.87448075
#> 817     survival       National        9                  0 0.86999269
#> 818     survival       National       11                  0 0.86915195
#> 819     survival       National       70                  0 0.82989376
#> 820     survival       National       37                  0 0.85270307
#> 821     survival       National       48                  0 0.85027155
#> 822     survival       National       36                  0 0.84558016
#> 823     survival       National        3                  0 0.87835704
#> 824     survival       National       35                  0 0.85358157
#> 825     survival       National        8                  0 0.87118796
#> 826     survival       National       67                  0 0.83117804
#> 827     survival       National       34                  0 0.85639780
#> 828     survival       National        1                  0 0.87864385
#> 829     survival       National       12                  0 0.86495902
#> 830     survival       National       71                  0 0.83344398
#> 831     survival       National       38                  0 0.85121590
#> 832     survival       National       49                  0 0.84407933
#> 833     survival       National       16                  0 0.86961689
#> 834     survival       National       75                  0 0.82348979
#> 835     survival       National       42                  0 0.84610578
#> 836     survival       National       53                  0 0.84296853
#> 837     survival       National       20                  0 0.86056188
#> 838     survival       National       79                  0 0.82360130
#> 839     survival       National       46                  0 0.84434769
#> 840     survival       National       13                  0 0.86631972
#> 841     survival       National       72                  0 0.83192837
#> 842     survival       National       39                  0 0.84871919
#> 843     survival       National       50                  0 0.84028776
#> 844     survival       National       17                  0 0.86646715
#> 845     survival       National       76                  0 0.82587884
#> 846     survival       National       43                  0 0.84661808
#> 847     survival       National       54                  0 0.84139966
#> 848     survival       National       21                  0 0.86560864
#> 849     survival       National       80                  0 0.81901238
#> 850     survival       National       47                  0 0.84565325
#> 851     survival       National       14                  0 0.86636567
#> 852     survival       National       73                  0 0.82988538
#> 853     survival       National       40                  0 0.84564436
#> 854     survival       National       51                  0 0.84343111
#> 855     survival       National       18                  0 0.86196106
#> 856     survival       National       77                  0 0.82605685
#> 857     survival       National       44                  0 0.85273625
#> 858     survival       National       55                  0 0.84210046
#> 859     survival       National       22                  0 0.86181903
#> 860     survival       National       81                  0 0.81994851
#> 861     survival       National       92                  0 0.81485195
#> 862     survival       National       15                  0 0.86710339
#> 863     survival       National       74                  0 0.82937472
#> 864     survival       National       41                  0 0.84912022
#> 865     survival       National       52                  0 0.84086343
#> 866     survival       National       19                  0 0.86549900
#> 867     survival       National       78                  0 0.82465042
#> 868     survival       National       45                  0 0.84296857
#> 869     survival       National       56                  0 0.83389579
#> 870     survival       National       23                  0 0.85846258
#> 871     survival       National       82                  0 0.82241856
#> 872     survival       National       93                  0 0.81455372
#> 873     survival       National       60                  0 0.83419261
#> 874     survival       National       27                  0 0.86165657
#> 875     survival       National       86                  0 0.81637236
#> 876     survival       National       97                  0 0.81338586
#> 877     survival       National       64                  0 0.83580273
#> 878     survival       National       31                  0 0.85496465
#> 879     survival       National       90                  0 0.81812062
#> 880     survival       National       57                  0 0.83995041
#> 881     survival       National       24                  0 0.85650053
#> 882     survival       National       83                  0 0.82206184
#> 883     survival       National       94                  0 0.81429234
#> 884     survival       National       61                  0 0.83337676
#> 885     survival       National       28                  0 0.85915358
#> 886     survival       National       87                  0 0.81811600
#> 887     survival       National       98                  0 0.81402348
#> 888     survival       National       65                  0 0.83576304
#> 889     survival       National       32                  0 0.85116536
#> 890     survival       National       91                  0 0.81754617
#> 891     survival       National       58                  0 0.83528607
#> 892     survival       National       25                  0 0.86023699
#> 893     survival       National       84                  0 0.81887331
#> 894     survival       National       95                  0 0.81258562
#> 895     survival       National       62                  0 0.83508212
#> 896     survival       National       29                  0 0.85851059
#> 897     survival       National       88                  0 0.82313299
#> 898     survival       National       99                  0 0.81098232
#> 899     survival       National       66                  0 0.83272871
#> 900     survival       National       33                  0 0.84933073
#> 901     survival       National        0                  0 0.88079279
#> 902     survival       National       59                  0 0.83573232
#> 903     survival       National       26                  0 0.86015616
#> 904     survival       National       85                  0 0.81723813
#> 905     survival       National       96                  0 0.81144136
#> 906     survival       National       63                  0 0.83612361
#> 907     survival       National       30                  0 0.85642552
#> 908     survival       National       89                  0 0.80888539
#> 909     survival       National      100                  0 0.81118764
#>            lower     upper probViable              Parameter
#> 102 0.1079292021 0.4992519      0.000   Expected recruitment
#> 103 0.0917897836 0.4697089      0.000   Expected recruitment
#> 104 0.0964830381 0.4817734      0.000   Expected recruitment
#> 105 0.1024410980 0.5108371      0.000   Expected recruitment
#> 106 0.1027134926 0.4636047      0.000   Expected recruitment
#> 107 0.1008570633 0.4753443      0.000   Expected recruitment
#> 108 0.0285835299 0.3289770      0.000   Expected recruitment
#> 109 0.1467354335 0.5498888      0.000   Expected recruitment
#> 110 0.0269936506 0.3134213      0.000   Expected recruitment
#> 111 0.0716852070 0.4306647      0.000   Expected recruitment
#> 112 0.0863207477 0.4653865      0.000   Expected recruitment
#> 113 0.0590021407 0.4095643      0.000   Expected recruitment
#> 114 0.1528989790 0.5659441      0.000   Expected recruitment
#> 115 0.0613396130 0.4106504      0.000   Expected recruitment
#> 116 0.1514599743 0.5644679      0.000   Expected recruitment
#> 117 0.0621509931 0.4003080      0.000   Expected recruitment
#> 118 0.1419461425 0.5403457      0.000   Expected recruitment
#> 119 0.0269273028 0.3153632      0.000   Expected recruitment
#> 120 0.1399293799 0.5787774      0.000   Expected recruitment
#> 121 0.1578342261 0.5901809      0.000   Expected recruitment
#> 122 0.0688254102 0.4327314      0.000   Expected recruitment
#> 123 0.0815482716 0.4578238      0.000   Expected recruitment
#> 124 0.0534263001 0.4141650      0.000   Expected recruitment
#> 125 0.1600303963 0.5641602      0.000   Expected recruitment
#> 126 0.0303942198 0.3294860      0.000   Expected recruitment
#> 127 0.0767278918 0.4262127      0.000   Expected recruitment
#> 128 0.0542271667 0.3931102      0.000   Expected recruitment
#> 129 0.1301989847 0.5283652      0.000   Expected recruitment
#> 130 0.0223546342 0.3042736      0.000   Expected recruitment
#> 131 0.0134585997 0.2808699      0.000   Expected recruitment
#> 132 0.1572961331 0.5727446      0.000   Expected recruitment
#> 133 0.0765932575 0.4230097      0.000   Expected recruitment
#> 134 0.0765572931 0.4377009      0.000   Expected recruitment
#> 135 0.0542180758 0.4074933      0.000   Expected recruitment
#> 136 0.0399771839 0.3447133      0.000   Expected recruitment
#> 137 0.0043308379 0.2108679      0.000   Expected recruitment
#> 138 0.0179408973 0.2863454      0.000   Expected recruitment
#> 139 0.0521196208 0.3894449      0.000   Expected recruitment
#> 140 0.1184423483 0.5234144      0.000   Expected recruitment
#> 141 0.0209437282 0.3139605      0.000   Expected recruitment
#> 142 0.0137429611 0.2710828      0.000   Expected recruitment
#> 143 0.0422835872 0.3681245      0.000   Expected recruitment
#> 144 0.1081604077 0.5111332      0.000   Expected recruitment
#> 145 0.0857453404 0.4443990      0.000   Expected recruitment
#> 146 0.0553532697 0.3934951      0.000   Expected recruitment
#> 147 0.0368188550 0.3558516      0.000   Expected recruitment
#> 148 0.0022527271 0.2180726      0.000   Expected recruitment
#> 149 0.0153163679 0.2787104      0.000   Expected recruitment
#> 150 0.0539670048 0.4129100      0.000   Expected recruitment
#> 151 0.1312261478 0.5367737      0.000   Expected recruitment
#> 152 0.0230813449 0.3036697      0.000   Expected recruitment
#> 153 0.0106256658 0.2714266      0.000   Expected recruitment
#> 154 0.0478921067 0.3702977      0.000   Expected recruitment
#> 155 0.1086552245 0.5003140      0.000   Expected recruitment
#> 156 0.0805400479 0.4475746      0.000   Expected recruitment
#> 157 0.0509641502 0.3912621      0.000   Expected recruitment
#> 158 0.0340193033 0.3474552      0.000   Expected recruitment
#> 159 0.0026534846 0.2136226      0.000   Expected recruitment
#> 160 0.0158285123 0.2740653      0.000   Expected recruitment
#> 161 0.0520628272 0.3840246      0.000   Expected recruitment
#> 162 0.1134617393 0.5109438      0.000   Expected recruitment
#> 163 0.0236271748 0.3088848      0.000   Expected recruitment
#> 164 0.0129076344 0.2653648      0.000   Expected recruitment
#> 165 0.0414574260 0.3582559      0.000   Expected recruitment
#> 166 0.0043692231 0.2202099      0.000   Expected recruitment
#> 167 0.0194086473 0.2800404      0.000   Expected recruitment
#> 168 0.0114636040 0.2678823      0.000   Expected recruitment
#> 169 0.0369441107 0.3446481      0.000   Expected recruitment
#> 170 0.0028107248 0.2075875      0.000   Expected recruitment
#> 171 0.0301923820 0.3264404      0.000   Expected recruitment
#> 172 0.0527002464 0.3744919      0.000   Expected recruitment
#> 173 0.1274965092 0.4964689      0.000   Expected recruitment
#> 174 0.0184171506 0.3088619      0.000   Expected recruitment
#> 175 0.0132799550 0.2671229      0.000   Expected recruitment
#> 176 0.0062599861 0.2201083      0.000   Expected recruitment
#> 177 0.0031973520 0.2377223      0.000   Expected recruitment
#> 178 0.0310123370 0.3456386      0.000   Expected recruitment
#> 179 0.0110170031 0.2610748      0.000   Expected recruitment
#> 180 0.0352358432 0.3648340      0.000   Expected recruitment
#> 181 0.0079698293 0.2442029      0.000   Expected recruitment
#> 182 0.0290483496 0.3175212      0.000   Expected recruitment
#> 183 0.0073663445 0.2525420      0.000   Expected recruitment
#> 184 0.0134116285 0.2645967      0.000   Expected recruitment
#> 185 0.0213128612 0.2968342      0.000   Expected recruitment
#> 186 0.0127971020 0.2595228      0.000   Expected recruitment
#> 187 0.0077481254 0.2293643      0.000   Expected recruitment
#> 188 0.0090232223 0.2466575      0.000   Expected recruitment
#> 189 0.0325875682 0.3446916      0.000   Expected recruitment
#> 190 0.0103094348 0.2476015      0.000   Expected recruitment
#> 191 0.0094250518 0.2500032      0.000   Expected recruitment
#> 192 0.0065244567 0.2331529      0.000   Expected recruitment
#> 193 0.0052214046 0.2316132      0.000   Expected recruitment
#> 194 0.0067845059 0.2256576      0.000   Expected recruitment
#> 195 0.0052232574 0.2244337      0.000   Expected recruitment
#> 196 0.0178918348 0.2890211      0.000   Expected recruitment
#> 197 0.0055072090 0.2478620      0.000   Expected recruitment
#> 198 0.0062444851 0.2392479      0.000   Expected recruitment
#> 199 0.0049059126 0.2283134      0.000   Expected recruitment
#> 200 0.0298888840 0.3455622      0.000   Expected recruitment
#> 201 0.0050795511 0.2380352      0.000   Expected recruitment
#> 202 0.0047398310 0.2132379      0.000   Expected recruitment
#> 203 0.7764765312 0.9368256      0.000      Expected survival
#> 204 0.7677626596 0.9321966      0.000      Expected survival
#> 205 0.7670938434 0.9326586      0.000      Expected survival
#> 206 0.7649682282 0.9362196      0.000      Expected survival
#> 207 0.7607393655 0.9322969      0.000      Expected survival
#> 208 0.7575536837 0.9302476      0.000      Expected survival
#> 209 0.7742973786 0.9342918      0.000      Expected survival
#> 210 0.7793200130 0.9462082      0.000      Expected survival
#> 211 0.7750512838 0.9390920      0.000      Expected survival
#> 212 0.7887786628 0.9414727      0.000      Expected survival
#> 213 0.7664628966 0.9367892      0.000      Expected survival
#> 214 0.7751239607 0.9373358      0.000      Expected survival
#> 215 0.7590425859 0.9328201      0.000      Expected survival
#> 216 0.7804853328 0.9384732      0.000      Expected survival
#> 217 0.7771332023 0.9473148      0.000      Expected survival
#> 218 0.7771819556 0.9385618      0.000      Expected survival
#> 219 0.7546818682 0.9325140      0.000      Expected survival
#> 220 0.7736149149 0.9386350      0.000      Expected survival
#> 221 0.7478853979 0.9212247      0.000      Expected survival
#> 222 0.7780365433 0.9421894      0.000      Expected survival
#> 223 0.7593945749 0.9310624      0.000      Expected survival
#> 224 0.7386633595 0.9152090      0.000      Expected survival
#> 225 0.7701665818 0.9392505      0.000      Expected survival
#> 226 0.7660207050 0.9391414      0.000      Expected survival
#> 227 0.7737066662 0.9391721      0.000      Expected survival
#> 228 0.7403094645 0.9223487      0.000      Expected survival
#> 229 0.7496116449 0.9307859      0.000      Expected survival
#> 230 0.7699817570 0.9333351      0.000      Expected survival
#> 231 0.7374525370 0.9180682      0.000      Expected survival
#> 232 0.7356680570 0.9204262      0.000      Expected survival
#> 233 0.7257090596 0.9160972      0.000      Expected survival
#> 234 0.7498584632 0.9316418      0.000      Expected survival
#> 235 0.7354666675 0.9139755      0.000      Expected survival
#> 236 0.7756184435 0.9345028      0.000      Expected survival
#> 237 0.7665994821 0.9317062      0.000      Expected survival
#> 238 0.7476784394 0.9246239      0.000      Expected survival
#> 239 0.7421245634 0.9254019      0.000      Expected survival
#> 240 0.7556439428 0.9237395      0.000      Expected survival
#> 241 0.7618155732 0.9338709      0.000      Expected survival
#> 242 0.7254895345 0.9144225      0.000      Expected survival
#> 243 0.7478100353 0.9219176      0.000      Expected survival
#> 244 0.7344940943 0.9215556      0.000      Expected survival
#> 245 0.7503584487 0.9303778      0.000      Expected survival
#> 246 0.7964237757 0.9432318      0.000      Expected survival
#> 247 0.7375663767 0.9192601      0.000      Expected survival
#> 248 0.7303578776 0.9092833      0.000      Expected survival
#> 249 0.7603588287 0.9291921      0.000      Expected survival
#> 250 0.7359938936 0.9177685      0.000      Expected survival
#> 251 0.7031782567 0.9009843      0.000      Expected survival
#> 252 0.7581641987 0.9316208      0.000      Expected survival
#> 253 0.7372239886 0.9214828      0.000      Expected survival
#> 254 0.7374297273 0.9227291      0.000      Expected survival
#> 255 0.7228482518 0.9142055      0.000      Expected survival
#> 256 0.7198246793 0.9033486      0.000      Expected survival
#> 257 0.7803662429 0.9451407      0.000      Expected survival
#> 258 0.7056438727 0.8980982      0.000      Expected survival
#> 259 0.7226650445 0.9162240      0.000      Expected survival
#> 260 0.7460035365 0.9251543      0.000      Expected survival
#> 261 0.7062084967 0.9009861      0.000      Expected survival
#> 262 0.7021711343 0.8989137      0.000      Expected survival
#> 263 0.7342839171 0.9173306      0.000      Expected survival
#> 264 0.7868749427 0.9468515      0.000      Expected survival
#> 265 0.7400049134 0.9215816      0.000      Expected survival
#> 266 0.7263904966 0.9170449      0.000      Expected survival
#> 267 0.7270289318 0.9100210      0.000      Expected survival
#> 268 0.7075233423 0.8953954      0.000      Expected survival
#> 269 0.7074391733 0.9050399      0.000      Expected survival
#> 270 0.7240074818 0.9151109      0.000      Expected survival
#> 271 0.7505121496 0.9277160      0.000      Expected survival
#> 272 0.7048231963 0.9041788      0.000      Expected survival
#> 273 0.7514275818 0.9233391      0.000      Expected survival
#> 274 0.7244607972 0.9092558      0.000      Expected survival
#> 275 0.7804902661 0.9413631      0.000      Expected survival
#> 276 0.7472258355 0.9170104      0.000      Expected survival
#> 277 0.7320890852 0.9100785      0.000      Expected survival
#> 278 0.7121010518 0.9067889      0.000      Expected survival
#> 279 0.7023164487 0.9073356      0.000      Expected survival
#> 280 0.7039493756 0.8958275      0.000      Expected survival
#> 281 0.7205750148 0.9113394      0.000      Expected survival
#> 282 0.7481205107 0.9253422      0.000      Expected survival
#> 283 0.7137651321 0.9004940      0.000      Expected survival
#> 284 0.7428566823 0.9237137      0.000      Expected survival
#> 285 0.7204971378 0.9085328      0.000      Expected survival
#> 286 0.7001651959 0.9029360      0.000      Expected survival
#> 287 0.7047271231 0.9033797      0.000      Expected survival
#> 288 0.7238851361 0.9102440      0.000      Expected survival
#> 289 0.7172922477 0.9106215      0.000      Expected survival
#> 290 0.7060233608 0.9003150      0.000      Expected survival
#> 291 0.7435235080 0.9179784      0.000      Expected survival
#> 292 0.7173297519 0.9049676      0.000      Expected survival
#> 293 0.7506322498 0.9287135      0.000      Expected survival
#> 294 0.7181390032 0.9106699      0.000      Expected survival
#> 295 0.7076567818 0.9060399      0.000      Expected survival
#> 296 0.7228368631 0.9102990      0.000      Expected survival
#> 297 0.6973424673 0.9029509      0.000      Expected survival
#> 298 0.7437913797 0.9241463      0.000      Expected survival
#> 299 0.7116612184 0.9033326      0.000      Expected survival
#> 300 0.7194786603 0.9126057      0.000      Expected survival
#> 301 0.7111699262 0.9058747      0.000      Expected survival
#> 302 0.7134365332 0.9064586      0.000      Expected survival
#> 303 0.7093022334 0.9075587      0.000      Expected survival
#> 304 0.0291365269 0.3259661      0.000   Adjusted recruitment
#> 305 0.0372886949 0.3614452      0.000   Adjusted recruitment
#> 306 0.0374588173 0.3549290      0.000   Adjusted recruitment
#> 307 0.0338796798 0.3430862      0.000   Adjusted recruitment
#> 308 0.0401521590 0.3523841      0.000   Adjusted recruitment
#> 309 0.0147428374 0.2680595      0.000   Adjusted recruitment
#> 310 0.0318274357 0.3431285      0.000   Adjusted recruitment
#> 311 0.0344438286 0.3522607      0.000   Adjusted recruitment
#> 312 0.0184106706 0.2512852      0.000   Adjusted recruitment
#> 313 0.0258479161 0.3051383      0.000   Adjusted recruitment
#> 314 0.0330123841 0.3231513      0.000   Adjusted recruitment
#> 315 0.0318784915 0.3474716      0.000   Adjusted recruitment
#> 316 0.0182341431 0.2831825      0.000   Adjusted recruitment
#> 317 0.0368643080 0.3368525      0.000   Adjusted recruitment
#> 318 0.0276344174 0.3301052      0.000   Adjusted recruitment
#> 319 0.0156033892 0.2404567      0.000   Adjusted recruitment
#> 320 0.0262906444 0.3263392      0.000   Adjusted recruitment
#> 321 0.0306665453 0.3515423      0.000   Adjusted recruitment
#> 322 0.0431991071 0.3753494      0.000   Adjusted recruitment
#> 323 0.0154699732 0.2547762      0.000   Adjusted recruitment
#> 324 0.0360544176 0.3666219      0.000   Adjusted recruitment
#> 325 0.0256496128 0.3235808      0.000   Adjusted recruitment
#> 326 0.0109792054 0.2358665      0.000   Adjusted recruitment
#> 327 0.0416944556 0.3678866      0.000   Adjusted recruitment
#> 328 0.0100119852 0.2087433      0.000   Adjusted recruitment
#> 329 0.0307883092 0.3383189      0.000   Adjusted recruitment
#> 330 0.0173862960 0.2409065      0.000   Adjusted recruitment
#> 331 0.0214916022 0.3118042      0.000   Adjusted recruitment
#> 332 0.0259681979 0.3406018      0.000   Adjusted recruitment
#> 333 0.0337056498 0.3546736      0.000   Adjusted recruitment
#> 334 0.0142508552 0.2364753      0.000   Adjusted recruitment
#> 335 0.0108164028 0.2083372      0.000   Adjusted recruitment
#> 336 0.0067029735 0.1726442      0.000   Adjusted recruitment
#> 337 0.0099160033 0.2141138      0.000   Adjusted recruitment
#> 338 0.0114674850 0.2393278      0.000   Adjusted recruitment
#> 339 0.0082378817 0.2074148      0.000   Adjusted recruitment
#> 340 0.0232183191 0.2965336      0.000   Adjusted recruitment
#> 341 0.0147475319 0.2614416      0.000   Adjusted recruitment
#> 342 0.0222812884 0.3156669      0.000   Adjusted recruitment
#> 343 0.0075516990 0.2017516      0.000   Adjusted recruitment
#> 344 0.0112326772 0.2278311      0.000   Adjusted recruitment
#> 345 0.0131775701 0.2437441      0.000   Adjusted recruitment
#> 346 0.0103287626 0.2225118      0.000   Adjusted recruitment
#> 347 0.0075940879 0.1816959      0.000   Adjusted recruitment
#> 348 0.0182729147 0.2843950      0.000   Adjusted recruitment
#> 349 0.0121558142 0.2371929      0.000   Adjusted recruitment
#> 350 0.0094759246 0.2169815      0.000   Adjusted recruitment
#> 351 0.0247493000 0.2937443      0.000   Adjusted recruitment
#> 352 0.0036464844 0.1676213      0.000   Adjusted recruitment
#> 353 0.0016103896 0.1412499      0.000   Adjusted recruitment
#> 354 0.0088345086 0.1761022      0.000   Adjusted recruitment
#> 355 0.0112367747 0.2361410      0.000   Adjusted recruitment
#> 356 0.0145342240 0.2461458      0.000   Adjusted recruitment
#> 357 0.0104918984 0.2203660      0.000   Adjusted recruitment
#> 358 0.0055465626 0.1832451      0.000   Adjusted recruitment
#> 359 0.0207993611 0.2717176      0.000   Adjusted recruitment
#> 360 0.0109494980 0.2312599      0.000   Adjusted recruitment
#> 361 0.0089103495 0.2152105      0.000   Adjusted recruitment
#> 362 0.0230267126 0.2931694      0.000   Adjusted recruitment
#> 363 0.0037716695 0.1548838      0.000   Adjusted recruitment
#> 364 0.0021765448 0.1325365      0.000   Adjusted recruitment
#> 365 0.0076466715 0.2061147      0.000   Adjusted recruitment
#> 366 0.0207022357 0.2923856      0.000   Adjusted recruitment
#> 367 0.0136544836 0.2472347      0.000   Adjusted recruitment
#> 368 0.0017168044 0.1338593      0.000   Adjusted recruitment
#> 369 0.0061125935 0.1722843      0.000   Adjusted recruitment
#> 370 0.0172808816 0.2809766      0.000   Adjusted recruitment
#> 371 0.0021932805 0.1344441      0.000   Adjusted recruitment
#> 372 0.0087976551 0.2130619      0.000   Adjusted recruitment
#> 373 0.0207331605 0.2876902      0.000   Adjusted recruitment
#> 374 0.0031437280 0.1551021      0.000   Adjusted recruitment
#> 375 0.0017309587 0.1476325      0.000   Adjusted recruitment
#> 376 0.0026961761 0.1426451      0.000   Adjusted recruitment
#> 377 0.0188571624 0.2775938      0.000   Adjusted recruitment
#> 378 0.0021350739 0.1382968      0.000   Adjusted recruitment
#> 379 0.0013833171 0.1375811      0.000   Adjusted recruitment
#> 380 0.0062709700 0.1712377      0.000   Adjusted recruitment
#> 381 0.0039660508 0.1724432      0.000   Adjusted recruitment
#> 382 0.0026064470 0.1414033      0.000   Adjusted recruitment
#> 383 0.0016932690 0.1295730      0.000   Adjusted recruitment
#> 384 0.0223064030 0.2964540      0.000   Adjusted recruitment
#> 385 0.0028658492 0.1530848      0.000   Adjusted recruitment
#> 386 0.0018791123 0.1473358      0.000   Adjusted recruitment
#> 387 0.0025229909 0.1571206      0.000   Adjusted recruitment
#> 388 0.0048702852 0.1563059      0.000   Adjusted recruitment
#> 389 0.0048159185 0.1636690      0.000   Adjusted recruitment
#> 390 0.0009625138 0.1298750      0.000   Adjusted recruitment
#> 391 0.0055003941 0.1750218      0.000   Adjusted recruitment
#> 392 0.0024200709 0.1643280      0.000   Adjusted recruitment
#> 393 0.0050942607 0.1485761      0.000   Adjusted recruitment
#> 394 0.0030968345 0.1553783      0.000   Adjusted recruitment
#> 395 0.0218210790 0.2712075      0.000   Adjusted recruitment
#> 396 0.0052665551 0.1800770      0.000   Adjusted recruitment
#> 397 0.0016530640 0.1278925      0.000   Adjusted recruitment
#> 398 0.0053941563 0.1767940      0.000   Adjusted recruitment
#> 399 0.0043456861 0.1761507      0.000   Adjusted recruitment
#> 400 0.0052409665 0.1610585      0.000   Adjusted recruitment
#> 401 0.0011907048 0.1198052      0.000   Adjusted recruitment
#> 402 0.0039654822 0.1595759      0.000   Adjusted recruitment
#> 403 0.0010247933 0.1255365      0.000   Adjusted recruitment
#> 404 0.0040631111 0.1750140      0.000   Adjusted recruitment
#> 405 1.0000000000 1.0000000      1.000                      c
#> 406 1.0000000000 1.0000000      1.000                      c
#> 407 1.0000000000 1.0000000      1.000                      c
#> 408 1.0000000000 1.0000000      1.000                      c
#> 409 1.0000000000 1.0000000      1.000                      c
#> 410 1.0000000000 1.0000000      1.000                      c
#> 411 1.0000000000 1.0000000      1.000                      c
#> 412 1.0000000000 1.0000000      1.000                      c
#> 413 1.0000000000 1.0000000      1.000                      c
#> 414 1.0000000000 1.0000000      1.000                      c
#> 415 1.0000000000 1.0000000      1.000                      c
#> 416 1.0000000000 1.0000000      1.000                      c
#> 417 1.0000000000 1.0000000      1.000                      c
#> 418 1.0000000000 1.0000000      1.000                      c
#> 419 1.0000000000 1.0000000      1.000                      c
#> 420 1.0000000000 1.0000000      1.000                      c
#> 421 1.0000000000 1.0000000      1.000                      c
#> 422 1.0000000000 1.0000000      1.000                      c
#> 423 1.0000000000 1.0000000      1.000                      c
#> 424 1.0000000000 1.0000000      1.000                      c
#> 425 1.0000000000 1.0000000      1.000                      c
#> 426 1.0000000000 1.0000000      1.000                      c
#> 427 1.0000000000 1.0000000      1.000                      c
#> 428 1.0000000000 1.0000000      1.000                      c
#> 429 1.0000000000 1.0000000      1.000                      c
#> 430 1.0000000000 1.0000000      1.000                      c
#> 431 1.0000000000 1.0000000      1.000                      c
#> 432 1.0000000000 1.0000000      1.000                      c
#> 433 1.0000000000 1.0000000      1.000                      c
#> 434 1.0000000000 1.0000000      1.000                      c
#> 435 1.0000000000 1.0000000      1.000                      c
#> 436 1.0000000000 1.0000000      1.000                      c
#> 437 1.0000000000 1.0000000      1.000                      c
#> 438 1.0000000000 1.0000000      1.000                      c
#> 439 1.0000000000 1.0000000      1.000                      c
#> 440 1.0000000000 1.0000000      1.000                      c
#> 441 1.0000000000 1.0000000      1.000                      c
#> 442 1.0000000000 1.0000000      1.000                      c
#> 443 1.0000000000 1.0000000      1.000                      c
#> 444 1.0000000000 1.0000000      1.000                      c
#> 445 1.0000000000 1.0000000      1.000                      c
#> 446 1.0000000000 1.0000000      1.000                      c
#> 447 1.0000000000 1.0000000      1.000                      c
#> 448 1.0000000000 1.0000000      1.000                      c
#> 449 1.0000000000 1.0000000      1.000                      c
#> 450 1.0000000000 1.0000000      1.000                      c
#> 451 1.0000000000 1.0000000      1.000                      c
#> 452 1.0000000000 1.0000000      1.000                      c
#> 453 1.0000000000 1.0000000      1.000                      c
#> 454 1.0000000000 1.0000000      1.000                      c
#> 455 1.0000000000 1.0000000      1.000                      c
#> 456 1.0000000000 1.0000000      1.000                      c
#> 457 1.0000000000 1.0000000      1.000                      c
#> 458 1.0000000000 1.0000000      1.000                      c
#> 459 1.0000000000 1.0000000      1.000                      c
#> 460 1.0000000000 1.0000000      1.000                      c
#> 461 1.0000000000 1.0000000      1.000                      c
#> 462 1.0000000000 1.0000000      1.000                      c
#> 463 1.0000000000 1.0000000      1.000                      c
#> 464 1.0000000000 1.0000000      1.000                      c
#> 465 1.0000000000 1.0000000      1.000                      c
#> 466 1.0000000000 1.0000000      1.000                      c
#> 467 1.0000000000 1.0000000      1.000                      c
#> 468 1.0000000000 1.0000000      1.000                      c
#> 469 1.0000000000 1.0000000      1.000                      c
#> 470 1.0000000000 1.0000000      1.000                      c
#> 471 1.0000000000 1.0000000      1.000                      c
#> 472 1.0000000000 1.0000000      1.000                      c
#> 473 1.0000000000 1.0000000      1.000                      c
#> 474 1.0000000000 1.0000000      1.000                      c
#> 475 1.0000000000 1.0000000      1.000                      c
#> 476 1.0000000000 1.0000000      1.000                      c
#> 477 1.0000000000 1.0000000      1.000                      c
#> 478 1.0000000000 1.0000000      1.000                      c
#> 479 1.0000000000 1.0000000      1.000                      c
#> 480 1.0000000000 1.0000000      1.000                      c
#> 481 1.0000000000 1.0000000      1.000                      c
#> 482 1.0000000000 1.0000000      1.000                      c
#> 483 1.0000000000 1.0000000      1.000                      c
#> 484 1.0000000000 1.0000000      1.000                      c
#> 485 1.0000000000 1.0000000      1.000                      c
#> 486 1.0000000000 1.0000000      1.000                      c
#> 487 1.0000000000 1.0000000      1.000                      c
#> 488 1.0000000000 1.0000000      1.000                      c
#> 489 1.0000000000 1.0000000      1.000                      c
#> 490 1.0000000000 1.0000000      1.000                      c
#> 491 1.0000000000 1.0000000      1.000                      c
#> 492 1.0000000000 1.0000000      1.000                      c
#> 493 1.0000000000 1.0000000      1.000                      c
#> 494 1.0000000000 1.0000000      1.000                      c
#> 495 1.0000000000 1.0000000      1.000                      c
#> 496 1.0000000000 1.0000000      1.000                      c
#> 497 1.0000000000 1.0000000      1.000                      c
#> 498 1.0000000000 1.0000000      1.000                      c
#> 499 1.0000000000 1.0000000      1.000                      c
#> 500 1.0000000000 1.0000000      1.000                      c
#> 501 1.0000000000 1.0000000      1.000                      c
#> 502 1.0000000000 1.0000000      1.000                      c
#> 503 1.0000000000 1.0000000      1.000                      c
#> 504 1.0000000000 1.0000000      1.000                      c
#> 505 1.0000000000 1.0000000      1.000                      c
#> 506 0.7389250000 1.2260500      0.502 Population growth rate
#> 507 0.7850000000 1.2570250      0.587 Population growth rate
#> 508 0.7619250000 1.2330250      0.488 Population growth rate
#> 509 0.7529750000 1.2100500      0.449 Population growth rate
#> 510 0.7279500000 1.1620250      0.333 Population growth rate
#> 511 0.7670000000 1.2230000      0.504 Population growth rate
#> 512 0.7729750000 1.2620000      0.582 Population growth rate
#> 513 0.6960000000 1.0860750      0.174 Population growth rate
#> 514 0.7449750000 1.2150500      0.467 Population growth rate
#> 515 0.7389750000 1.2320000      0.508 Population growth rate
#> 516 0.7668250000 1.2391750      0.542 Population growth rate
#> 517 0.7489750000 1.1890250      0.427 Population growth rate
#> 518 0.7279750000 1.1610000      0.385 Population growth rate
#> 519 0.7500000000 1.2190250      0.450 Population growth rate
#> 520 0.7219500000 1.1700250      0.340 Population growth rate
#> 521 0.7770000000 1.2800750      0.601 Population growth rate
#> 522 0.7619750000 1.2072000      0.511 Population growth rate
#> 523 0.7680000000 1.2482750      0.574 Population growth rate
#> 524 0.7030000000 1.0820500      0.172 Population growth rate
#> 525 0.7480000000 1.1630000      0.373 Population growth rate
#> 526 0.7879750000 1.2970750      0.643 Population growth rate
#> 527 0.7529750000 1.2261250      0.527 Population growth rate
#> 528 0.7519250000 1.1990250      0.399 Population growth rate
#> 529 0.7169500000 1.1690500      0.308 Population growth rate
#> 530 0.7719750000 1.2670250      0.613 Population growth rate
#> 531 0.7180000000 1.1510250      0.342 Population growth rate
#> 532 0.7800000000 1.2540750      0.599 Population growth rate
#> 533 0.7560000000 1.2340500      0.496 Population growth rate
#> 534 0.7759250000 1.2580250      0.538 Population growth rate
#> 535 0.6940000000 1.0731000      0.161 Population growth rate
#> 536 0.6919000000 1.0580250      0.124 Population growth rate
#> 537 0.7880000000 1.2830000      0.630 Population growth rate
#> 538 0.7009250000 1.0961500      0.194 Population growth rate
#> 539 0.7319500000 1.1750750      0.410 Population growth rate
#> 540 0.7150000000 1.1781000      0.333 Population growth rate
#> 541 0.7060000000 1.1110000      0.234 Population growth rate
#> 542 0.7019750000 1.0790500      0.158 Population growth rate
#> 543 0.7448250000 1.1790250      0.394 Population growth rate
#> 544 0.7089000000 1.0900000      0.179 Population growth rate
#> 545 0.7509500000 1.2371500      0.534 Population growth rate
#> 546 0.6990000000 1.0880250      0.151 Population growth rate
#> 547 0.6909750000 1.0760250      0.121 Population growth rate
#> 548 0.7170000000 1.1180750      0.262 Population growth rate
#> 549 0.6939750000 1.1050250      0.198 Population growth rate
#> 550 0.7429000000 1.1650250      0.390 Population growth rate
#> 551 0.7289250000 1.1440750      0.309 Population growth rate
#> 552 0.7130000000 1.1170250      0.240 Population growth rate
#> 553 0.6709250000 1.0230250      0.054 Population growth rate
#> 554 0.6899750000 1.0530500      0.116 Population growth rate
#> 555 0.7209500000 1.1610500      0.284 Population growth rate
#> 556 0.6959750000 1.1242500      0.216 Population growth rate
#> 557 0.7000000000 1.0660250      0.138 Population growth rate
#> 558 0.6819000000 1.0490000      0.121 Population growth rate
#> 559 0.7179500000 1.1300250      0.252 Population growth rate
#> 560 0.7059250000 1.0940500      0.195 Population growth rate
#> 561 0.7549000000 1.1841000      0.396 Population growth rate
#> 562 0.7309750000 1.1370250      0.269 Population growth rate
#> 563 0.7040000000 1.1210000      0.227 Population growth rate
#> 564 0.6660000000 1.0030000      0.041 Population growth rate
#> 565 0.6849750000 1.0640250      0.118 Population growth rate
#> 566 0.7130000000 1.1320500      0.270 Population growth rate
#> 567 0.7128500000 1.1020000      0.204 Population growth rate
#> 568 0.6939500000 1.0770250      0.159 Population growth rate
#> 569 0.6750000000 1.0430250      0.108 Population growth rate
#> 570 0.7090000000 1.1190250      0.284 Population growth rate
#> 571 0.6559500000 1.0130000      0.047 Population growth rate
#> 572 0.7479250000 1.1690250      0.378 Population growth rate
#> 573 0.7190000000 1.1530500      0.302 Population growth rate
#> 574 0.7180000000 1.1240000      0.237 Population growth rate
#> 575 0.6639500000 1.0070250      0.039 Population growth rate
#> 576 0.6609000000 1.0250500      0.057 Population growth rate
#> 577 0.7169750000 1.1150000      0.278 Population growth rate
#> 578 0.6629750000 1.0101750      0.059 Population growth rate
#> 579 0.6890000000 1.0711000      0.137 Population growth rate
#> 580 0.6909500000 1.0500000      0.125 Population growth rate
#> 581 0.6759750000 1.0300250      0.066 Population growth rate
#> 582 0.6639750000 1.0080250      0.045 Population growth rate
#> 583 0.6909750000 1.0530000      0.110 Population growth rate
#> 584 0.7099750000 1.1451000      0.303 Population growth rate
#> 585 0.7049750000 1.1120250      0.241 Population growth rate
#> 586 0.6580000000 1.0100250      0.050 Population growth rate
#> 587 0.6649750000 1.0200000      0.057 Population growth rate
#> 588 0.6819750000 1.0350250      0.079 Population growth rate
#> 589 0.6810000000 1.0440250      0.095 Population growth rate
#> 590 0.6829750000 1.0790500      0.145 Population growth rate
#> 591 0.6839500000 1.0410000      0.080 Population growth rate
#> 592 0.6770000000 1.0230500      0.069 Population growth rate
#> 593 0.6718750000 1.0300250      0.075 Population growth rate
#> 594 0.6759750000 1.0140500      0.052 Population growth rate
#> 595 0.6740000000 1.0340250      0.090 Population growth rate
#> 596 0.6750000000 1.0210250      0.058 Population growth rate
#> 597 0.6820000000 1.0280750      0.081 Population growth rate
#> 598 0.6750000000 1.0520750      0.108 Population growth rate
#> 599 0.6749750000 1.0410750      0.091 Population growth rate
#> 600 0.6879750000 1.0370750      0.089 Population growth rate
#> 601 0.6969500000 1.0820000      0.141 Population growth rate
#> 602 0.6729750000 1.0280250      0.073 Population growth rate
#> 603 0.6649250000 1.0241000      0.072 Population growth rate
#> 604 0.6789250000 1.0280500      0.070 Population growth rate
#> 605 0.6569750000 1.0210250      0.051 Population growth rate
#> 606 0.6889500000 1.0600000      0.125 Population growth rate
#> 607 0.8497623501 1.1063926      0.396   Expected growth rate
#> 608 0.8499963761 1.1066572      0.452   Expected growth rate
#> 609 0.8498369157 1.0998957      0.437   Expected growth rate
#> 610 0.8765287857 1.1468541      0.659   Expected growth rate
#> 611 0.8334478460 1.0924659      0.383   Expected growth rate
#> 612 0.8772507543 1.1233477      0.542   Expected growth rate
#> 613 0.8446472522 1.0954428      0.362   Expected growth rate
#> 614 0.8653067529 1.1210310      0.577   Expected growth rate
#> 615 0.8712898126 1.1367136      0.560   Expected growth rate
#> 616 0.7705641365 1.0017345      0.037   Expected growth rate
#> 617 0.8473028721 1.0965020      0.376   Expected growth rate
#> 618 0.8892564501 1.1622660      0.731   Expected growth rate
#> 619 0.8749688742 1.1258111      0.519   Expected growth rate
#> 620 0.8357601949 1.0858981      0.321   Expected growth rate
#> 621 0.8902074475 1.1344661      0.623   Expected growth rate
#> 622 0.8788748849 1.1437362      0.657   Expected growth rate
#> 623 0.8597603706 1.1185574      0.489   Expected growth rate
#> 624 0.8933944733 1.1523723      0.693   Expected growth rate
#> 625 0.8117391835 1.0697568      0.191   Expected growth rate
#> 626 0.8736554180 1.1373699      0.563   Expected growth rate
#> 627 0.7685435367 1.0087162      0.047   Expected growth rate
#> 628 0.8175970439 1.0665082      0.242   Expected growth rate
#> 629 0.8996740043 1.1648498      0.703   Expected growth rate
#> 630 0.8611966425 1.1169267      0.504   Expected growth rate
#> 631 0.8332474165 1.0956301      0.350   Expected growth rate
#> 632 0.8069158635 1.0568854      0.179   Expected growth rate
#> 633 0.7946400471 1.0263261      0.110   Expected growth rate
#> 634 0.7663586317 1.0110055      0.055   Expected growth rate
#> 635 0.8368504475 1.0741273      0.274   Expected growth rate
#> 636 0.8073518025 1.0601743      0.173   Expected growth rate
#> 637 0.7754736107 1.0294159      0.079   Expected growth rate
#> 638 0.7637995272 1.0051864      0.040   Expected growth rate
#> 639 0.8118227777 1.0704157      0.203   Expected growth rate
#> 640 0.8877960422 1.1504288      0.697   Expected growth rate
#> 641 0.8605938915 1.1168016      0.510   Expected growth rate
#> 642 0.8470533561 1.0835754      0.325   Expected growth rate
#> 643 0.8041146050 1.0530396      0.169   Expected growth rate
#> 644 0.7955893736 1.0418377      0.096   Expected growth rate
#> 645 0.7823120515 1.0036874      0.052   Expected growth rate
#> 646 0.8330876063 1.0774340      0.255   Expected growth rate
#> 647 0.7956419940 1.0461848      0.151   Expected growth rate
#> 648 0.7924236778 1.0315776      0.080   Expected growth rate
#> 649 0.7732384340 1.0113616      0.053   Expected growth rate
#> 650 0.8071423346 1.0686565      0.230   Expected growth rate
#> 651 0.9014365932 1.1661562      0.739   Expected growth rate
#> 652 0.8565920098 1.1123744      0.499   Expected growth rate
#> 653 0.8316013545 1.0749620      0.311   Expected growth rate
#> 654 0.8207276052 1.0617557      0.191   Expected growth rate
#> 655 0.7740535589 1.0251429      0.079   Expected growth rate
#> 656 0.7189169007 0.9432136      0.003   Expected growth rate
#> 657 0.8244690929 1.0674089      0.232   Expected growth rate
#> 658 0.8053907056 1.0385995      0.123   Expected growth rate
#> 659 0.7842398649 1.0258746      0.077   Expected growth rate
#> 660 0.7645938339 0.9993018      0.034   Expected growth rate
#> 661 0.7543756141 0.9800034      0.014   Expected growth rate
#> 662 0.7856696262 1.0308428      0.113   Expected growth rate
#> 663 0.7764632617 1.0136866      0.057   Expected growth rate
#> 664 0.8304813491 1.0754242      0.285   Expected growth rate
#> 665 0.8091906186 1.0503799      0.157   Expected growth rate
#> 666 0.7947267742 1.0372381      0.103   Expected growth rate
#> 667 0.7202455635 0.9410478      0.004   Expected growth rate
#> 668 0.7560349029 0.9880468      0.025   Expected growth rate
#> 669 0.7976619881 1.0336279      0.122   Expected growth rate
#> 670 0.7805603339 1.0213868      0.078   Expected growth rate
#> 671 0.7593226329 0.9985797      0.036   Expected growth rate
#> 672 0.7536079362 0.9716610      0.010   Expected growth rate
#> 673 0.7332397855 0.9492793      0.004   Expected growth rate
#> 674 0.7285208220 0.9536920      0.005   Expected growth rate
#> 675 0.7532365939 0.9923791      0.028   Expected growth rate
#> 676 0.7434910327 0.9729378      0.011   Expected growth rate
#> 677 0.7268953122 0.9576080      0.007   Expected growth rate
#> 678 0.8003470505 1.0358372      0.111   Expected growth rate
#> 679 0.7562091633 0.9778812      0.019   Expected growth rate
#> 680 0.8047078991 1.0458430      0.137   Expected growth rate
#> 681 0.7887268105 1.0138535      0.065   Expected growth rate
#> 682 0.7609647551 0.9939631      0.029   Expected growth rate
#> 683 0.7465176605 0.9742778      0.013   Expected growth rate
#> 684 0.7263821177 0.9591822      0.008   Expected growth rate
#> 685 0.7202347020 0.9494289      0.002   Expected growth rate
#> 686 0.7452568958 0.9882367      0.023   Expected growth rate
#> 687 0.7362986491 0.9653131      0.011   Expected growth rate
#> 688 0.7259455633 0.9526612      0.005   Expected growth rate
#> 689 0.7326576019 0.9587883      0.006   Expected growth rate
#> 690 0.7513774153 0.9838269      0.018   Expected growth rate
#> 691 0.7552586673 0.9918062      0.027   Expected growth rate
#> 692 0.7779107349 1.0112597      0.051   Expected growth rate
#> 693 0.7673820736 0.9918635      0.026   Expected growth rate
#> 694 0.7404464371 0.9832994      0.018   Expected growth rate
#> 695 0.7285964732 0.9621522      0.003   Expected growth rate
#> 696 0.7429513533 0.9685295      0.009   Expected growth rate
#> 697 0.7555692955 0.9807608      0.018   Expected growth rate
#> 698 0.7404328139 0.9618269      0.010   Expected growth rate
#> 699 0.7315312147 0.9467174      0.003   Expected growth rate
#> 700 0.7273139339 0.9513384      0.004   Expected growth rate
#> 701 0.7508455276 0.9730833      0.010   Expected growth rate
#> 702 0.7352199451 0.9530471      0.007   Expected growth rate
#> 703 0.7214264750 0.9475747      0.005   Expected growth rate
#> 704 0.7151784199 0.9579284      0.006   Expected growth rate
#> 705 0.7422216436 0.9732051      0.014   Expected growth rate
#> 706 0.7380614412 0.9614774      0.002   Expected growth rate
#> 707 0.7222541600 0.9517494      0.005   Expected growth rate
#> 708 0.0677593595 0.6861724      0.000            Recruitment
#> 709 0.0516958321 0.6102766      0.000            Recruitment
#> 710 0.0749176347 0.7098579      0.000            Recruitment
#> 711 0.0803043180 0.7047682      0.000            Recruitment
#> 712 0.0525812888 0.6526783      0.000            Recruitment
#> 713 0.0368213413 0.5025705      0.000            Recruitment
#> 714 0.0745773898 0.7228905      0.000            Recruitment
#> 715 0.0636548715 0.6862570      0.000            Recruitment
#> 716 0.0494986000 0.5874886      0.000            Recruitment
#> 717 0.0364682861 0.5663650      0.000            Recruitment
#> 718 0.0863982142 0.7506987      0.000            Recruitment
#> 719 0.0660247682 0.6463025      0.000            Recruitment
#> 720 0.0637569830 0.6949432      0.000            Recruitment
#> 721 0.0615766183 0.6766379      0.000            Recruitment
#> 722 0.0737286161 0.6737049      0.000            Recruitment
#> 723 0.0429832045 0.6236084      0.000            Recruitment
#> 724 0.0312067783 0.4809133      0.000            Recruitment
#> 725 0.0674112996 0.7093472      0.000            Recruitment
#> 726 0.0613330905 0.7030845      0.000            Recruitment
#> 727 0.0460534251 0.5863387      0.000            Recruitment
#> 728 0.0309399465 0.5095525      0.000            Recruitment
#> 729 0.0198320066 0.4282276      0.000            Recruitment
#> 730 0.0512992256 0.6471616      0.000            Recruitment
#> 731 0.0219584109 0.4717331      0.000            Recruitment
#> 732 0.0552688348 0.6602103      0.000            Recruitment
#> 733 0.0582730538 0.6519322      0.000            Recruitment
#> 734 0.0445625768 0.6313337      0.000            Recruitment
#> 735 0.0347725920 0.4818130      0.000            Recruitment
#> 736 0.0224653544 0.4556622      0.000            Recruitment
#> 737 0.0519363957 0.6812036      0.000            Recruitment
#> 738 0.0414663209 0.5753804      0.000            Recruitment
#> 739 0.0285017104 0.4729505      0.000            Recruitment
#> 740 0.0216328055 0.4166745      0.000            Recruitment
#> 741 0.0134059469 0.3452884      0.000            Recruitment
#> 742 0.0294950638 0.5228833      0.000            Recruitment
#> 743 0.0721088351 0.7332438      0.000            Recruitment
#> 744 0.0688876572 0.7045213      0.000            Recruitment
#> 745 0.0464366382 0.5930672      0.000            Recruitment
#> 746 0.0294856747 0.5361190      0.000            Recruitment
#> 747 0.0224735494 0.4722821      0.000            Recruitment
#> 748 0.0151033981 0.4035032      0.000            Recruitment
#> 749 0.0446128061 0.5929080      0.000            Recruitment
#> 750 0.0263551402 0.4874882      0.000            Recruitment
#> 751 0.0206575252 0.4450235      0.000            Recruitment
#> 752 0.0151881758 0.3633917      0.000            Recruitment
#> 753 0.0365458294 0.5687900      0.000            Recruitment
#> 754 0.0243116284 0.4743857      0.000            Recruitment
#> 755 0.0189518493 0.4339629      0.000            Recruitment
#> 756 0.0110007883 0.3500436      0.000            Recruitment
#> 757 0.0072929689 0.3352427      0.000            Recruitment
#> 758 0.0032207791 0.2824998      0.000            Recruitment
#> 759 0.0176690172 0.3522043      0.000            Recruitment
#> 760 0.0436421579 0.5424149      0.000            Recruitment
#> 761 0.0290684480 0.4922916      0.000            Recruitment
#> 762 0.0209837968 0.4407319      0.000            Recruitment
#> 763 0.0110931252 0.3664903      0.000            Recruitment
#> 764 0.0415987222 0.5434353      0.000            Recruitment
#> 765 0.0218989960 0.4625198      0.000            Recruitment
#> 766 0.0178206990 0.4304210      0.000            Recruitment
#> 767 0.0096318370 0.3273381      0.000            Recruitment
#> 768 0.0075433391 0.3097676      0.000            Recruitment
#> 769 0.0043530896 0.2650730      0.000            Recruitment
#> 770 0.0152933431 0.4122295      0.000            Recruitment
#> 771 0.0414044713 0.5847711      0.000            Recruitment
#> 772 0.0273089671 0.4944695      0.000            Recruitment
#> 773 0.0200239703 0.4174865      0.000            Recruitment
#> 774 0.0122251871 0.3445686      0.000            Recruitment
#> 775 0.0345617632 0.5619531      0.000            Recruitment
#> 776 0.0833889111 0.7357732      0.000            Recruitment
#> 777 0.0175953101 0.4261238      0.000            Recruitment
#> 778 0.0105331102 0.3601540      0.000            Recruitment
#> 779 0.0062874561 0.3102042      0.000            Recruitment
#> 780 0.0034619174 0.2952650      0.000            Recruitment
#> 781 0.0053923523 0.2852903      0.000            Recruitment
#> 782 0.0377143248 0.5551876      0.000            Recruitment
#> 783 0.0229349701 0.4786557      0.000            Recruitment
#> 784 0.0164757634 0.4148295      0.000            Recruitment
#> 785 0.0125419399 0.3424755      0.000            Recruitment
#> 786 0.0079321016 0.3448864      0.000            Recruitment
#> 787 0.0052128941 0.2828066      0.000            Recruitment
#> 788 0.0033865380 0.2591460      0.000            Recruitment
#> 789 0.0107883127 0.3535880      0.000            Recruitment
#> 790 0.0057316983 0.3061697      0.000            Recruitment
#> 791 0.0037582246 0.2946715      0.000            Recruitment
#> 792 0.0050459818 0.3142413      0.000            Recruitment
#> 793 0.0097405703 0.3126119      0.000            Recruitment
#> 794 0.0101885214 0.2971523      0.000            Recruitment
#> 795 0.0019250277 0.2597499      0.000            Recruitment
#> 796 0.0104819330 0.3221170      0.000            Recruitment
#> 797 0.0042701477 0.2765936      0.000            Recruitment
#> 798 0.0033061280 0.2557850      0.000            Recruitment
#> 799 0.0061936690 0.3107566      0.000            Recruitment
#> 800 0.0081262223 0.3500280      0.000            Recruitment
#> 801 0.0048401417 0.3286559      0.000            Recruitment
#> 802 0.0023814095 0.2396104      0.000            Recruitment
#> 803 0.0043865611 0.2688882      0.000            Recruitment
#> 804 0.0086913722 0.3523015      0.000            Recruitment
#> 805 0.0034336087 0.2677186      0.000            Recruitment
#> 806 0.0027666343 0.2751622      0.000            Recruitment
#> 807 0.0020495866 0.2510730      0.000            Recruitment
#> 808 0.0079309644 0.3191519      0.000            Recruitment
#> 809 0.6887388994 0.9964980      0.050  Adult female survival
#> 810 0.6557443421 0.9704901      0.013  Adult female survival
#> 811 0.6619151620 0.9771581      0.006  Adult female survival
#> 812 0.6910394011 0.9979414      0.062  Adult female survival
#> 813 0.6919377124 0.9991677      0.058  Adult female survival
#> 814 0.6860739344 0.9958181      0.044  Adult female survival
#> 815 0.6838846397 0.9955492      0.048  Adult female survival
#> 816 0.6901248306 0.9949489      0.043  Adult female survival
#> 817 0.6807512454 0.9959683      0.049  Adult female survival
#> 818 0.6861139109 0.9980030      0.050  Adult female survival
#> 819 0.6561359292 0.9796098      0.011  Adult female survival
#> 820 0.6712425319 0.9852403      0.020  Adult female survival
#> 821 0.6793095193 0.9841360      0.012  Adult female survival
#> 822 0.6739463777 0.9858376      0.015  Adult female survival
#> 823 0.6931096307 0.9971782      0.057  Adult female survival
#> 824 0.6750279788 0.9857396      0.021  Adult female survival
#> 825 0.7028072994 0.9970593      0.049  Adult female survival
#> 826 0.6680598167 0.9790395      0.010  Adult female survival
#> 827 0.6709802462 0.9896601      0.024  Adult female survival
#> 828 0.6921609785 0.9969374      0.064  Adult female survival
#> 829 0.6893432635 0.9947779      0.035  Adult female survival
#> 830 0.6663882189 0.9828828      0.014  Adult female survival
#> 831 0.6701055587 0.9847051      0.014  Adult female survival
#> 832 0.6652229976 0.9864214      0.020  Adult female survival
#> 833 0.6741077451 0.9962609      0.055  Adult female survival
#> 834 0.6569072402 0.9702591      0.007  Adult female survival
#> 835 0.6785932830 0.9848583      0.012  Adult female survival
#> 836 0.6640815157 0.9825903      0.015  Adult female survival
#> 837 0.6767942893 0.9902445      0.028  Adult female survival
#> 838 0.6548840411 0.9725043      0.007  Adult female survival
#> 839 0.6733235465 0.9863100      0.020  Adult female survival
#> 840 0.6945356264 0.9948746      0.039  Adult female survival
#> 841 0.6646140253 0.9811703      0.015  Adult female survival
#> 842 0.6768942658 0.9863872      0.022  Adult female survival
#> 843 0.6595527619 0.9828196      0.014  Adult female survival
#> 844 0.7031219980 0.9935430      0.031  Adult female survival
#> 845 0.6588869937 0.9723476      0.009  Adult female survival
#> 846 0.6684965363 0.9822561      0.013  Adult female survival
#> 847 0.6633060194 0.9774225      0.008  Adult female survival
#> 848 0.6750264487 0.9947772      0.039  Adult female survival
#> 849 0.6480571225 0.9643577      0.002  Adult female survival
#> 850 0.6650491320 0.9925149      0.032  Adult female survival
#> 851 0.6834613585 0.9945560      0.035  Adult female survival
#> 852 0.6538594327 0.9755963      0.006  Adult female survival
#> 853 0.6718791593 0.9853216      0.018  Adult female survival
#> 854 0.6776467710 0.9872884      0.024  Adult female survival
#> 855 0.6831560613 0.9906978      0.027  Adult female survival
#> 856 0.6595139309 0.9723967      0.008  Adult female survival
#> 857 0.6673040106 0.9860840      0.018  Adult female survival
#> 858 0.6650532561 0.9826286      0.017  Adult female survival
#> 859 0.6819907010 0.9937219      0.037  Adult female survival
#> 860 0.6488395841 0.9676752      0.007  Adult female survival
#> 861 0.6479607028 0.9610719      0.004  Adult female survival
#> 862 0.6791781489 0.9946171      0.036  Adult female survival
#> 863 0.6625521662 0.9768561      0.009  Adult female survival
#> 864 0.6782487044 0.9870295      0.020  Adult female survival
#> 865 0.6689398981 0.9861052      0.017  Adult female survival
#> 866 0.6804900589 0.9960237      0.043  Adult female survival
#> 867 0.6469954370 0.9768115      0.007  Adult female survival
#> 868 0.6673131435 0.9867086      0.019  Adult female survival
#> 869 0.6668415042 0.9818637      0.019  Adult female survival
#> 870 0.6733212150 0.9943382      0.034  Adult female survival
#> 871 0.6562534491 0.9600548      0.001  Adult female survival
#> 872 0.6395554227 0.9712876      0.006  Adult female survival
#> 873 0.6665584029 0.9781594      0.014  Adult female survival
#> 874 0.6793841837 0.9933915      0.031  Adult female survival
#> 875 0.6509801245 0.9727598      0.003  Adult female survival
#> 876 0.6479199254 0.9682128      0.009  Adult female survival
#> 877 0.6522602593 0.9793868      0.014  Adult female survival
#> 878 0.6729596738 0.9870978      0.022  Adult female survival
#> 879 0.6573623991 0.9710874      0.007  Adult female survival
#> 880 0.6664614748 0.9787575      0.014  Adult female survival
#> 881 0.6621171472 0.9937011      0.033  Adult female survival
#> 882 0.6558474755 0.9711337      0.008  Adult female survival
#> 883 0.6562430924 0.9670402      0.009  Adult female survival
#> 884 0.6634280759 0.9766587      0.006  Adult female survival
#> 885 0.6775381112 0.9914822      0.030  Adult female survival
#> 886 0.6485201089 0.9651176      0.005  Adult female survival
#> 887 0.6495690174 0.9607484      0.004  Adult female survival
#> 888 0.6634375204 0.9821995      0.015  Adult female survival
#> 889 0.6675672377 0.9915870      0.030  Adult female survival
#> 890 0.6552058641 0.9656703      0.006  Adult female survival
#> 891 0.6655362349 0.9816444      0.008  Adult female survival
#> 892 0.6834258174 0.9919564      0.034  Adult female survival
#> 893 0.6616134341 0.9681242      0.010  Adult female survival
#> 894 0.6440545280 0.9611671      0.007  Adult female survival
#> 895 0.6591251504 0.9766966      0.012  Adult female survival
#> 896 0.6874171125 0.9873277      0.020  Adult female survival
#> 897 0.6593950082 0.9679416      0.007  Adult female survival
#> 898 0.6431004239 0.9693388      0.004  Adult female survival
#> 899 0.6705007173 0.9758619      0.007  Adult female survival
#> 900 0.6636034180 0.9884321      0.024  Adult female survival
#> 901 0.7000989641 0.9974106      0.058  Adult female survival
#> 902 0.6607501292 0.9788118      0.010  Adult female survival
#> 903 0.6812448691 0.9923942      0.029  Adult female survival
#> 904 0.6530051631 0.9709060      0.006  Adult female survival
#> 905 0.6428687368 0.9572917      0.004  Adult female survival
#> 906 0.6667997881 0.9810117      0.015  Adult female survival
#> 907 0.6791183569 0.9896423      0.024  Adult female survival
#> 908 0.6433634605 0.9577030      0.004  Adult female survival
#> 909 0.6479633509 0.9601764      0.002  Adult female survival
#> 
```

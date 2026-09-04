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
  cPars = nationalTrajectoryDefaults(),
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

  Number or vector of numbers. Initial population size for one or more
  sample populations. If NA then population growth rate is
  \$\_t=S_t\*(1+cR_t)/s\$.

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

  data frame with Anthro, Fire_excl_anthro and Year numeric columns.
  Anthro and Fire_excl_anthro are vectors of numbers between 0 and 100
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
[`addN0Variation()`](https://landscitech.github.io/caribouMetrics/dev/reference/addN0Variation.md),
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
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
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotTrajectories.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md),
[`trajectoriesFromSummaryForApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummaryForApp.md)

## Examples

``` r
trajectoriesFromNational()
#> Using saved object
#> $summary
#>     MetricTypeID PopulationName AnthroID Fire_excl_anthroID       Mean
#> 102         Rbar       National       34                  0 0.19679350
#> 103         Rbar       National       11                  0 0.29357737
#> 104         Rbar       National       19                  0 0.26009108
#> 105         Rbar       National       13                  0 0.28478098
#> 106         Rbar       National       12                  0 0.29421411
#> 107         Rbar       National       35                  0 0.19955520
#> 108         Rbar       National       23                  0 0.23800233
#> 109         Rbar       National       17                  0 0.27241829
#> 110         Rbar       National       22                  0 0.24768016
#> 111         Rbar       National       16                  0 0.27903705
#> 112         Rbar       National       10                  0 0.30370208
#> 113         Rbar       National       21                  0 0.25001602
#> 114         Rbar       National        9                  0 0.30088941
#> 115         Rbar       National       20                  0 0.25649665
#> 116         Rbar       National       14                  0 0.28586862
#> 117         Rbar       National       25                  0 0.23712377
#> 118         Rbar       National       36                  0 0.19441474
#> 119         Rbar       National       24                  0 0.24131242
#> 120         Rbar       National       18                  0 0.27050806
#> 121         Rbar       National        6                  0 0.32888429
#> 122         Rbar       National        0                  0 0.35898393
#> 123         Rbar       National       51                  0 0.15272072
#> 124         Rbar       National       62                  0 0.12870596
#> 125         Rbar       National       33                  0 0.20841433
#> 126         Rbar       National        4                  0 0.34098484
#> 127         Rbar       National       15                  0 0.28075868
#> 128         Rbar       National       26                  0 0.22859349
#> 129         Rbar       National       37                  0 0.19381262
#> 130         Rbar       National        8                  0 0.31122450
#> 131         Rbar       National       59                  0 0.13275263
#> 132         Rbar       National        7                  0 0.32340242
#> 133         Rbar       National        1                  0 0.35106220
#> 134         Rbar       National       52                  0 0.15214231
#> 135         Rbar       National       63                  0 0.12298281
#> 136         Rbar       National       74                  0 0.10462211
#> 137         Rbar       National        5                  0 0.32830395
#> 138         Rbar       National       56                  0 0.13976940
#> 139         Rbar       National       27                  0 0.22645547
#> 140         Rbar       National       38                  0 0.18829702
#> 141         Rbar       National       49                  0 0.16009776
#> 142         Rbar       National       60                  0 0.13368340
#> 143         Rbar       National       31                  0 0.20939574
#> 144         Rbar       National        2                  0 0.34397306
#> 145         Rbar       National       53                  0 0.15118493
#> 146         Rbar       National       64                  0 0.12344514
#> 147         Rbar       National       75                  0 0.10012485
#> 148         Rbar       National       46                  0 0.16289198
#> 149         Rbar       National       57                  0 0.14073510
#> 150         Rbar       National       28                  0 0.22368033
#> 151         Rbar       National       39                  0 0.19088183
#> 152         Rbar       National       50                  0 0.15559732
#> 153         Rbar       National       61                  0 0.12708447
#> 154         Rbar       National       32                  0 0.20660502
#> 155         Rbar       National        3                  0 0.34109073
#> 156         Rbar       National       54                  0 0.14370664
#> 157         Rbar       National       65                  0 0.11602402
#> 158         Rbar       National       76                  0 0.10006736
#> 159         Rbar       National       47                  0 0.16164290
#> 160         Rbar       National       58                  0 0.13974786
#> 161         Rbar       National       29                  0 0.22281638
#> 162         Rbar       National       40                  0 0.18552012
#> 163         Rbar       National       91                  0 0.07738643
#> 164         Rbar       National       45                  0 0.16970622
#> 165         Rbar       National       73                  0 0.10054610
#> 166         Rbar       National       44                  0 0.16951668
#> 167         Rbar       National       55                  0 0.14301123
#> 168         Rbar       National       66                  0 0.12221012
#> 169         Rbar       National       77                  0 0.09855718
#> 170         Rbar       National       48                  0 0.15694473
#> 171         Rbar       National       99                  0 0.06744629
#> 172         Rbar       National       30                  0 0.21547180
#> 173         Rbar       National       41                  0 0.18064203
#> 174         Rbar       National       92                  0 0.07873161
#> 175         Rbar       National       86                  0 0.08314368
#> 176         Rbar       National       97                  0 0.07173279
#> 177         Rbar       National       68                  0 0.11055399
#> 178         Rbar       National       96                  0 0.07282137
#> 179         Rbar       National       67                  0 0.11553106
#> 180         Rbar       National       78                  0 0.09731134
#> 181         Rbar       National       89                  0 0.07912311
#> 182         Rbar       National      100                  0 0.07125766
#> 183         Rbar       National       71                  0 0.10847504
#> 184         Rbar       National       42                  0 0.17573647
#> 185         Rbar       National       93                  0 0.07795523
#> 186         Rbar       National       87                  0 0.08542600
#> 187         Rbar       National       98                  0 0.06971465
#> 188         Rbar       National       69                  0 0.10950107
#> 189         Rbar       National       80                  0 0.09196845
#> 190         Rbar       National       72                  0 0.10551180
#> 191         Rbar       National       79                  0 0.09601256
#> 192         Rbar       National       90                  0 0.08242244
#> 193         Rbar       National       84                  0 0.09085084
#> 194         Rbar       National       82                  0 0.09249596
#> 195         Rbar       National       43                  0 0.17728560
#> 196         Rbar       National       94                  0 0.07392549
#> 197         Rbar       National       88                  0 0.08205178
#> 198         Rbar       National       70                  0 0.10933756
#> 199         Rbar       National       81                  0 0.08939434
#> 200         Rbar       National       83                  0 0.08899228
#> 201         Rbar       National       95                  0 0.07513466
#> 202         Rbar       National       85                  0 0.08251886
#> 203         Sbar       National       11                  0 0.86922803
#> 204         Sbar       National       23                  0 0.85774295
#> 205         Sbar       National       12                  0 0.86642043
#> 206         Sbar       National       10                  0 0.87016771
#> 207         Sbar       National        0                  0 0.87692804
#> 208         Sbar       National       15                  0 0.86550382
#> 209         Sbar       National       22                  0 0.86104807
#> 210         Sbar       National       14                  0 0.86626464
#> 211         Sbar       National       25                  0 0.86015041
#> 212         Sbar       National       13                  0 0.86709509
#> 213         Sbar       National       24                  0 0.86062866
#> 214         Sbar       National        1                  0 0.87736759
#> 215         Sbar       National       52                  0 0.83985803
#> 216         Sbar       National       63                  0 0.83252134
#> 217         Sbar       National       51                  0 0.84061528
#> 218         Sbar       National        5                  0 0.87473548
#> 219         Sbar       National       16                  0 0.86803921
#> 220         Sbar       National       27                  0 0.85788156
#> 221         Sbar       National       38                  0 0.84882519
#> 222         Sbar       National       26                  0 0.85816905
#> 223         Sbar       National       20                  0 0.86021641
#> 224         Sbar       National        8                  0 0.86805750
#> 225         Sbar       National        2                  0 0.87494729
#> 226         Sbar       National       53                  0 0.84007790
#> 227         Sbar       National       64                  0 0.83332730
#> 228         Sbar       National       35                  0 0.85104913
#> 229         Sbar       National        6                  0 0.87108224
#> 230         Sbar       National       17                  0 0.86471750
#> 231         Sbar       National       28                  0 0.85781494
#> 232         Sbar       National       39                  0 0.84817345
#> 233         Sbar       National       50                  0 0.83863882
#> 234         Sbar       National       21                  0 0.86228989
#> 235         Sbar       National        9                  0 0.87010020
#> 236         Sbar       National        3                  0 0.87589869
#> 237         Sbar       National       54                  0 0.83765190
#> 238         Sbar       National       65                  0 0.82987436
#> 239         Sbar       National       36                  0 0.84835380
#> 240         Sbar       National        7                  0 0.87063490
#> 241         Sbar       National       18                  0 0.86189931
#> 242         Sbar       National       29                  0 0.85552526
#> 243         Sbar       National       40                  0 0.84659032
#> 244         Sbar       National       91                  0 0.81435371
#> 245         Sbar       National       62                  0 0.83362138
#> 246         Sbar       National       33                  0 0.85050674
#> 247         Sbar       National        4                  0 0.87443906
#> 248         Sbar       National       55                  0 0.83715373
#> 249         Sbar       National       66                  0 0.82982449
#> 250         Sbar       National       37                  0 0.85256561
#> 251         Sbar       National       48                  0 0.84145666
#> 252         Sbar       National       19                  0 0.86436629
#> 253         Sbar       National       30                  0 0.85390640
#> 254         Sbar       National       41                  0 0.84454085
#> 255         Sbar       National       92                  0 0.81303063
#> 256         Sbar       National       46                  0 0.84391407
#> 257         Sbar       National       34                  0 0.85319380
#> 258         Sbar       National       45                  0 0.84525949
#> 259         Sbar       National       56                  0 0.83726963
#> 260         Sbar       National       67                  0 0.82754863
#> 261         Sbar       National       78                  0 0.82199986
#> 262         Sbar       National       49                  0 0.84276542
#> 263         Sbar       National       60                  0 0.83175395
#> 264         Sbar       National       31                  0 0.85550721
#> 265         Sbar       National       42                  0 0.84500553
#> 266         Sbar       National       93                  0 0.81142020
#> 267         Sbar       National       47                  0 0.84287618
#> 268         Sbar       National       75                  0 0.82494778
#> 269         Sbar       National       69                  0 0.82571582
#> 270         Sbar       National       57                  0 0.83717271
#> 271         Sbar       National       68                  0 0.82769307
#> 272         Sbar       National       79                  0 0.82109160
#> 273         Sbar       National       90                  0 0.81328640
#> 274         Sbar       National       61                  0 0.83291594
#> 275         Sbar       National       32                  0 0.85276280
#> 276         Sbar       National       43                  0 0.84718569
#> 277         Sbar       National       94                  0 0.81172672
#> 278         Sbar       National       88                  0 0.81579443
#> 279         Sbar       National       76                  0 0.82244107
#> 280         Sbar       National       70                  0 0.82682127
#> 281         Sbar       National       58                  0 0.83565925
#> 282         Sbar       National       98                  0 0.80907928
#> 283         Sbar       National       80                  0 0.82022927
#> 284         Sbar       National       74                  0 0.82461905
#> 285         Sbar       National       85                  0 0.81676191
#> 286         Sbar       National       73                  0 0.82618549
#> 287         Sbar       National       44                  0 0.84689296
#> 288         Sbar       National       95                  0 0.81118648
#> 289         Sbar       National       89                  0 0.81250274
#> 290         Sbar       National       77                  0 0.82567810
#> 291         Sbar       National       71                  0 0.83042618
#> 292         Sbar       National       59                  0 0.83311193
#> 293         Sbar       National       99                  0 0.80938199
#> 294         Sbar       National       81                  0 0.81884366
#> 295         Sbar       National       96                  0 0.81114819
#> 296         Sbar       National       86                  0 0.81500666
#> 297         Sbar       National       97                  0 0.81197776
#> 298         Sbar       National       72                  0 0.82672456
#> 299         Sbar       National      100                  0 0.80842074
#> 300         Sbar       National       83                  0 0.82116606
#> 301         Sbar       National       84                  0 0.81901860
#> 302         Sbar       National       87                  0 0.81433497
#> 303         Sbar       National       82                  0 0.82057447
#> 304            X       National       11                  0 0.14157688
#> 305            X       National       16                  0 0.13709730
#> 306            X       National       13                  0 0.13847764
#> 307            X       National       12                  0 0.14492214
#> 308            X       National       14                  0 0.14171031
#> 309            X       National       52                  0 0.07569796
#> 310            X       National        0                  0 0.16791547
#> 311            X       National       17                  0 0.13011144
#> 312            X       National        1                  0 0.16533176
#> 313            X       National       39                  0 0.09514465
#> 314            X       National       27                  0 0.10926246
#> 315            X       National       15                  0 0.13490805
#> 316            X       National       28                  0 0.11014933
#> 317            X       National        3                  0 0.16169538
#> 318            X       National       54                  0 0.07053353
#> 319            X       National        2                  0 0.16931117
#> 320            X       National       53                  0 0.07620150
#> 321            X       National       24                  0 0.11608607
#> 322            X       National       18                  0 0.12970623
#> 323            X       National       29                  0 0.11030067
#> 324            X       National       40                  0 0.09220507
#> 325            X       National       51                  0 0.07559507
#> 326            X       National       22                  0 0.12378060
#> 327            X       National       10                  0 0.14980806
#> 328            X       National        4                  0 0.16599289
#> 329            X       National       55                  0 0.07115050
#> 330            X       National       26                  0 0.11219420
#> 331            X       National       37                  0 0.09538539
#> 332            X       National       25                  0 0.11790635
#> 333            X       National       19                  0 0.13080211
#> 334            X       National       30                  0 0.10677547
#> 335            X       National       41                  0 0.08948101
#> 336            X       National       92                  0 0.03867391
#> 337            X       National       23                  0 0.11770491
#> 338            X       National       34                  0 0.09706254
#> 339            X       National        5                  0 0.15994256
#> 340            X       National       56                  0 0.07014122
#> 341            X       National       67                  0 0.05812598
#> 342            X       National       38                  0 0.09396343
#> 343            X       National        9                  0 0.14737775
#> 344            X       National       20                  0 0.12581003
#> 345            X       National       31                  0 0.10523840
#> 346            X       National       42                  0 0.08636844
#> 347            X       National       93                  0 0.04099213
#> 348            X       National       64                  0 0.05979695
#> 349            X       National       35                  0 0.09719679
#> 350            X       National        6                  0 0.16276045
#> 351            X       National       57                  0 0.07165826
#> 352            X       National       68                  0 0.05473115
#> 353            X       National       79                  0 0.04818485
#> 354            X       National       50                  0 0.07904066
#> 355            X       National       21                  0 0.12245213
#> 356            X       National       32                  0 0.10368981
#> 357            X       National       43                  0 0.08867022
#> 358            X       National       94                  0 0.03700714
#> 359            X       National       65                  0 0.05750284
#> 360            X       National       36                  0 0.09866554
#> 361            X       National        7                  0 0.15534079
#> 362            X       National       58                  0 0.06860039
#> 363            X       National       69                  0 0.05482511
#> 364            X       National       80                  0 0.04656174
#> 365            X       National       91                  0 0.03913878
#> 366            X       National       62                  0 0.06625571
#> 367            X       National       33                  0 0.10439409
#> 368            X       National       44                  0 0.08503031
#> 369            X       National       95                  0 0.03775301
#> 370            X       National       66                  0 0.06101637
#> 371            X       National       77                  0 0.04960729
#> 372            X       National        8                  0 0.15353744
#> 373            X       National       59                  0 0.06666621
#> 374            X       National       70                  0 0.05450682
#> 375            X       National       81                  0 0.04613909
#> 376            X       National       75                  0 0.05090747
#> 377            X       National       63                  0 0.06385875
#> 378            X       National       74                  0 0.05297048
#> 379            X       National       45                  0 0.08316977
#> 380            X       National       96                  0 0.03435668
#> 381            X       National       90                  0 0.04222365
#> 382            X       National       78                  0 0.04921643
#> 383            X       National       49                  0 0.08136226
#> 384            X       National       60                  0 0.06625139
#> 385            X       National       71                  0 0.05316760
#> 386            X       National       82                  0 0.04527759
#> 387            X       National       76                  0 0.04821798
#> 388            X       National       47                  0 0.08013279
#> 389            X       National       98                  0 0.03390860
#> 390            X       National       46                  0 0.08046601
#> 391            X       National       97                  0 0.03491251
#> 392            X       National       72                  0 0.05327130
#> 393            X       National       85                  0 0.04146172
#> 394            X       National       73                  0 0.05115017
#> 395            X       National       61                  0 0.06468771
#> 396            X       National       99                  0 0.03327897
#> 397            X       National       83                  0 0.04394301
#> 398            X       National      100                  0 0.03487037
#> 399            X       National       48                  0 0.07876249
#> 400            X       National       86                  0 0.04157952
#> 401            X       National       88                  0 0.04085655
#> 402            X       National       87                  0 0.04299027
#> 403            X       National       84                  0 0.04542641
#> 404            X       National       89                  0 0.03972046
#> 405            c       National       18                  0 1.00000000
#> 406            c       National       13                  0 1.00000000
#> 407            c       National       16                  0 1.00000000
#> 408            c       National       17                  0 1.00000000
#> 409            c       National        0                  0 1.00000000
#> 410            c       National       28                  0 1.00000000
#> 411            c       National        3                  0 1.00000000
#> 412            c       National       14                  0 1.00000000
#> 413            c       National        4                  0 1.00000000
#> 414            c       National       19                  0 1.00000000
#> 415            c       National        1                  0 1.00000000
#> 416            c       National       41                  0 1.00000000
#> 417            c       National       29                  0 1.00000000
#> 418            c       National       23                  0 1.00000000
#> 419            c       National       11                  0 1.00000000
#> 420            c       National        5                  0 1.00000000
#> 421            c       National       56                  0 1.00000000
#> 422            c       National       27                  0 1.00000000
#> 423            c       National       15                  0 1.00000000
#> 424            c       National       26                  0 1.00000000
#> 425            c       National       20                  0 1.00000000
#> 426            c       National        2                  0 1.00000000
#> 427            c       National       42                  0 1.00000000
#> 428            c       National       30                  0 1.00000000
#> 429            c       National       24                  0 1.00000000
#> 430            c       National       12                  0 1.00000000
#> 431            c       National        6                  0 1.00000000
#> 432            c       National       57                  0 1.00000000
#> 433            c       National       68                  0 1.00000000
#> 434            c       National       39                  0 1.00000000
#> 435            c       National       10                  0 1.00000000
#> 436            c       National       21                  0 1.00000000
#> 437            c       National       32                  0 1.00000000
#> 438            c       National       43                  0 1.00000000
#> 439            c       National       31                  0 1.00000000
#> 440            c       National       25                  0 1.00000000
#> 441            c       National       53                  0 1.00000000
#> 442            c       National        7                  0 1.00000000
#> 443            c       National       58                  0 1.00000000
#> 444            c       National       69                  0 1.00000000
#> 445            c       National       40                  0 1.00000000
#> 446            c       National       51                  0 1.00000000
#> 447            c       National       22                  0 1.00000000
#> 448            c       National       33                  0 1.00000000
#> 449            c       National       44                  0 1.00000000
#> 450            c       National       55                  0 1.00000000
#> 451            c       National       66                  0 1.00000000
#> 452            c       National       54                  0 1.00000000
#> 453            c       National        8                  0 1.00000000
#> 454            c       National       59                  0 1.00000000
#> 455            c       National       70                  0 1.00000000
#> 456            c       National       81                  0 1.00000000
#> 457            c       National       52                  0 1.00000000
#> 458            c       National       63                  0 1.00000000
#> 459            c       National       34                  0 1.00000000
#> 460            c       National       45                  0 1.00000000
#> 461            c       National       96                  0 1.00000000
#> 462            c       National       67                  0 1.00000000
#> 463            c       National       38                  0 1.00000000
#> 464            c       National        9                  0 1.00000000
#> 465            c       National       60                  0 1.00000000
#> 466            c       National       71                  0 1.00000000
#> 467            c       National       82                  0 1.00000000
#> 468            c       National       93                  0 1.00000000
#> 469            c       National       64                  0 1.00000000
#> 470            c       National       35                  0 1.00000000
#> 471            c       National       46                  0 1.00000000
#> 472            c       National       97                  0 1.00000000
#> 473            c       National       91                  0 1.00000000
#> 474            c       National       79                  0 1.00000000
#> 475            c       National       50                  0 1.00000000
#> 476            c       National       61                  0 1.00000000
#> 477            c       National       72                  0 1.00000000
#> 478            c       National       83                  0 1.00000000
#> 479            c       National       94                  0 1.00000000
#> 480            c       National       65                  0 1.00000000
#> 481            c       National       36                  0 1.00000000
#> 482            c       National       47                  0 1.00000000
#> 483            c       National       98                  0 1.00000000
#> 484            c       National       92                  0 1.00000000
#> 485            c       National       80                  0 1.00000000
#> 486            c       National       74                  0 1.00000000
#> 487            c       National       62                  0 1.00000000
#> 488            c       National       73                  0 1.00000000
#> 489            c       National       84                  0 1.00000000
#> 490            c       National       95                  0 1.00000000
#> 491            c       National       49                  0 1.00000000
#> 492            c       National       37                  0 1.00000000
#> 493            c       National       48                  0 1.00000000
#> 494            c       National       99                  0 1.00000000
#> 495            c       National       76                  0 1.00000000
#> 496            c       National       87                  0 1.00000000
#> 497            c       National       75                  0 1.00000000
#> 498            c       National       86                  0 1.00000000
#> 499            c       National       78                  0 1.00000000
#> 500            c       National       85                  0 1.00000000
#> 501            c       National      100                  0 1.00000000
#> 502            c       National       90                  0 1.00000000
#> 503            c       National       88                  0 1.00000000
#> 504            c       National       77                  0 1.00000000
#> 505            c       National       89                  0 1.00000000
#> 506       lambda       National       15                  0 0.98577000
#> 507       lambda       National        5                  0 1.01319000
#> 508       lambda       National       17                  0 0.98067400
#> 509       lambda       National       19                  0 0.97901000
#> 510       lambda       National       30                  0 0.94867000
#> 511       lambda       National       12                  0 0.99225300
#> 512       lambda       National        6                  0 1.01313500
#> 513       lambda       National       57                  0 0.89719800
#> 514       lambda       National       18                  0 0.97245900
#> 515       lambda       National       16                  0 0.98295000
#> 516       lambda       National       10                  0 1.00225700
#> 517       lambda       National       21                  0 0.96890100
#> 518       lambda       National       28                  0 0.95847000
#> 519       lambda       National       20                  0 0.97057200
#> 520       lambda       National       31                  0 0.94942200
#> 521       lambda       National        2                  0 1.02407700
#> 522       lambda       National       13                  0 0.98895600
#> 523       lambda       National        7                  0 1.00910500
#> 524       lambda       National       58                  0 0.89781200
#> 525       lambda       National       29                  0 0.95401500
#> 526       lambda       National        0                  0 1.02239200
#> 527       lambda       National       11                  0 0.99456100
#> 528       lambda       National       22                  0 0.96872200
#> 529       lambda       National       33                  0 0.94235700
#> 530       lambda       National        4                  0 1.02455400
#> 531       lambda       National       32                  0 0.94720300
#> 532       lambda       National        3                  0 1.01667100
#> 533       lambda       National       14                  0 0.99314700
#> 534       lambda       National        8                  0 1.01002400
#> 535       lambda       National       59                  0 0.88944200
#> 536       lambda       National       70                  0 0.87312400
#> 537       lambda       National        1                  0 1.02751800
#> 538       lambda       National       52                  0 0.90320500
#> 539       lambda       National       23                  0 0.96223300
#> 540       lambda       National       34                  0 0.94069100
#> 541       lambda       National       45                  0 0.91587800
#> 542       lambda       National       56                  0 0.89559700
#> 543       lambda       National       27                  0 0.95505700
#> 544       lambda       National       55                  0 0.89728700
#> 545       lambda       National        9                  0 0.99947100
#> 546       lambda       National       60                  0 0.88981500
#> 547       lambda       National       71                  0 0.87620200
#> 548       lambda       National       42                  0 0.91951000
#> 549       lambda       National       53                  0 0.90658900
#> 550       lambda       National       24                  0 0.96216800
#> 551       lambda       National       35                  0 0.93225300
#> 552       lambda       National       46                  0 0.91208400
#> 553       lambda       National       97                  0 0.83773600
#> 554       lambda       National       68                  0 0.87581600
#> 555       lambda       National       39                  0 0.93030700
#> 556       lambda       National       50                  0 0.90952000
#> 557       lambda       National       61                  0 0.89154400
#> 558       lambda       National       72                  0 0.87179000
#> 559       lambda       National       43                  0 0.92291600
#> 560       lambda       National       54                  0 0.89971800
#> 561       lambda       National       25                  0 0.96761200
#> 562       lambda       National       36                  0 0.93689200
#> 563       lambda       National       47                  0 0.91374600
#> 564       lambda       National       98                  0 0.83821900
#> 565       lambda       National       69                  0 0.87637000
#> 566       lambda       National       40                  0 0.92299600
#> 567       lambda       National       51                  0 0.90745900
#> 568       lambda       National       62                  0 0.89039400
#> 569       lambda       National       73                  0 0.86709200
#> 570       lambda       National       44                  0 0.91942600
#> 571       lambda       National       95                  0 0.84464400
#> 572       lambda       National       26                  0 0.95570300
#> 573       lambda       National       37                  0 0.93476600
#> 574       lambda       National       48                  0 0.91595800
#> 575       lambda       National       99                  0 0.83979800
#> 576       lambda       National       93                  0 0.84514600
#> 577       lambda       National       41                  0 0.92107600
#> 578       lambda       National       92                  0 0.84332000
#> 579       lambda       National       63                  0 0.88690600
#> 580       lambda       National       74                  0 0.86812400
#> 581       lambda       National       85                  0 0.84657100
#> 582       lambda       National       96                  0 0.84184800
#> 583       lambda       National       67                  0 0.87960400
#> 584       lambda       National       38                  0 0.92936300
#> 585       lambda       National       49                  0 0.90598700
#> 586       lambda       National      100                  0 0.83843300
#> 587       lambda       National       94                  0 0.84139500
#> 588       lambda       National       82                  0 0.86253100
#> 589       lambda       National       76                  0 0.86400800
#> 590       lambda       National       64                  0 0.88713700
#> 591       lambda       National       75                  0 0.86738500
#> 592       lambda       National       86                  0 0.85080500
#> 593       lambda       National       80                  0 0.86075500
#> 594       lambda       National       91                  0 0.84960500
#> 595       lambda       National       79                  0 0.86233700
#> 596       lambda       National       90                  0 0.85156900
#> 597       lambda       National       84                  0 0.85738300
#> 598       lambda       National       78                  0 0.86611500
#> 599       lambda       National       83                  0 0.85952200
#> 600       lambda       National       77                  0 0.86945400
#> 601       lambda       National       65                  0 0.88044300
#> 602       lambda       National       88                  0 0.85443200
#> 603       lambda       National       87                  0 0.85230600
#> 604       lambda       National       81                  0 0.85819600
#> 605       lambda       National       89                  0 0.84593000
#> 606       lambda       National       66                  0 0.88428600
#> 607   lambda_bar       National       19                  0 0.97672419
#> 608   lambda_bar       National       17                  0 0.98259360
#> 609   lambda_bar       National       18                  0 0.97848062
#> 610   lambda_bar       National        6                  0 1.01441863
#> 611   lambda_bar       National       21                  0 0.97004368
#> 612   lambda_bar       National       11                  0 0.99691637
#> 613   lambda_bar       National       22                  0 0.96765475
#> 614   lambda_bar       National       10                  0 1.00236643
#> 615   lambda_bar       National        8                  0 1.00313729
#> 616   lambda_bar       National       59                  0 0.88848467
#> 617   lambda_bar       National       20                  0 0.97049699
#> 618   lambda_bar       National        1                  0 1.03135777
#> 619   lambda_bar       National       12                  0 0.99397791
#> 620   lambda_bar       National       23                  0 0.95983921
#> 621   lambda_bar       National        7                  0 1.01134006
#> 622   lambda_bar       National        5                  0 1.01841500
#> 623   lambda_bar       National       16                  0 0.98925986
#> 624   lambda_bar       National        4                  0 1.02340467
#> 625   lambda_bar       National       34                  0 0.93711479
#> 626   lambda_bar       National        9                  0 1.00098165
#> 627   lambda_bar       National       60                  0 0.88732244
#> 628   lambda_bar       National       31                  0 0.94503276
#> 629   lambda_bar       National        2                  0 1.02531148
#> 630   lambda_bar       National       13                  0 0.99057333
#> 631   lambda_bar       National       24                  0 0.96439307
#> 632   lambda_bar       National       35                  0 0.93584138
#> 633   lambda_bar       National       46                  0 0.91267267
#> 634   lambda_bar       National       57                  0 0.89618663
#> 635   lambda_bar       National       28                  0 0.95375732
#> 636   lambda_bar       National       39                  0 0.92907188
#> 637   lambda_bar       National       50                  0 0.90387306
#> 638   lambda_bar       National       61                  0 0.88576560
#> 639   lambda_bar       National       32                  0 0.94072728
#> 640   lambda_bar       National        3                  0 1.02529206
#> 641   lambda_bar       National       14                  0 0.99002421
#> 642   lambda_bar       National       25                  0 0.96213939
#> 643   lambda_bar       National       36                  0 0.93088570
#> 644   lambda_bar       National       47                  0 0.91102883
#> 645   lambda_bar       National       58                  0 0.89408254
#> 646   lambda_bar       National       29                  0 0.95088218
#> 647   lambda_bar       National       40                  0 0.92515549
#> 648   lambda_bar       National       51                  0 0.90478181
#> 649   lambda_bar       National       62                  0 0.88725551
#> 650   lambda_bar       National       33                  0 0.93901688
#> 651   lambda_bar       National        0                  0 1.03429719
#> 652   lambda_bar       National       15                  0 0.98704517
#> 653   lambda_bar       National       26                  0 0.95625417
#> 654   lambda_bar       National       37                  0 0.93511172
#> 655   lambda_bar       National       48                  0 0.90759528
#> 656   lambda_bar       National       99                  0 0.83667308
#> 657   lambda_bar       National       30                  0 0.94596865
#> 658   lambda_bar       National       41                  0 0.92082207
#> 659   lambda_bar       National       52                  0 0.90384962
#> 660   lambda_bar       National       63                  0 0.88377912
#> 661   lambda_bar       National       74                  0 0.86783282
#> 662   lambda_bar       National       45                  0 0.91699227
#> 663   lambda_bar       National       56                  0 0.89573747
#> 664   lambda_bar       National       27                  0 0.95506261
#> 665   lambda_bar       National       38                  0 0.92878566
#> 666   lambda_bar       National       49                  0 0.91038549
#> 667   lambda_bar       National      100                  0 0.83721228
#> 668   lambda_bar       National       71                  0 0.87549791
#> 669   lambda_bar       National       42                  0 0.91924917
#> 670   lambda_bar       National       53                  0 0.90365639
#> 671   lambda_bar       National       64                  0 0.88472754
#> 672   lambda_bar       National       75                  0 0.86619951
#> 673   lambda_bar       National       86                  0 0.84888595
#> 674   lambda_bar       National       97                  0 0.84110749
#> 675   lambda_bar       National       68                  0 0.87350959
#> 676   lambda_bar       National       79                  0 0.86058682
#> 677   lambda_bar       National       90                  0 0.84669958
#> 678   lambda_bar       National       44                  0 0.91861128
#> 679   lambda_bar       National       72                  0 0.87031933
#> 680   lambda_bar       National       43                  0 0.92227859
#> 681   lambda_bar       National       54                  0 0.89789329
#> 682   lambda_bar       National       65                  0 0.87807409
#> 683   lambda_bar       National       76                  0 0.86362167
#> 684   lambda_bar       National       87                  0 0.84908228
#> 685   lambda_bar       National       98                  0 0.83726629
#> 686   lambda_bar       National       69                  0 0.87089586
#> 687   lambda_bar       National       80                  0 0.85791720
#> 688   lambda_bar       National       91                  0 0.84590009
#> 689   lambda_bar       National       85                  0 0.85039637
#> 690   lambda_bar       National       73                  0 0.86764527
#> 691   lambda_bar       National       67                  0 0.87529042
#> 692   lambda_bar       National       55                  0 0.89703529
#> 693   lambda_bar       National       66                  0 0.88048253
#> 694   lambda_bar       National       77                  0 0.86630473
#> 695   lambda_bar       National       88                  0 0.84920977
#> 696   lambda_bar       National       82                  0 0.85850427
#> 697   lambda_bar       National       70                  0 0.87205032
#> 698   lambda_bar       National       81                  0 0.85546406
#> 699   lambda_bar       National       92                  0 0.84502077
#> 700   lambda_bar       National       94                  0 0.84174146
#> 701   lambda_bar       National       78                  0 0.86204796
#> 702   lambda_bar       National       89                  0 0.84461589
#> 703   lambda_bar       National       96                  0 0.84070183
#> 704   lambda_bar       National       93                  0 0.84309810
#> 705   lambda_bar       National       83                  0 0.85772507
#> 706   lambda_bar       National       84                  0 0.85629850
#> 707   lambda_bar       National       95                  0 0.84165671
#> 708  recruitment       National        8                  0 0.30707488
#> 709  recruitment       National       19                  0 0.26160422
#> 710  recruitment       National        6                  0 0.32552090
#> 711  recruitment       National        4                  0 0.33198577
#> 712  recruitment       National       20                  0 0.25162006
#> 713  recruitment       National       35                  0 0.19439358
#> 714  recruitment       National        7                  0 0.31068158
#> 715  recruitment       National        9                  0 0.29475551
#> 716  recruitment       National       24                  0 0.23217213
#> 717  recruitment       National       39                  0 0.19028931
#> 718  recruitment       National        2                  0 0.33862234
#> 719  recruitment       National       13                  0 0.27695528
#> 720  recruitment       National        5                  0 0.31988513
#> 721  recruitment       National       16                  0 0.27419460
#> 722  recruitment       National       10                  0 0.29961611
#> 723  recruitment       National       21                  0 0.24490426
#> 724  recruitment       National       36                  0 0.19733109
#> 725  recruitment       National        3                  0 0.32339075
#> 726  recruitment       National       14                  0 0.28342062
#> 727  recruitment       National       25                  0 0.23581271
#> 728  recruitment       National       40                  0 0.18441013
#> 729  recruitment       National       51                  0 0.15119015
#> 730  recruitment       National       18                  0 0.25941246
#> 731  recruitment       National       50                  0 0.15808133
#> 732  recruitment       National       17                  0 0.26022288
#> 733  recruitment       National       11                  0 0.28315377
#> 734  recruitment       National       22                  0 0.24756119
#> 735  recruitment       National       37                  0 0.19077078
#> 736  recruitment       National       48                  0 0.15752499
#> 737  recruitment       National       15                  0 0.26981610
#> 738  recruitment       National       26                  0 0.22438840
#> 739  recruitment       National       41                  0 0.17896202
#> 740  recruitment       National       52                  0 0.15139592
#> 741  recruitment       National       63                  0 0.12771750
#> 742  recruitment       National       34                  0 0.19412508
#> 743  recruitment       National        1                  0 0.33066352
#> 744  recruitment       National       12                  0 0.28984428
#> 745  recruitment       National       23                  0 0.23540981
#> 746  recruitment       National       38                  0 0.18792686
#> 747  recruitment       National       49                  0 0.16272451
#> 748  recruitment       National       60                  0 0.13250279
#> 749  recruitment       National       27                  0 0.21852491
#> 750  recruitment       National       42                  0 0.17273688
#> 751  recruitment       National       53                  0 0.15240301
#> 752  recruitment       National       64                  0 0.11959389
#> 753  recruitment       National       31                  0 0.21047679
#> 754  recruitment       National       46                  0 0.16093202
#> 755  recruitment       National       57                  0 0.14331653
#> 756  recruitment       National       68                  0 0.10946230
#> 757  recruitment       National       79                  0 0.09636969
#> 758  recruitment       National       90                  0 0.08444729
#> 759  recruitment       National       61                  0 0.12937541
#> 760  recruitment       National       28                  0 0.22029865
#> 761  recruitment       National       43                  0 0.17734045
#> 762  recruitment       National       54                  0 0.14106706
#> 763  recruitment       National       65                  0 0.11500568
#> 764  recruitment       National       32                  0 0.20737962
#> 765  recruitment       National       47                  0 0.16026557
#> 766  recruitment       National       58                  0 0.13720078
#> 767  recruitment       National       69                  0 0.10965021
#> 768  recruitment       National       80                  0 0.09312348
#> 769  recruitment       National       91                  0 0.07827757
#> 770  recruitment       National       62                  0 0.13251142
#> 771  recruitment       National       29                  0 0.22060134
#> 772  recruitment       National       44                  0 0.17006062
#> 773  recruitment       National       55                  0 0.14230100
#> 774  recruitment       National       66                  0 0.12203273
#> 775  recruitment       National       33                  0 0.20878819
#> 776  recruitment       National        0                  0 0.33583094
#> 777  recruitment       National       59                  0 0.13333241
#> 778  recruitment       National       70                  0 0.10901364
#> 779  recruitment       National       81                  0 0.09227818
#> 780  recruitment       National       92                  0 0.07734781
#> 781  recruitment       National       86                  0 0.08315903
#> 782  recruitment       National       30                  0 0.21355094
#> 783  recruitment       National       45                  0 0.16633954
#> 784  recruitment       National       56                  0 0.14028244
#> 785  recruitment       National       67                  0 0.11625196
#> 786  recruitment       National       78                  0 0.09843286
#> 787  recruitment       National       89                  0 0.07944093
#> 788  recruitment       National      100                  0 0.06974074
#> 789  recruitment       National       71                  0 0.10633520
#> 790  recruitment       National       82                  0 0.09055517
#> 791  recruitment       National       93                  0 0.08198426
#> 792  recruitment       National       87                  0 0.08598054
#> 793  recruitment       National       75                  0 0.10181494
#> 794  recruitment       National       73                  0 0.10230034
#> 795  recruitment       National       97                  0 0.06982502
#> 796  recruitment       National       74                  0 0.10594096
#> 797  recruitment       National       85                  0 0.08292345
#> 798  recruitment       National       94                  0 0.07401427
#> 799  recruitment       National       84                  0 0.09085283
#> 800  recruitment       National       72                  0 0.10654259
#> 801  recruitment       National       83                  0 0.08788602
#> 802  recruitment       National       98                  0 0.06781720
#> 803  recruitment       National       88                  0 0.08171310
#> 804  recruitment       National       76                  0 0.09643597
#> 805  recruitment       National       95                  0 0.07550602
#> 806  recruitment       National       96                  0 0.06871335
#> 807  recruitment       National       99                  0 0.06655794
#> 808  recruitment       National       77                  0 0.09921457
#> 809     survival       National        4                  0 0.88001696
#> 810     survival       National       69                  0 0.82968943
#> 811     survival       National       68                  0 0.83024213
#> 812     survival       National        2                  0 0.87673857
#> 813     survival       National        5                  0 0.87338693
#> 814     survival       National        7                  0 0.87367488
#> 815     survival       National       10                  0 0.87128856
#> 816     survival       National        6                  0 0.87172706
#> 817     survival       National        9                  0 0.87109774
#> 818     survival       National       11                  0 0.87143199
#> 819     survival       National       70                  0 0.82743933
#> 820     survival       National       37                  0 0.85277519
#> 821     survival       National       48                  0 0.84916888
#> 822     survival       National       36                  0 0.85262298
#> 823     survival       National        3                  0 0.87489706
#> 824     survival       National       35                  0 0.84932293
#> 825     survival       National        8                  0 0.87679430
#> 826     survival       National       67                  0 0.83150570
#> 827     survival       National       34                  0 0.85768093
#> 828     survival       National        1                  0 0.88249059
#> 829     survival       National       12                  0 0.86724988
#> 830     survival       National       71                  0 0.83260064
#> 831     survival       National       38                  0 0.84996450
#> 832     survival       National       49                  0 0.83812708
#> 833     survival       National       16                  0 0.86436934
#> 834     survival       National       75                  0 0.82582264
#> 835     survival       National       42                  0 0.84614665
#> 836     survival       National       53                  0 0.84202560
#> 837     survival       National       20                  0 0.86244687
#> 838     survival       National       79                  0 0.82186844
#> 839     survival       National       46                  0 0.84452861
#> 840     survival       National       13                  0 0.86915373
#> 841     survival       National       72                  0 0.82752235
#> 842     survival       National       39                  0 0.84948051
#> 843     survival       National       50                  0 0.84258051
#> 844     survival       National       17                  0 0.86779674
#> 845     survival       National       76                  0 0.82408851
#> 846     survival       National       43                  0 0.84781604
#> 847     survival       National       54                  0 0.84026026
#> 848     survival       National       21                  0 0.86302362
#> 849     survival       National       80                  0 0.82268664
#> 850     survival       National       47                  0 0.84601772
#> 851     survival       National       14                  0 0.86996979
#> 852     survival       National       73                  0 0.82489854
#> 853     survival       National       40                  0 0.84438588
#> 854     survival       National       51                  0 0.84368429
#> 855     survival       National       18                  0 0.86072651
#> 856     survival       National       77                  0 0.82840644
#> 857     survival       National       44                  0 0.84707654
#> 858     survival       National       55                  0 0.83766791
#> 859     survival       National       22                  0 0.86215241
#> 860     survival       National       81                  0 0.82054516
#> 861     survival       National       92                  0 0.81270929
#> 862     survival       National       15                  0 0.86975988
#> 863     survival       National       74                  0 0.82493585
#> 864     survival       National       41                  0 0.84530085
#> 865     survival       National       52                  0 0.84019130
#> 866     survival       National       19                  0 0.86650111
#> 867     survival       National       78                  0 0.82484722
#> 868     survival       National       45                  0 0.84533241
#> 869     survival       National       56                  0 0.83688880
#> 870     survival       National       23                  0 0.86126152
#> 871     survival       National       82                  0 0.82487467
#> 872     survival       National       93                  0 0.81131617
#> 873     survival       National       60                  0 0.83444448
#> 874     survival       National       27                  0 0.86155165
#> 875     survival       National       86                  0 0.81606726
#> 876     survival       National       97                  0 0.80988486
#> 877     survival       National       64                  0 0.83689506
#> 878     survival       National       31                  0 0.85909774
#> 879     survival       National       90                  0 0.81679743
#> 880     survival       National       57                  0 0.83760983
#> 881     survival       National       24                  0 0.86192291
#> 882     survival       National       83                  0 0.82309669
#> 883     survival       National       94                  0 0.81088966
#> 884     survival       National       61                  0 0.83706469
#> 885     survival       National       28                  0 0.86304184
#> 886     survival       National       87                  0 0.81691221
#> 887     survival       National       98                  0 0.81041529
#> 888     survival       National       65                  0 0.83263546
#> 889     survival       National       32                  0 0.85773020
#> 890     survival       National       91                  0 0.81751526
#> 891     survival       National       58                  0 0.83938139
#> 892     survival       National       25                  0 0.86589699
#> 893     survival       National       84                  0 0.82087391
#> 894     survival       National       95                  0 0.81418362
#> 895     survival       National       62                  0 0.83540223
#> 896     survival       National       29                  0 0.85915522
#> 897     survival       National       88                  0 0.82077985
#> 898     survival       National       99                  0 0.81249195
#> 899     survival       National       66                  0 0.83359695
#> 900     survival       National       33                  0 0.85382491
#> 901     survival       National        0                  0 0.87610150
#> 902     survival       National       59                  0 0.83425808
#> 903     survival       National       26                  0 0.86041257
#> 904     survival       National       85                  0 0.81337307
#> 905     survival       National       96                  0 0.81402854
#> 906     survival       National       63                  0 0.83393234
#> 907     survival       National       30                  0 0.85702460
#> 908     survival       National       89                  0 0.81461149
#> 909     survival       National      100                  0 0.80936885
#>            lower     upper probViable                 Metric
#> 102 0.0531343927 0.4042838      0.000   Expected recruitment
#> 103 0.1116219558 0.5109438      0.000   Expected recruitment
#> 104 0.0956029348 0.4817734      0.000   Expected recruitment
#> 105 0.1081604077 0.5031764      0.000   Expected recruitment
#> 106 0.1271782628 0.4972864      0.000   Expected recruitment
#> 107 0.0592024556 0.3963140      0.000   Expected recruitment
#> 108 0.0765572931 0.4376818      0.000   Expected recruitment
#> 109 0.1016268727 0.4708860      0.000   Expected recruitment
#> 110 0.0814185492 0.4605223      0.000   Expected recruitment
#> 111 0.1003288456 0.5108371      0.000   Expected recruitment
#> 112 0.1331902473 0.5259641      0.000   Expected recruitment
#> 113 0.0863207477 0.4610402      0.000   Expected recruitment
#> 114 0.1183461491 0.5248124      0.000   Expected recruitment
#> 115 0.0938703733 0.4697089      0.000   Expected recruitment
#> 116 0.1085380084 0.5025087      0.000   Expected recruitment
#> 117 0.0805400479 0.4475746      0.000   Expected recruitment
#> 118 0.0512357720 0.3888855      0.000   Expected recruitment
#> 119 0.0900330776 0.4448334      0.000   Expected recruitment
#> 120 0.1022899541 0.4636047      0.000   Expected recruitment
#> 121 0.1399293799 0.5774776      0.000   Expected recruitment
#> 122 0.1578680504 0.5900724      0.000   Expected recruitment
#> 123 0.0325875682 0.3446916      0.000   Expected recruitment
#> 124 0.0236271748 0.3088848      0.000   Expected recruitment
#> 125 0.0563439709 0.4069945      0.000   Expected recruitment
#> 126 0.1600303963 0.5641602      0.000   Expected recruitment
#> 127 0.1079292021 0.4987087      0.000   Expected recruitment
#> 128 0.0767278918 0.4256221      0.000   Expected recruitment
#> 129 0.0542393522 0.3938058      0.000   Expected recruitment
#> 130 0.1304944970 0.5247123      0.000   Expected recruitment
#> 131 0.0223546342 0.3056294      0.000   Expected recruitment
#> 132 0.1419461425 0.5395569      0.000   Expected recruitment
#> 133 0.1564095506 0.5733543      0.000   Expected recruitment
#> 134 0.0298888840 0.3461829      0.000   Expected recruitment
#> 135 0.0193253734 0.3088619      0.000   Expected recruitment
#> 136 0.0132799550 0.2663303      0.000   Expected recruitment
#> 137 0.1467354335 0.5501956      0.000   Expected recruitment
#> 138 0.0243928145 0.3134213      0.000   Expected recruitment
#> 139 0.0728635804 0.4311046      0.000   Expected recruitment
#> 140 0.0523870025 0.3897727      0.000   Expected recruitment
#> 141 0.0361687082 0.3621542      0.000   Expected recruitment
#> 142 0.0209437282 0.3139605      0.000   Expected recruitment
#> 143 0.0610288304 0.4145533      0.000   Expected recruitment
#> 144 0.1514599743 0.5640116      0.000   Expected recruitment
#> 145 0.0301923820 0.3264388      0.000   Expected recruitment
#> 146 0.0210061530 0.3015278      0.000   Expected recruitment
#> 147 0.0127442972 0.2628834      0.000   Expected recruitment
#> 148 0.0368188550 0.3525461      0.000   Expected recruitment
#> 149 0.0289282168 0.3289770      0.000   Expected recruitment
#> 150 0.0688254102 0.4312571      0.000   Expected recruitment
#> 151 0.0539670048 0.4049691      0.000   Expected recruitment
#> 152 0.0310123370 0.3456386      0.000   Expected recruitment
#> 153 0.0230813449 0.3030542      0.000   Expected recruitment
#> 154 0.0601479407 0.4121221      0.000   Expected recruitment
#> 155 0.1543160280 0.5683337      0.000   Expected recruitment
#> 156 0.0292487870 0.3175212      0.000   Expected recruitment
#> 157 0.0178918348 0.2872474      0.000   Expected recruitment
#> 158 0.0136847056 0.2582561      0.000   Expected recruitment
#> 159 0.0340193033 0.3474552      0.000   Expected recruitment
#> 160 0.0269273028 0.3153632      0.000   Expected recruitment
#> 161 0.0765932575 0.4221466      0.000   Expected recruitment
#> 162 0.0505599980 0.3845220      0.000   Expected recruitment
#> 163 0.0056755127 0.2230533      0.000   Expected recruitment
#> 164 0.0399771839 0.3470584      0.000   Expected recruitment
#> 165 0.0129076344 0.2575694      0.000   Expected recruitment
#> 166 0.0414574260 0.3582559      0.000   Expected recruitment
#> 167 0.0303942198 0.3294860      0.000   Expected recruitment
#> 168 0.0194086473 0.2785746      0.000   Expected recruitment
#> 169 0.0114636040 0.2678823      0.000   Expected recruitment
#> 170 0.0368174515 0.3464185      0.000   Expected recruitment
#> 171 0.0028772403 0.2075875      0.000   Expected recruitment
#> 172 0.0613915790 0.3990953      0.000   Expected recruitment
#> 173 0.0527002464 0.3736949      0.000   Expected recruitment
#> 174 0.0051551452 0.2316132      0.000   Expected recruitment
#> 175 0.0081414265 0.2293643      0.000   Expected recruitment
#> 176 0.0022527271 0.2108627      0.000   Expected recruitment
#> 177 0.0153786286 0.2836552      0.000   Expected recruitment
#> 178 0.0034798559 0.2279330      0.000   Expected recruitment
#> 179 0.0174998170 0.2822789      0.000   Expected recruitment
#> 180 0.0117606829 0.2614358      0.000   Expected recruitment
#> 181 0.0069600853 0.2264994      0.000   Expected recruitment
#> 182 0.0042759867 0.2108679      0.000   Expected recruitment
#> 183 0.0137429611 0.2761293      0.000   Expected recruitment
#> 184 0.0426373807 0.3695786      0.000   Expected recruitment
#> 185 0.0049317310 0.2344996      0.000   Expected recruitment
#> 186 0.0061064248 0.2401856      0.000   Expected recruitment
#> 187 0.0026393946 0.2165766      0.000   Expected recruitment
#> 188 0.0158285123 0.2740653      0.000   Expected recruitment
#> 189 0.0091135919 0.2466575      0.000   Expected recruitment
#> 190 0.0109703234 0.2709054      0.000   Expected recruitment
#> 191 0.0103094348 0.2476015      0.000   Expected recruitment
#> 192 0.0055304317 0.2510843      0.000   Expected recruitment
#> 193 0.0065244567 0.2422595      0.000   Expected recruitment
#> 194 0.0095295129 0.2459630      0.000   Expected recruitment
#> 195 0.0478921067 0.3677934      0.000   Expected recruitment
#> 196 0.0046754453 0.2132808      0.000   Expected recruitment
#> 197 0.0048703719 0.2335431      0.000   Expected recruitment
#> 198 0.0133712881 0.2808699      0.000   Expected recruitment
#> 199 0.0090628558 0.2553122      0.000   Expected recruitment
#> 200 0.0087967622 0.2318047      0.000   Expected recruitment
#> 201 0.0043215496 0.2156538      0.000   Expected recruitment
#> 202 0.0065060171 0.2179664      0.000   Expected recruitment
#> 203 0.7732218484 0.9352439      0.000      Expected survival
#> 204 0.7522205384 0.9299177      0.000      Expected survival
#> 205 0.7717329104 0.9379788      0.000      Expected survival
#> 206 0.7708581820 0.9406360      0.000      Expected survival
#> 207 0.7949553200 0.9446645      0.000      Expected survival
#> 208 0.7727848622 0.9377910      0.000      Expected survival
#> 209 0.7650282426 0.9334503      0.000      Expected survival
#> 210 0.7742883430 0.9362277      0.000      Expected survival
#> 211 0.7657004419 0.9314726      0.000      Expected survival
#> 212 0.7683397164 0.9404786      0.000      Expected survival
#> 213 0.7656589759 0.9377041      0.000      Expected survival
#> 214 0.7867238922 0.9465088      0.000      Expected survival
#> 215 0.7382423515 0.9227797      0.000      Expected survival
#> 216 0.7293363759 0.9151857      0.000      Expected survival
#> 217 0.7500177669 0.9233867      0.000      Expected survival
#> 218 0.7800659319 0.9476354      0.000      Expected survival
#> 219 0.7792143313 0.9381489      0.000      Expected survival
#> 220 0.7585005610 0.9331984      0.000      Expected survival
#> 221 0.7527059885 0.9281738      0.000      Expected survival
#> 222 0.7573485108 0.9285205      0.000      Expected survival
#> 223 0.7656575896 0.9347472      0.000      Expected survival
#> 224 0.7779537017 0.9410517      0.000      Expected survival
#> 225 0.7795816892 0.9455721      0.000      Expected survival
#> 226 0.7350414521 0.9183723      0.000      Expected survival
#> 227 0.7264386399 0.9194991      0.000      Expected survival
#> 228 0.7577082370 0.9250550      0.000      Expected survival
#> 229 0.7829012445 0.9394527      0.000      Expected survival
#> 230 0.7722109612 0.9379229      0.000      Expected survival
#> 231 0.7644154954 0.9327486      0.000      Expected survival
#> 232 0.7447700400 0.9278095      0.000      Expected survival
#> 233 0.7327369928 0.9179318      0.000      Expected survival
#> 234 0.7643820493 0.9353862      0.000      Expected survival
#> 235 0.7705001658 0.9402270      0.000      Expected survival
#> 236 0.7818248825 0.9462017      0.000      Expected survival
#> 237 0.7388344394 0.9138155      0.000      Expected survival
#> 238 0.7218619304 0.9140383      0.000      Expected survival
#> 239 0.7474985959 0.9238597      0.000      Expected survival
#> 240 0.7853583000 0.9396912      0.000      Expected survival
#> 241 0.7633694020 0.9376032      0.000      Expected survival
#> 242 0.7657436800 0.9303765      0.000      Expected survival
#> 243 0.7502211409 0.9297629      0.000      Expected survival
#> 244 0.7137651321 0.9037580      0.000      Expected survival
#> 245 0.7321369112 0.9193236      0.000      Expected survival
#> 246 0.7452274759 0.9324492      0.000      Expected survival
#> 247 0.7784969233 0.9440757      0.000      Expected survival
#> 248 0.7385312256 0.9179329      0.000      Expected survival
#> 249 0.7299338284 0.9086562      0.000      Expected survival
#> 250 0.7628559269 0.9314792      0.000      Expected survival
#> 251 0.7390156135 0.9190277      0.000      Expected survival
#> 252 0.7801254309 0.9335553      0.000      Expected survival
#> 253 0.7574685872 0.9306024      0.000      Expected survival
#> 254 0.7465983630 0.9255395      0.000      Expected survival
#> 255 0.7060233608 0.8999408      0.000      Expected survival
#> 256 0.7475674168 0.9205015      0.000      Expected survival
#> 257 0.7434912911 0.9321024      0.000      Expected survival
#> 258 0.7435549381 0.9229355      0.000      Expected survival
#> 259 0.7304098067 0.9197375      0.000      Expected survival
#> 260 0.7204794897 0.9164797      0.000      Expected survival
#> 261 0.7199427820 0.9077695      0.000      Expected survival
#> 262 0.7449476501 0.9191905      0.000      Expected survival
#> 263 0.7308424930 0.9156635      0.000      Expected survival
#> 264 0.7549262091 0.9346699      0.000      Expected survival
#> 265 0.7431936049 0.9178390      0.000      Expected survival
#> 266 0.6973424673 0.9029509      0.000      Expected survival
#> 267 0.7418481351 0.9236249      0.000      Expected survival
#> 268 0.7252180949 0.9104136      0.000      Expected survival
#> 269 0.7176334297 0.9114073      0.000      Expected survival
#> 270 0.7256712231 0.9156373      0.000      Expected survival
#> 271 0.7213074460 0.9136882      0.000      Expected survival
#> 272 0.7154038364 0.9093259      0.000      Expected survival
#> 273 0.7048231963 0.9041788      0.000      Expected survival
#> 274 0.7279168746 0.9160079      0.000      Expected survival
#> 275 0.7495582921 0.9312466      0.000      Expected survival
#> 276 0.7538289424 0.9237826      0.000      Expected survival
#> 277 0.7076567818 0.9060399      0.000      Expected survival
#> 278 0.7006645418 0.9030260      0.000      Expected survival
#> 279 0.7189271544 0.9085279      0.000      Expected survival
#> 280 0.7211165189 0.9112812      0.000      Expected survival
#> 281 0.7338917813 0.9204431      0.000      Expected survival
#> 282 0.7028326843 0.8958275      0.000      Expected survival
#> 283 0.7143773602 0.9071711      0.000      Expected survival
#> 284 0.7172696298 0.9077104      0.000      Expected survival
#> 285 0.7115698800 0.9055099      0.000      Expected survival
#> 286 0.7223870890 0.9092640      0.000      Expected survival
#> 287 0.7506501529 0.9255955      0.000      Expected survival
#> 288 0.7017765536 0.9033797      0.000      Expected survival
#> 289 0.7062084967 0.9009861      0.000      Expected survival
#> 290 0.7202731460 0.9103699      0.000      Expected survival
#> 291 0.7336683717 0.9139495      0.000      Expected survival
#> 292 0.7328932262 0.9146341      0.000      Expected survival
#> 293 0.7040278507 0.9009843      0.000      Expected survival
#> 294 0.7086484662 0.9077684      0.000      Expected survival
#> 295 0.7063618919 0.8980982      0.000      Expected survival
#> 296 0.7077519543 0.8980725      0.000      Expected survival
#> 297 0.7091187931 0.9050399      0.000      Expected survival
#> 298 0.7252195313 0.9129598      0.000      Expected survival
#> 299 0.7035243269 0.9001785      0.000      Expected survival
#> 300 0.7158121801 0.9103363      0.000      Expected survival
#> 301 0.7132413408 0.9049926      0.000      Expected survival
#> 302 0.7065660715 0.9035824      0.000      Expected survival
#> 303 0.7175110444 0.9053626      0.000      Expected survival
#> 304 0.0321735389 0.3366164      0.000   Adjusted recruitment
#> 305 0.0292063128 0.3242002      0.000   Adjusted recruitment
#> 306 0.0278694488 0.3281632      0.000   Adjusted recruitment
#> 307 0.0324138938 0.3459531      0.000   Adjusted recruitment
#> 308 0.0285822885 0.3402281      0.000   Adjusted recruitment
#> 309 0.0101774252 0.2225546      0.000   Adjusted recruitment
#> 310 0.0374196456 0.3729945      0.000   Adjusted recruitment
#> 311 0.0269786767 0.3211466      0.000   Adjusted recruitment
#> 312 0.0389145133 0.3675657      0.000   Adjusted recruitment
#> 313 0.0171568098 0.2621045      0.000   Adjusted recruitment
#> 314 0.0201226674 0.2910190      0.000   Adjusted recruitment
#> 315 0.0306010378 0.3181265      0.000   Adjusted recruitment
#> 316 0.0209677599 0.2957152      0.000   Adjusted recruitment
#> 317 0.0340427777 0.3569081      0.000   Adjusted recruitment
#> 318 0.0094401676 0.2027218      0.000   Adjusted recruitment
#> 319 0.0447392711 0.3646461      0.000   Adjusted recruitment
#> 320 0.0108357058 0.2318503      0.000   Adjusted recruitment
#> 321 0.0241973213 0.2885077      0.000   Adjusted recruitment
#> 322 0.0303166480 0.3026923      0.000   Adjusted recruitment
#> 323 0.0206414564 0.2880254      0.000   Adjusted recruitment
#> 324 0.0152775800 0.2408822      0.000   Adjusted recruitment
#> 325 0.0107131096 0.2070858      0.000   Adjusted recruitment
#> 326 0.0213980373 0.3131796      0.000   Adjusted recruitment
#> 327 0.0339597556 0.3520612      0.000   Adjusted recruitment
#> 328 0.0399290797 0.3573161      0.000   Adjusted recruitment
#> 329 0.0092915036 0.2166781      0.000   Adjusted recruitment
#> 330 0.0205034937 0.2937772      0.000   Adjusted recruitment
#> 331 0.0153972586 0.2619478      0.000   Adjusted recruitment
#> 332 0.0230162860 0.3053503      0.000   Adjusted recruitment
#> 333 0.0271843597 0.3316981      0.000   Adjusted recruitment
#> 334 0.0206096989 0.2847776      0.000   Adjusted recruitment
#> 335 0.0166820073 0.2217627      0.000   Adjusted recruitment
#> 336 0.0018167402 0.1246639      0.000   Adjusted recruitment
#> 337 0.0235305142 0.3010318      0.000   Adjusted recruitment
#> 338 0.0166093987 0.2595889      0.000   Adjusted recruitment
#> 339 0.0375151863 0.3503686      0.000   Adjusted recruitment
#> 340 0.0078922531 0.2090131      0.000   Adjusted recruitment
#> 341 0.0061318933 0.1772524      0.000   Adjusted recruitment
#> 342 0.0166064999 0.2594740      0.000   Adjusted recruitment
#> 343 0.0320591361 0.3450418      0.000   Adjusted recruitment
#> 344 0.0248679537 0.3114203      0.000   Adjusted recruitment
#> 345 0.0165571423 0.2776860      0.000   Adjusted recruitment
#> 346 0.0119496560 0.2382067      0.000   Adjusted recruitment
#> 347 0.0021374052 0.1450186      0.000   Adjusted recruitment
#> 348 0.0062549697 0.1901536      0.000   Adjusted recruitment
#> 349 0.0181349734 0.2607805      0.000   Adjusted recruitment
#> 350 0.0396235641 0.3568864      0.000   Adjusted recruitment
#> 351 0.0106669947 0.2125299      0.000   Adjusted recruitment
#> 352 0.0060284716 0.1780960      0.000   Adjusted recruitment
#> 353 0.0037650779 0.1645509      0.000   Adjusted recruitment
#> 354 0.0105576891 0.2185518      0.000   Adjusted recruitment
#> 355 0.0240644860 0.3100132      0.000   Adjusted recruitment
#> 356 0.0187172228 0.2710773      0.000   Adjusted recruitment
#> 357 0.0140423609 0.2397324      0.000   Adjusted recruitment
#> 358 0.0016708416 0.1444160      0.000   Adjusted recruitment
#> 359 0.0060607120 0.1751310      0.000   Adjusted recruitment
#> 360 0.0147156272 0.2645807      0.000   Adjusted recruitment
#> 361 0.0362629988 0.3554705      0.000   Adjusted recruitment
#> 362 0.0104791014 0.2118674      0.000   Adjusted recruitment
#> 363 0.0056395710 0.1702886      0.000   Adjusted recruitment
#> 364 0.0032264718 0.1695680      0.000   Adjusted recruitment
#> 365 0.0019915940 0.1329715      0.000   Adjusted recruitment
#> 366 0.0070415103 0.2137062      0.000   Adjusted recruitment
#> 367 0.0183427246 0.2729477      0.000   Adjusted recruitment
#> 368 0.0138232296 0.2334781      0.000   Adjusted recruitment
#> 369 0.0015052640 0.1342821      0.000   Adjusted recruitment
#> 370 0.0053167399 0.1817718      0.000   Adjusted recruitment
#> 371 0.0043375141 0.1744025      0.000   Adjusted recruitment
#> 372 0.0306114972 0.3424903      0.000   Adjusted recruitment
#> 373 0.0079941750 0.2086809      0.000   Adjusted recruitment
#> 374 0.0052615471 0.1681402      0.000   Adjusted recruitment
#> 375 0.0032333774 0.1637945      0.000   Adjusted recruitment
#> 376 0.0050861791 0.1699462      0.000   Adjusted recruitment
#> 377 0.0063669096 0.2063902      0.000   Adjusted recruitment
#> 378 0.0048458867 0.1691529      0.000   Adjusted recruitment
#> 379 0.0120732185 0.2310371      0.000   Adjusted recruitment
#> 380 0.0014383719 0.1232021      0.000   Adjusted recruitment
#> 381 0.0021224692 0.1604594      0.000   Adjusted recruitment
#> 382 0.0045115091 0.1702972      0.000   Adjusted recruitment
#> 383 0.0104680885 0.2440307      0.000   Adjusted recruitment
#> 384 0.0084188592 0.1956433      0.000   Adjusted recruitment
#> 385 0.0044002119 0.1681555      0.000   Adjusted recruitment
#> 386 0.0037862743 0.1432308      0.000   Adjusted recruitment
#> 387 0.0051061912 0.1508803      0.000   Adjusted recruitment
#> 388 0.0117475402 0.2166163      0.000   Adjusted recruitment
#> 389 0.0010940486 0.1226220      0.000   Adjusted recruitment
#> 390 0.0123809910 0.2299037      0.000   Adjusted recruitment
#> 391 0.0009417492 0.1290865      0.000   Adjusted recruitment
#> 392 0.0043352878 0.1693930      0.000   Adjusted recruitment
#> 393 0.0023800118 0.1498169      0.000   Adjusted recruitment
#> 394 0.0050552898 0.1730505      0.000   Adjusted recruitment
#> 395 0.0083537309 0.2123416      0.000   Adjusted recruitment
#> 396 0.0009213603 0.1185853      0.000   Adjusted recruitment
#> 397 0.0029695203 0.1619735      0.000   Adjusted recruitment
#> 398 0.0017459689 0.1242309      0.000   Adjusted recruitment
#> 399 0.0107114435 0.2292046      0.000   Adjusted recruitment
#> 400 0.0029608209 0.1459668      0.000   Adjusted recruitment
#> 401 0.0021132505 0.1444254      0.000   Adjusted recruitment
#> 402 0.0022617569 0.1514023      0.000   Adjusted recruitment
#> 403 0.0025271628 0.1601807      0.000   Adjusted recruitment
#> 404 0.0023887014 0.1380607      0.000   Adjusted recruitment
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
#> 506 0.7499750000 1.2130500      0.488 Population growth rate
#> 507 0.7520000000 1.2590750      0.570 Population growth rate
#> 508 0.7549750000 1.2130750      0.478 Population growth rate
#> 509 0.7400000000 1.2072000      0.466 Population growth rate
#> 510 0.7360000000 1.1770000      0.340 Population growth rate
#> 511 0.7629500000 1.2560500      0.505 Population growth rate
#> 512 0.7749750000 1.2650000      0.590 Population growth rate
#> 513 0.7049750000 1.1010000      0.175 Population growth rate
#> 514 0.7460000000 1.2090250      0.443 Population growth rate
#> 515 0.7599250000 1.2140500      0.467 Population growth rate
#> 516 0.7599750000 1.2600500      0.555 Population growth rate
#> 517 0.7400000000 1.2121250      0.430 Population growth rate
#> 518 0.7380000000 1.1820250      0.372 Population growth rate
#> 519 0.7420000000 1.2070250      0.444 Population growth rate
#> 520 0.7389750000 1.1600250      0.360 Population growth rate
#> 521 0.7929750000 1.2750000      0.606 Population growth rate
#> 522 0.7569750000 1.2330500      0.512 Population growth rate
#> 523 0.7799750000 1.2651000      0.565 Population growth rate
#> 524 0.7089000000 1.0920500      0.177 Population growth rate
#> 525 0.7379250000 1.1840000      0.385 Population growth rate
#> 526 0.7989750000 1.2731250      0.623 Population growth rate
#> 527 0.7619250000 1.2320250      0.513 Population growth rate
#> 528 0.7327250000 1.1990250      0.426 Population growth rate
#> 529 0.7200000000 1.1620000      0.339 Population growth rate
#> 530 0.7899500000 1.2740750      0.612 Population growth rate
#> 531 0.7509750000 1.1500500      0.350 Population growth rate
#> 532 0.7789000000 1.2653250      0.599 Population growth rate
#> 533 0.7620000000 1.2330250      0.518 Population growth rate
#> 534 0.7679750000 1.2430500      0.593 Population growth rate
#> 535 0.6959750000 1.0681000      0.168 Population growth rate
#> 536 0.6828750000 1.0500750      0.107 Population growth rate
#> 537 0.7919250000 1.2780500      0.626 Population growth rate
#> 538 0.7089500000 1.1081250      0.204 Population growth rate
#> 539 0.7549250000 1.1950500      0.417 Population growth rate
#> 540 0.7329750000 1.1510750      0.307 Population growth rate
#> 541 0.7069750000 1.1180250      0.267 Population growth rate
#> 542 0.7010000000 1.0870250      0.177 Population growth rate
#> 543 0.7430000000 1.1660000      0.383 Population growth rate
#> 544 0.6959750000 1.0860750      0.174 Population growth rate
#> 545 0.7517500000 1.2390000      0.541 Population growth rate
#> 546 0.6979750000 1.0760250      0.152 Population growth rate
#> 547 0.6939500000 1.0520250      0.114 Population growth rate
#> 548 0.7219750000 1.1250000      0.259 Population growth rate
#> 549 0.7049750000 1.1081000      0.206 Population growth rate
#> 550 0.7419750000 1.1880000      0.412 Population growth rate
#> 551 0.7219750000 1.1600750      0.290 Population growth rate
#> 552 0.7170000000 1.1100500      0.243 Population growth rate
#> 553 0.6639500000 1.0020750      0.040 Population growth rate
#> 554 0.6700000000 1.0610750      0.116 Population growth rate
#> 555 0.7199500000 1.1280250      0.299 Population growth rate
#> 556 0.7099750000 1.1040000      0.212 Population growth rate
#> 557 0.7039250000 1.0720250      0.153 Population growth rate
#> 558 0.6750000000 1.0390500      0.118 Population growth rate
#> 559 0.7169750000 1.1380500      0.274 Population growth rate
#> 560 0.6859500000 1.0990500      0.176 Population growth rate
#> 561 0.7430000000 1.1910000      0.427 Population growth rate
#> 562 0.7099750000 1.1581000      0.323 Population growth rate
#> 563 0.6970000000 1.1090250      0.252 Population growth rate
#> 564 0.6599500000 1.0020000      0.039 Population growth rate
#> 565 0.6919000000 1.0480750      0.105 Population growth rate
#> 566 0.7049750000 1.1360500      0.284 Population growth rate
#> 567 0.7099500000 1.1080250      0.214 Population growth rate
#> 568 0.6959500000 1.0790250      0.160 Population growth rate
#> 569 0.6849750000 1.0490250      0.096 Population growth rate
#> 570 0.7089500000 1.1180250      0.247 Population growth rate
#> 571 0.6638500000 1.0290500      0.064 Population growth rate
#> 572 0.7349500000 1.1761250      0.379 Population growth rate
#> 573 0.7309500000 1.1350250      0.307 Population growth rate
#> 574 0.7139750000 1.1180250      0.233 Population growth rate
#> 575 0.6609750000 1.0110000      0.046 Population growth rate
#> 576 0.6558750000 1.0230250      0.057 Population growth rate
#> 577 0.7140000000 1.1260750      0.264 Population growth rate
#> 578 0.6689250000 1.0130000      0.048 Population growth rate
#> 579 0.6808750000 1.0850750      0.144 Population growth rate
#> 580 0.6880000000 1.0460000      0.099 Population growth rate
#> 581 0.6580000000 1.0250750      0.062 Population growth rate
#> 582 0.6730000000 1.0110000      0.044 Population growth rate
#> 583 0.6909250000 1.0610750      0.121 Population growth rate
#> 584 0.7099500000 1.1380000      0.316 Population growth rate
#> 585 0.6890000000 1.1100000      0.225 Population growth rate
#> 586 0.6679750000 1.0010500      0.038 Population growth rate
#> 587 0.6589250000 1.0270000      0.056 Population growth rate
#> 588 0.6779750000 1.0260000      0.070 Population growth rate
#> 589 0.6730000000 1.0411000      0.087 Population growth rate
#> 590 0.6809750000 1.0741750      0.137 Population growth rate
#> 591 0.6849250000 1.0500250      0.096 Population growth rate
#> 592 0.6789750000 1.0180250      0.057 Population growth rate
#> 593 0.6729250000 1.0370750      0.081 Population growth rate
#> 594 0.6810000000 1.0170250      0.054 Population growth rate
#> 595 0.6769500000 1.0450250      0.101 Population growth rate
#> 596 0.6820000000 1.0170250      0.060 Population growth rate
#> 597 0.6730000000 1.0330000      0.082 Population growth rate
#> 598 0.6900000000 1.0620500      0.100 Population growth rate
#> 599 0.6769500000 1.0380250      0.078 Population growth rate
#> 600 0.6790000000 1.0380500      0.093 Population growth rate
#> 601 0.6950000000 1.0710000      0.118 Population growth rate
#> 602 0.6739250000 1.0150750      0.056 Population growth rate
#> 603 0.6769500000 1.0330500      0.068 Population growth rate
#> 604 0.6709500000 1.0310250      0.089 Population growth rate
#> 605 0.6719500000 1.0160000      0.057 Population growth rate
#> 606 0.6909750000 1.0660500      0.138 Population growth rate
#> 607 0.8519262982 1.1040938      0.405   Expected growth rate
#> 608 0.8490449680 1.1113474      0.462   Expected growth rate
#> 609 0.8554360625 1.1085225      0.438   Expected growth rate
#> 610 0.8840980415 1.1509566      0.631   Expected growth rate
#> 611 0.8456105553 1.0972394      0.366   Expected growth rate
#> 612 0.8646727244 1.1323476      0.545   Expected growth rate
#> 613 0.8388654875 1.0964243      0.354   Expected growth rate
#> 614 0.8703630690 1.1330570      0.576   Expected growth rate
#> 615 0.8772628378 1.1382397      0.562   Expected growth rate
#> 616 0.7721351149 1.0071148      0.045   Expected growth rate
#> 617 0.8466933442 1.0961397      0.382   Expected growth rate
#> 618 0.8927327316 1.1711593      0.729   Expected growth rate
#> 619 0.8583555760 1.1212037      0.529   Expected growth rate
#> 620 0.8331902984 1.0843186      0.322   Expected growth rate
#> 621 0.8858681807 1.1389842      0.624   Expected growth rate
#> 622 0.8885621513 1.1471546      0.662   Expected growth rate
#> 623 0.8619320471 1.1261643      0.493   Expected growth rate
#> 624 0.8913867016 1.1592897      0.714   Expected growth rate
#> 625 0.8101747534 1.0556850      0.205   Expected growth rate
#> 626 0.8682164170 1.1311137      0.564   Expected growth rate
#> 627 0.7752603835 1.0094486      0.043   Expected growth rate
#> 628 0.8175860876 1.0692730      0.230   Expected growth rate
#> 629 0.8934332218 1.1588102      0.725   Expected growth rate
#> 630 0.8553521803 1.1253093      0.486   Expected growth rate
#> 631 0.8405715587 1.0884500      0.351   Expected growth rate
#> 632 0.8196813335 1.0535101      0.171   Expected growth rate
#> 633 0.7937091779 1.0240379      0.102   Expected growth rate
#> 634 0.7695181747 1.0183293      0.069   Expected growth rate
#> 635 0.8414349440 1.0785352      0.276   Expected growth rate
#> 636 0.8005232342 1.0553627      0.160   Expected growth rate
#> 637 0.7832283577 1.0246652      0.071   Expected growth rate
#> 638 0.7670770243 1.0012637      0.040   Expected growth rate
#> 639 0.8176355049 1.0588823      0.222   Expected growth rate
#> 640 0.8960049370 1.1645151      0.684   Expected growth rate
#> 641 0.8586762867 1.1160153      0.505   Expected growth rate
#> 642 0.8362541701 1.0839019      0.338   Expected growth rate
#> 643 0.8105264713 1.0539014      0.173   Expected growth rate
#> 644 0.7871712729 1.0395553      0.112   Expected growth rate
#> 645 0.7781933217 1.0071903      0.050   Expected growth rate
#> 646 0.8310176415 1.0769562      0.267   Expected growth rate
#> 647 0.8004584863 1.0514990      0.155   Expected growth rate
#> 648 0.7918785054 1.0158414      0.071   Expected growth rate
#> 649 0.7723528630 0.9991521      0.037   Expected growth rate
#> 650 0.8153017833 1.0616033      0.199   Expected growth rate
#> 651 0.9010820690 1.1677218      0.741   Expected growth rate
#> 652 0.8539245391 1.1236282      0.481   Expected growth rate
#> 653 0.8320953250 1.0797383      0.285   Expected growth rate
#> 654 0.8198904421 1.0605672      0.176   Expected growth rate
#> 655 0.7815855881 1.0277390      0.088   Expected growth rate
#> 656 0.7190526188 0.9432136      0.003   Expected growth rate
#> 657 0.8148221363 1.0742262      0.227   Expected growth rate
#> 658 0.7920385291 1.0416467      0.132   Expected growth rate
#> 659 0.7774348424 1.0206409      0.071   Expected growth rate
#> 660 0.7621555981 0.9951630      0.030   Expected growth rate
#> 661 0.7504306114 0.9824657      0.014   Expected growth rate
#> 662 0.7919031631 1.0259442      0.097   Expected growth rate
#> 663 0.7798342769 1.0155970      0.061   Expected growth rate
#> 664 0.8254325607 1.0749981      0.288   Expected growth rate
#> 665 0.8067377211 1.0604118      0.153   Expected growth rate
#> 666 0.7867758100 1.0414909      0.096   Expected growth rate
#> 667 0.7217005057 0.9410478      0.004   Expected growth rate
#> 668 0.7640584342 0.9944242      0.030   Expected growth rate
#> 669 0.7983357261 1.0382117      0.108   Expected growth rate
#> 670 0.7804199174 1.0158412      0.073   Expected growth rate
#> 671 0.7641895052 0.9941983      0.032   Expected growth rate
#> 672 0.7535579975 0.9790636      0.016   Expected growth rate
#> 673 0.7275739213 0.9512094      0.004   Expected growth rate
#> 674 0.7285208220 0.9576338      0.005   Expected growth rate
#> 675 0.7529287712 0.9839417      0.023   Expected growth rate
#> 676 0.7366449377 0.9720666      0.012   Expected growth rate
#> 677 0.7281653298 0.9543661      0.008   Expected growth rate
#> 678 0.8057129329 1.0401763      0.115   Expected growth rate
#> 679 0.7511392785 0.9847515      0.020   Expected growth rate
#> 680 0.8074217348 1.0456170      0.120   Expected growth rate
#> 681 0.7792016709 1.0061239      0.051   Expected growth rate
#> 682 0.7570221765 0.9958527      0.029   Expected growth rate
#> 683 0.7456087178 0.9719802      0.011   Expected growth rate
#> 684 0.7330210699 0.9595825      0.011   Expected growth rate
#> 685 0.7201323644 0.9483986      0.002   Expected growth rate
#> 686 0.7465735750 0.9839515      0.022   Expected growth rate
#> 687 0.7413965974 0.9723522      0.010   Expected growth rate
#> 688 0.7259455633 0.9539558      0.005   Expected growth rate
#> 689 0.7364656121 0.9558698      0.003   Expected growth rate
#> 690 0.7539351998 0.9721576      0.009   Expected growth rate
#> 691 0.7565201202 0.9908746      0.026   Expected growth rate
#> 692 0.7739681657 1.0148965      0.060   Expected growth rate
#> 693 0.7685415488 0.9871079      0.022   Expected growth rate
#> 694 0.7434961469 0.9713776      0.011   Expected growth rate
#> 695 0.7301114011 0.9544649      0.010   Expected growth rate
#> 696 0.7439409259 0.9735959      0.013   Expected growth rate
#> 697 0.7535238971 0.9858628      0.022   Expected growth rate
#> 698 0.7294781589 0.9664497      0.007   Expected growth rate
#> 699 0.7315312147 0.9467174      0.003   Expected growth rate
#> 700 0.7273891434 0.9513384      0.004   Expected growth rate
#> 701 0.7480686017 0.9834813      0.020   Expected growth rate
#> 702 0.7352199451 0.9530471      0.007   Expected growth rate
#> 703 0.7282163004 0.9504685      0.005   Expected growth rate
#> 704 0.7151784199 0.9579284      0.006   Expected growth rate
#> 705 0.7451122405 0.9663586      0.006   Expected growth rate
#> 706 0.7376227976 0.9733396      0.010   Expected growth rate
#> 707 0.7250072370 0.9445479      0.008   Expected growth rate
#> 708 0.0612229944 0.6849807      0.000            Recruitment
#> 709 0.0543687194 0.6633961      0.000            Recruitment
#> 710 0.0792471282 0.7137728      0.000            Recruitment
#> 711 0.0798581595 0.7146323      0.000            Recruitment
#> 712 0.0497359075 0.6228405      0.000            Recruitment
#> 713 0.0362699468 0.5215611      0.000            Recruitment
#> 714 0.0725259976 0.7109410      0.000            Recruitment
#> 715 0.0641182723 0.6900835      0.000            Recruitment
#> 716 0.0483946426 0.5770155      0.000            Recruitment
#> 717 0.0343136196 0.5242091      0.000            Recruitment
#> 718 0.0894785423 0.7292923      0.000            Recruitment
#> 719 0.0557388976 0.6563265      0.000            Recruitment
#> 720 0.0750303725 0.7007372      0.000            Recruitment
#> 721 0.0584126255 0.6484005      0.000            Recruitment
#> 722 0.0679195112 0.7041224      0.000            Recruitment
#> 723 0.0481289719 0.6200265      0.000            Recruitment
#> 724 0.0294312545 0.5291613      0.000            Recruitment
#> 725 0.0680855553 0.7138161      0.000            Recruitment
#> 726 0.0571645770 0.6804562      0.000            Recruitment
#> 727 0.0460325720 0.6107006      0.000            Recruitment
#> 728 0.0305551600 0.4817645      0.000            Recruitment
#> 729 0.0214262193 0.4141717      0.000            Recruitment
#> 730 0.0606332960 0.6053846      0.000            Recruitment
#> 731 0.0211153781 0.4371036      0.000            Recruitment
#> 732 0.0539573534 0.6422932      0.000            Recruitment
#> 733 0.0643470777 0.6732327      0.000            Recruitment
#> 734 0.0427960747 0.6263592      0.000            Recruitment
#> 735 0.0307945171 0.5238955      0.000            Recruitment
#> 736 0.0214228870 0.4584092      0.000            Recruitment
#> 737 0.0612020755 0.6362529      0.000            Recruitment
#> 738 0.0410069874 0.5875544      0.000            Recruitment
#> 739 0.0333640147 0.4435254      0.000            Recruitment
#> 740 0.0203548503 0.4451092      0.000            Recruitment
#> 741 0.0127338192 0.4127804      0.000            Recruitment
#> 742 0.0332187975 0.5191779      0.000            Recruitment
#> 743 0.0778290266 0.7351315      0.000            Recruitment
#> 744 0.0648277877 0.6919062      0.000            Recruitment
#> 745 0.0470610285 0.6020637      0.000            Recruitment
#> 746 0.0332129998 0.5189480      0.000            Recruitment
#> 747 0.0209361769 0.4880614      0.000            Recruitment
#> 748 0.0168377184 0.3912865      0.000            Recruitment
#> 749 0.0402453348 0.5820380      0.000            Recruitment
#> 750 0.0238993121 0.4764134      0.000            Recruitment
#> 751 0.0216714116 0.4637007      0.000            Recruitment
#> 752 0.0125099394 0.3803072      0.000            Recruitment
#> 753 0.0331142845 0.5553721      0.000            Recruitment
#> 754 0.0247619820 0.4598073      0.000            Recruitment
#> 755 0.0213339895 0.4250598      0.000            Recruitment
#> 756 0.0120569432 0.3561919      0.000            Recruitment
#> 757 0.0075301558 0.3291018      0.000            Recruitment
#> 758 0.0042449384 0.3209188      0.000            Recruitment
#> 759 0.0167074617 0.4246833      0.000            Recruitment
#> 760 0.0419355198 0.5914304      0.000            Recruitment
#> 761 0.0280847219 0.4794649      0.000            Recruitment
#> 762 0.0188803353 0.4054436      0.000            Recruitment
#> 763 0.0121214239 0.3502621      0.000            Recruitment
#> 764 0.0374344455 0.5421546      0.000            Recruitment
#> 765 0.0234950805 0.4332325      0.000            Recruitment
#> 766 0.0209582027 0.4237347      0.000            Recruitment
#> 767 0.0112791421 0.3405771      0.000            Recruitment
#> 768 0.0064529436 0.3391360      0.000            Recruitment
#> 769 0.0039831880 0.2659430      0.000            Recruitment
#> 770 0.0140830205 0.4274125      0.000            Recruitment
#> 771 0.0412829128 0.5760508      0.000            Recruitment
#> 772 0.0276464592 0.4669562      0.000            Recruitment
#> 773 0.0185830072 0.4333561      0.000            Recruitment
#> 774 0.0106334797 0.3635436      0.000            Recruitment
#> 775 0.0366854491 0.5458954      0.000            Recruitment
#> 776 0.0748392911 0.7459891      0.000            Recruitment
#> 777 0.0159883501 0.4173619      0.000            Recruitment
#> 778 0.0105230942 0.3362804      0.000            Recruitment
#> 779 0.0064667548 0.3275890      0.000            Recruitment
#> 780 0.0036334804 0.2493279      0.000            Recruitment
#> 781 0.0059216418 0.2919336      0.000            Recruitment
#> 782 0.0412193978 0.5695551      0.000            Recruitment
#> 783 0.0241464370 0.4620743      0.000            Recruitment
#> 784 0.0157845063 0.4180261      0.000            Recruitment
#> 785 0.0122637866 0.3545048      0.000            Recruitment
#> 786 0.0090230183 0.3405943      0.000            Recruitment
#> 787 0.0047774028 0.2761214      0.000            Recruitment
#> 788 0.0034919378 0.2484617      0.000            Recruitment
#> 789 0.0088004238 0.3363110      0.000            Recruitment
#> 790 0.0075725485 0.2864617      0.000            Recruitment
#> 791 0.0042748104 0.2900372      0.000            Recruitment
#> 792 0.0045235138 0.3028046      0.000            Recruitment
#> 793 0.0101723583 0.3398924      0.000            Recruitment
#> 794 0.0101105795 0.3461009      0.000            Recruitment
#> 795 0.0018834983 0.2581731      0.000            Recruitment
#> 796 0.0096917735 0.3383058      0.000            Recruitment
#> 797 0.0047600237 0.2996338      0.000            Recruitment
#> 798 0.0033416833 0.2888319      0.000            Recruitment
#> 799 0.0050543256 0.3203614      0.000            Recruitment
#> 800 0.0086705756 0.3387860      0.000            Recruitment
#> 801 0.0059390405 0.3239470      0.000            Recruitment
#> 802 0.0021880972 0.2452439      0.000            Recruitment
#> 803 0.0042265010 0.2888508      0.000            Recruitment
#> 804 0.0102123824 0.3017606      0.000            Recruitment
#> 805 0.0030105279 0.2685641      0.000            Recruitment
#> 806 0.0028767437 0.2464041      0.000            Recruitment
#> 807 0.0018427205 0.2371706      0.000            Recruitment
#> 808 0.0086750281 0.3488049      0.000            Recruitment
#> 809 0.7004529166 0.9986809      0.072  Adult female survival
#> 810 0.6680744508 0.9718654      0.008  Adult female survival
#> 811 0.6562349431 0.9700473      0.003  Adult female survival
#> 812 0.7002128740 0.9959880      0.055  Adult female survival
#> 813 0.6744344931 0.9977777      0.057  Adult female survival
#> 814 0.6788871399 0.9958954      0.042  Adult female survival
#> 815 0.6730813524 0.9972438      0.051  Adult female survival
#> 816 0.6886202895 0.9954510      0.038  Adult female survival
#> 817 0.6820952000 0.9965031      0.051  Adult female survival
#> 818 0.6926783268 0.9958910      0.047  Adult female survival
#> 819 0.6576322446 0.9744411      0.006  Adult female survival
#> 820 0.6773970746 0.9891580      0.021  Adult female survival
#> 821 0.6706651030 0.9832498      0.009  Adult female survival
#> 822 0.6637748091 0.9845318      0.017  Adult female survival
#> 823 0.6909361777 0.9976784      0.062  Adult female survival
#> 824 0.6652324706 0.9864582      0.017  Adult female survival
#> 825 0.6956714195 0.9955864      0.048  Adult female survival
#> 826 0.6495665409 0.9735719      0.010  Adult female survival
#> 827 0.6821174169 0.9891626      0.024  Adult female survival
#> 828 0.6997459304 0.9982263      0.073  Adult female survival
#> 829 0.6839164542 0.9957217      0.046  Adult female survival
#> 830 0.6641704670 0.9776273      0.009  Adult female survival
#> 831 0.6659127721 0.9895443      0.024  Adult female survival
#> 832 0.6573957874 0.9864224      0.016  Adult female survival
#> 833 0.6873358439 0.9930816      0.034  Adult female survival
#> 834 0.6532950381 0.9790701      0.011  Adult female survival
#> 835 0.6776319098 0.9842778      0.016  Adult female survival
#> 836 0.6707435408 0.9849927      0.017  Adult female survival
#> 837 0.6857410671 0.9883011      0.024  Adult female survival
#> 838 0.6469211148 0.9743811      0.011  Adult female survival
#> 839 0.6747847887 0.9776885      0.008  Adult female survival
#> 840 0.6809910680 0.9967418      0.054  Adult female survival
#> 841 0.6457372664 0.9793204      0.008  Adult female survival
#> 842 0.6814957154 0.9894750      0.025  Adult female survival
#> 843 0.6692961655 0.9829290      0.012  Adult female survival
#> 844 0.6922440836 0.9942545      0.037  Adult female survival
#> 845 0.6538006288 0.9748738      0.010  Adult female survival
#> 846 0.6680851629 0.9849744      0.014  Adult female survival
#> 847 0.6564982582 0.9817332      0.013  Adult female survival
#> 848 0.6857461734 0.9900147      0.025  Adult female survival
#> 849 0.6531052884 0.9735214      0.007  Adult female survival
#> 850 0.6701443809 0.9855198      0.019  Adult female survival
#> 851 0.6973784511 0.9917656      0.030  Adult female survival
#> 852 0.6596845467 0.9692727      0.007  Adult female survival
#> 853 0.6626992276 0.9896217      0.024  Adult female survival
#> 854 0.6705724031 0.9828549      0.014  Adult female survival
#> 855 0.6836844804 0.9944839      0.030  Adult female survival
#> 856 0.6513704331 0.9750689      0.004  Adult female survival
#> 857 0.6747404541 0.9821991      0.015  Adult female survival
#> 858 0.6592036610 0.9780144      0.008  Adult female survival
#> 859 0.6690851526 0.9906431      0.026  Adult female survival
#> 860 0.6558603345 0.9703438      0.002  Adult female survival
#> 861 0.6529140218 0.9597749      0.004  Adult female survival
#> 862 0.6873552821 0.9973018      0.049  Adult female survival
#> 863 0.6695907253 0.9638342      0.004  Adult female survival
#> 864 0.6668075215 0.9874359      0.021  Adult female survival
#> 865 0.6632356350 0.9831685      0.011  Adult female survival
#> 866 0.6767919347 0.9972329      0.051  Adult female survival
#> 867 0.6627527504 0.9739464      0.008  Adult female survival
#> 868 0.6734056375 0.9840584      0.018  Adult female survival
#> 869 0.6563486361 0.9760971      0.012  Adult female survival
#> 870 0.6834248373 0.9888001      0.022  Adult female survival
#> 871 0.6511431905 0.9678314      0.006  Adult female survival
#> 872 0.6408376910 0.9641714      0.006  Adult female survival
#> 873 0.6621283144 0.9766726      0.010  Adult female survival
#> 874 0.6817748071 0.9914487      0.031  Adult female survival
#> 875 0.6577704924 0.9627421      0.002  Adult female survival
#> 876 0.6515727172 0.9638559      0.003  Adult female survival
#> 877 0.6582124774 0.9791388      0.013  Adult female survival
#> 878 0.6881187728 0.9883040      0.021  Adult female survival
#> 879 0.6582091137 0.9709627      0.004  Adult female survival
#> 880 0.6676610967 0.9805358      0.010  Adult female survival
#> 881 0.6864042895 0.9955305      0.042  Adult female survival
#> 882 0.6577090198 0.9714355      0.009  Adult female survival
#> 883 0.6449394687 0.9705025      0.007  Adult female survival
#> 884 0.6644643009 0.9788343      0.008  Adult female survival
#> 885 0.6909831361 0.9897023      0.024  Adult female survival
#> 886 0.6602339448 0.9659922      0.004  Adult female survival
#> 887 0.6427782696 0.9553078      0.003  Adult female survival
#> 888 0.6619596420 0.9771397      0.011  Adult female survival
#> 889 0.6875209350 0.9903587      0.027  Adult female survival
#> 890 0.6552832521 0.9693033      0.003  Adult female survival
#> 891 0.6725237977 0.9851960      0.018  Adult female survival
#> 892 0.6843157876 0.9915257      0.027  Adult female survival
#> 893 0.6525884332 0.9725420      0.009  Adult female survival
#> 894 0.6450892234 0.9706497      0.010  Adult female survival
#> 895 0.6536829341 0.9789506      0.014  Adult female survival
#> 896 0.6786919174 0.9885163      0.022  Adult female survival
#> 897 0.6538155493 0.9663520      0.007  Adult female survival
#> 898 0.6458362510 0.9637695      0.005  Adult female survival
#> 899 0.6634209248 0.9736430      0.005  Adult female survival
#> 900 0.6703185746 0.9902842      0.026  Adult female survival
#> 901 0.6947181598 0.9959137      0.054  Adult female survival
#> 902 0.6658093988 0.9795541      0.013  Adult female survival
#> 903 0.6793453241 0.9927221      0.033  Adult female survival
#> 904 0.6411340193 0.9677318      0.005  Adult female survival
#> 905 0.6567613809 0.9684318      0.005  Adult female survival
#> 906 0.6518106957 0.9781257      0.010  Adult female survival
#> 907 0.6780405059 0.9865346      0.021  Adult female survival
#> 908 0.6528363291 0.9654920      0.004  Adult female survival
#> 909 0.6499639406 0.9554999      0.004  Adult female survival
#> 
```

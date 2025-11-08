# caribouMetrics: Models and Metrics of Boreal Caribou Demography and Habitat Selection

caribouMetrics provides implementations of several models of Boreal
woodland caribou demography and habitat selection. A national
demographic model with density dependence and interannual variability
follows \[Johnson et. al. (2020)\](doi:10.1111/1365-2664.13637) with
modifications described in \[Dyson et al.
(2022)\](https://doi.org/10.1101/2022.06.01.494350). Demographic rates
vary with disturbance as estimated by \[Johnson et. al.
(2020)\](doi:10.1111/1365-2664.13637). The package also includes
Bayesian methods for integrating prior information from Johnson et al's
national analysis of demographic-disturbance relationships with
available local demographic data to reduce uncertainty in population
viability projections (\[Hughes et al.
2025\](https://doi.org/10.1016/j.ecoinf.2025.103095)). The national
model can be used to simulate example population trajectories, and
combined with a simple observation model and the Bayesian population
model to show how monitoring requirements depend on landscape condition
(\[Hughes et al. 2025\](https://doi.org/10.1016/j.ecoinf.2025.103095)).
Finally, caribouMetrics contains an implementation of \[Hornseth and
Rempel's (2016)\](https://doi.org/10.1139/cjz-2015-0101) Ontario boreal
caribou resource selection model described in \[Dyson et al.
(2022)\](https://doi.org/10.1101/2022.06.01.494350). Model
implementation is intended to be modular and flexible, allowing reuse of
components in a variety of contexts including projections of the
cumulative effects of disturbance and climate change \[(e.g. Stewart et
al. 2023)\](https://doi.org/10.1002/eap.2816).

## See also

Useful links:

- <https://landscitech.github.io/caribouMetrics>

- <https://github.com/LandSciTech/caribouMetrics>

## Author

**Maintainer**: Sarah Endicott <sarah.endicott@ec.gc.ca>
([ORCID](https://orcid.org/0000-0001-9644-5343))

Authors:

- Josie Hughes <josie.hughes@ec.gc.ca>
  ([ORCID](https://orcid.org/0000-0001-7875-9015))

- Yuko Shimoda <shimoday@gmail.com>

- Craig Simpkins <simpkinscraig063@gmail.com>

- Tati Michelleti <tati.micheletti@gmail.com> (Functions getCoefs and
  sampleRates are derived from code written by Tati Micheletti)
  \[copyright holder\]

- Eliot McIntire <eliot.mcintire@canada.ca>

Other contributors:

- Saralees Nadarajah (Author of truncdist package, which function rtrunc
  and qtrunc were modified from) \[copyright holder\]

- Frederick Novomestky (Author of truncdist package, which function
  rtrunc and qtrunc were modified from) \[copyright holder\]

- His Majesty the King in Right of Canada, as represented by the
  Minister of Environment and Climate Change \[copyright holder\]

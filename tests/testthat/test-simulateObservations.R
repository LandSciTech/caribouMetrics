test_that("default works", {
  scns <- getScenarioDefaults()
  expect_is(simulateObservations(scns,
                       freqStartsByYear = data.frame(Year = 2014:2023,
                                                     numStarts = 10),
                       cowCounts = data.frame(Year = 2014:2023,
                                              Count = 10,
                                              Class = "cow")),
            "list")
})

test_that("error messages are as expected", {
  scns <- getScenarioDefaults()
  expect_error(
    simulateObservations(scns,
                         freqStartsByYear = data.frame(Year = 2014:2023,
                                                       numStarts = 10),
                         cowCounts = data.frame(Year = 2014:2016,
                                                Count = 10,
                                                Class = "cow")),
    "Year is missing expected values")

  expect_error(
    simulateObservations(scns,
                         freqStartsByYear = data.frame(Year = 2014:2023,
                                                       numStarts = 10),
                         cowCounts = data.frame(Year = 2014:2023,
                                                Count = 10,
                                                Class = "cat")),
    "Class contains unexpected values")
})

test_that("multiple scenarios not allowed",{
  scns <- getScenarioDefaults(data.frame(iF = 1:2))
  expect_error(simulateObservations(scns,
                                 freqStartsByYear = data.frame(Year = 2014:2023,
                                                               numStarts = 10),
                                 cowCounts = data.frame(Year = 2014:2023,
                                                        Count = 10,
                                                        Class = "cow")),
            "must have length 1")
})

# TODO: add test for non-default popGrow table should use testPopGrowTable
# internally, do same in getPriors

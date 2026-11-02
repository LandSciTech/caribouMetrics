test_that("works", {
  disturb <- data.frame(Year = c(2000:2005, 2008:2010), Anthro = 1:9, Fire_excl_anthro = 2)

  year_dat <- data.frame(Year = 2000:2010)
  expect_warning(
    disturb_filled <- fillDisturbance(disturb, year_dat), "missing disturbance"
  )
  
  expect_equal(disturb_filled %>% filter(Year == 2006) %>% pull(Anthro), 6)
  expect_equal(disturb_filled %>% filter(Year == 2006) %>% pull(Fire_excl_anthro), 0)
  
  expect_warning(
    disturb_null <- fillDisturbance(disturb, data.frame(Year = 2020:2025)), 
    "None of the disturbance data"
  )
  
})

test_that("works with multipop", {
  disturb <- data.frame(Year = c(2000:2005, 2008:2010), Anthro = 1:9, Fire_excl_anthro = 2,
                        PopulationName = "A") %>% 
    bind_rows(
      data.frame(Year = c(2020:2025, 2028:2030), Anthro = 1:9, Fire_excl_anthro = 2,
                 PopulationName = "B")
    )
  
  year_dat <- bind_rows(
    data.frame(Year = 2000:2010, PopulationName = "A"),
    data.frame(Year = 2020:2030, PopulationName = "B")
    )
  expect_warning(
    disturb_filled <- fillDisturbance(disturb, year_dat), "missing disturbance"
  )
  
  expect_equal(disturb_filled %>% filter(Year == 2006, PopulationName == "A") %>% pull(Anthro), 6)
  expect_equal(disturb_filled %>% filter(Year == 2026, PopulationName == "A") %>% pull(Anthro), integer(0))
  expect_equal(disturb_filled %>% filter(Year == 2026, PopulationName == "B") %>% pull(Anthro), 6)
  expect_equal(disturb_filled %>% filter(Year == 2006, PopulationName == "B") %>% pull(Anthro), integer(0))
  
  expect_equal(disturb_filled %>% filter(Year == 2006) %>% pull(Fire_excl_anthro), 0)

  
})

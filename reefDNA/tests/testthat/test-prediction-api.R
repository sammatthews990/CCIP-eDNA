test_that("operational predictors are transformed consistently", {
  scale_spec <- list(
    edna_prop_model = c(center = 0.5, scale = 0.25),
    distance_model = c(center = log1p(500 / 200), scale = 0.5),
    lag_model = c(center = log1p(30), scale = 1)
  )
  input <- data.frame(
    perc_pos = c(50, 75),
    distance_m = c(500, 200),
    lag_days = c(30, 60)
  )
  result <- reefDNA:::.prepare_operational_predictors(input, scale_spec)

  expect_equal(result$edna_prop_z, c(0, 1))
  expect_equal(result$distance_z[1], 0)
  expect_equal(result$lag_z[1], 0)
  expect_equal(result$bottom_time, c(216, 216))
  expect_true(all(is.finite(result$log_effort)))
  expect_equal(result$Reef, c("new_reef_1", "new_reef_2"))
})

test_that("operational predictor validation catches invalid values", {
  scale_spec <- list(
    edna_prop_model = c(center = 0.5, scale = 0.25),
    distance_model = c(center = 1, scale = 0.5),
    lag_model = c(center = 1, scale = 0.5)
  )
  expect_error(
    reefDNA:::.prepare_operational_predictors(
      data.frame(perc_pos = 101, distance_m = 10, lag_days = 1),
      scale_spec
    ),
    "between 0 and 100"
  )
  expect_error(
    reefDNA:::.prepare_operational_predictors(
      data.frame(perc_pos = 50, distance_m = -1, lag_days = 1),
      scale_spec
    ),
    "non-negative"
  )
})

test_that("site-visit grouping columns are constructed or retained", {
  scale_spec <- list(
    edna_prop_model = c(center = 0.5, scale = 0.2),
    distance_model = c(center = 1, scale = 0.5),
    lag_model = c(center = 3, scale = 1)
  )
  prepared <- reefDNA:::.prepare_operational_predictors(
    data.frame(
      perc_pos = c(50, 60), distance_m = c(200, 500),
      lag_days = c(30, 60), Reef = "Example Reef",
      cull_site_name = c("A", "B"), edna_campaign_id = "campaign_1"
    ),
    scale_spec, default_effort = 480
  )
  expect_equal(prepared$site_id, c("Example Reef__A", "Example Reef__B"))
  expect_equal(prepared$edna_campaign_id, rep("campaign_1", 2))
  expect_equal(prepared$bottom_time, rep(480, 2))
})

test_that("sampling design is encoded for design-aware BRT models", {
  scale_spec <- list(
    edna_prop_model = c(center = 0.5, scale = 0.2),
    distance_model = c(center = 1, scale = 0.5),
    lag_model = c(center = 3, scale = 1)
  )
  prepared <- reefDNA:::.prepare_operational_predictors(
    data.frame(
      perc_pos = rep(50, 3), distance_m = rep(200, 3), lag_days = rep(30, 3),
      sampling_design = c("3x12", "4 sites x 6 reps", "other_or_mixed")
    ),
    scale_spec
  )
  expect_equal(prepared$design_4x6, c(0L, 1L, 0L))
  expect_equal(prepared$design_other_or_mixed, c(0L, 0L, 1L))
})

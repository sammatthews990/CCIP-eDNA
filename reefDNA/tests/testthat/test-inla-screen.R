test_that("event benchmark selects the most recent event without duplicating culls", {
  edna_raw <- data.frame(
    Date = as.Date(c("2025-01-01", "2025-01-01", "2025-01-08", "2025-01-08")),
    ReefName = "Eyrie",
    Site_name = c("old", "old", "recent", "recent"),
    Lat = c(-14.7000, -14.7001, -14.7010, -14.7011),
    Long = c(145.3700, 145.3701, 145.3710, 145.3711),
    LOD_sample_positive = c(0, 0, 1, 1),
    Conc_mean = c(0, 0, 10, 20)
  )
  cull_raw <- data.frame(
    CrownOfThornsStarfishCullDiveId = "cull_1",
    SurveyDate = as.Date("2025-01-10"),
    ReefName = "Eyrie",
    ReefLabel = "14-118",
    CullSiteName = "EYR_14-118_1",
    Latitude = -14.7005,
    Longitude = 145.3705,
    Bottomtime = 20,
    Depth = 5,
    Cohort1 = 1,
    Cohort2 = 0,
    Cohort3 = 0,
    Cohort4 = 0
  )

  edna <- prepare_inla_edna(edna_raw)
  culls <- prepare_inla_culls(cull_raw)
  links <- build_inla_space_time_links(edna, culls, max_distance_m = 2000)
  benchmark <- summarise_inla_event_benchmark(links, edna, culls)

  expect_equal(nrow(benchmark), 1)
  expect_equal(benchmark$edna_site_name, "recent")
  expect_equal(benchmark$edna_lag_days, 2)
  expect_equal(benchmark$edna_prop_positive, 1)
  expect_equal(benchmark$edna_conc_mean, 15)
  expect_equal(benchmark$edna_n_replicates, 2)
  expect_equal(anyDuplicated(benchmark$cull_id), 0)
})

test_that("additive candidate generator returns every predictor subset", {
  candidates <- build_inla_additive_candidates()
  expect_length(candidates, 16)
  expect_equal(candidates$null, character())
  expect_true(any(lengths(candidates) == 4))
  expect_true(all(vapply(candidates, function(x) all(x %in% c(
    "edna_prop_z", "edna_conc_z", "distance_z", "lag_z"
  )), logical(1))))
})

test_that("interaction candidates obey marginality", {
  candidates <- build_inla_interaction_candidates()
  expect_length(candidates, 7)
  expect_true(all(vapply(candidates, function(x) {
    all(c("edna_prop_z", "edna_conc_z", "distance_z", "lag_z") %in% x)
  }, logical(1))))
})

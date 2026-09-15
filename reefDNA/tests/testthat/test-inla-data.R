test_that("site identifiers are normalized consistently", {
  expect_equal(
    normalize_cull_site_id(c("EYR_14-118_1", "(eyr 14-118 1")),
    c("EYR141181", "EYR141181")
  )
})

test_that("raw eDNA links are aggregated to one row per cull dive", {
  edna_raw <- data.frame(
    Date = as.Date(c("2025-01-01", "2025-01-01", "2025-01-01")),
    ReefName = c("Eyrie", "Eyrie", "Other"),
    Site_name = c("A", "A", "B"),
    Lat = c(-14.70, -14.7001, -14.70),
    Long = c(145.37, 145.3701, 145.37),
    LOD_sample_positive = c(1, 0, 1),
    Conc_mean = c(10, 0, 5)
  )
  cull_raw <- data.frame(
    CrownOfThornsStarfishCullDiveId = "cull_1",
    SurveyDate = as.Date("2025-01-10"),
    ReefName = "Eyrie",
    ReefLabel = "14-118",
    CullSiteName = "EYR_14-118_1",
    Latitude = -14.70005,
    Longitude = 145.37005,
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
  benchmark <- summarise_inla_benchmark(links, edna, culls)

  expect_equal(nrow(links), 2)
  expect_equal(nrow(benchmark), 1)
  expect_equal(benchmark$n_edna, 2)
  expect_equal(benchmark$edna_prop_positive, 0.5)
  expect_equal(benchmark$cots_count, 1)
  expect_equal(benchmark$cpue, 0.05)
  expect_equal(anyDuplicated(benchmark$cull_id), 0)
})

test_that("metric link construction rejects longitude-latitude geometry", {
  point <- sf::st_as_sf(
    data.frame(
      edna_id = "e1", event_id = "event", Reef = "Eyrie",
      date_edna = as.Date("2025-01-01"), x = 145.37, y = -14.70
    ),
    coords = c("x", "y"), crs = 4326
  )
  cull <- sf::st_as_sf(
    data.frame(
      cull_id = "c1", Reef = "Eyrie", date_cull = as.Date("2025-01-02"),
      x = 145.37, y = -14.70
    ),
    coords = c("x", "y"), crs = 4326
  )
  expect_error(build_inla_space_time_links(point, cull), "projected CRS")
})

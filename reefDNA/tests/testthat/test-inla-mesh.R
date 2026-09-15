test_that("SPDE mesh helper returns auditable dimensions", {
  skip_if_not_installed("INLA")
  points <- sf::st_as_sf(
    data.frame(x = c(0, 10000, 0, 10000), y = c(0, 0, 10000, 10000)),
    coords = c("x", "y"), crs = 3112
  )
  mesh <- build_inla_spde_mesh(
    points, cutoff_m = 1000,
    max_edge_m = c(5000, 10000), offset_m = c(5000, 10000)
  )
  summary <- summarise_inla_mesh(mesh)
  expect_gt(summary$n_vertices, 4)
  expect_gt(summary$n_triangles, 1)
  expect_gt(summary$median_edge_m, 0)
})

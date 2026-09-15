#' Build a projected SPDE mesh for COTS CPUE modelling
#'
#' @param locations An `sf` point object in a projected metre-based CRS.
#' @param cutoff_m Minimum separation between retained mesh input locations.
#' @param max_edge_m Inner and outer maximum triangle edge lengths.
#' @param offset_m Inner and outer mesh extension distances.
#' @return An INLA two-dimensional mesh.
#' @export
build_inla_spde_mesh <- function(
    locations, cutoff_m = 5000, max_edge_m = c(20000, 80000),
    offset_m = c(20000, 100000)) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("The INLA package is required to build the SPDE mesh", call. = FALSE)
  }
  if (!inherits(locations, "sf")) stop("locations must be an sf object", call. = FALSE)
  if (sf::st_is_longlat(locations)) stop("locations must use a projected metre-based CRS", call. = FALSE)

  coords <- sf::st_coordinates(locations)
  coords <- unique(coords[, c("X", "Y"), drop = FALSE])
  if (nrow(coords) < 3L) stop("At least three unique locations are required", call. = FALSE)

  INLA::inla.mesh.2d(
    loc = coords,
    cutoff = cutoff_m,
    max.edge = max_edge_m,
    offset = offset_m
  )
}

#' Summarize an INLA SPDE mesh
#'
#' @param mesh An INLA mesh.
#' @return A one-row tibble of mesh dimensions and edge-length summaries.
#' @export
summarise_inla_mesh <- function(mesh) {
  triangles <- mesh$graph$tv
  loc <- mesh$loc[, 1:2, drop = FALSE]
  edges <- rbind(
    triangles[, c(1, 2), drop = FALSE],
    triangles[, c(2, 3), drop = FALSE],
    triangles[, c(1, 3), drop = FALSE]
  )
  edges <- unique(t(apply(edges, 1, sort)))
  edge_length <- sqrt(rowSums((loc[edges[, 1], , drop = FALSE] - loc[edges[, 2], , drop = FALSE])^2))

  tibble::tibble(
    n_vertices = mesh$n,
    n_triangles = nrow(triangles),
    min_edge_m = min(edge_length),
    median_edge_m = stats::median(edge_length),
    max_edge_m = max(edge_length)
  )
}

library(sf)

cat("Reading sector shapefile...\n")
sector_sf <- st_read("data/SectorShapefile", quiet = TRUE)
print(st_crs(sector_sf))
print(st_bbox(sector_sf))

cat("\nReading RRAP canonical GPKG...\n")
gbr_reefs <- st_read("data/rrap_canonical_2025-03-20-T15-18-17.gpkg", quiet = TRUE)
print(st_crs(gbr_reefs))
print(st_bbox(gbr_reefs))

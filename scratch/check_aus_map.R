library(rnaturalearth)
library(sf)

aus <- ne_countries(country = "australia", scale = "medium", returnclass = "sf")
print(st_bbox(aus))

# Check maps package
cat("maps package: ", requireNamespace("maps", quietly = TRUE), "\n")
if (requireNamespace("maps", quietly = TRUE)) {
    aus_map <- st_as_sf(maps::map("world", "Australia", plot = FALSE, fill = TRUE))
    print(st_bbox(aus_map))
}

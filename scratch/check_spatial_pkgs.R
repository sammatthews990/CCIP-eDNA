library(sf)
library(ggplot2)

cat("Checking available spatial packages and maps...\n")

packages <- c("ozmaps", "rnaturalearth", "rnaturalearthdata", "ggspatial", "ggrepel", "cowplot")
for (p in packages) {
    cat(p, ": ", requireNamespace(p, quietly = TRUE), "\n")
}

# Check if gpkg can be read
gpkg_path <- "data/Eotr_CotsCullSites_2025_11_19_1_58_PM.gpkg"
if (file.exists(gpkg_path)) {
    layers <- st_layers(gpkg_path)
    cat("GPKG layers:\n")
    print(layers)
}

library(sf)

gpkg <- st_read("data/Eotr_CotsCullSites_2025_11_19_1_58_PM.gpkg")
print(head(gpkg))
print(st_bbox(gpkg))

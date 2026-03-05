---
title: Spatial Vector Analysis and Geoprocessing
slug: /domain-knowledge/geospatial/use-cases/spatial-vector-analysis
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Use Case 2: Spatial Vector Analysis and Geoprocessing**

This analysis works with vector data — points, lines, and polygons with attribute tables — to perform spatial queries, overlays, and statistical analysis. Examples include analyzing infrastructure accessibility, calculating population exposure to hazards, or aggregating environmental indicators by administrative region.

## **Workflow Overview**

```
Data collection → Spatial operations → Analysis → Visualization
Vector layers (GeoPackage/Shapefile) + Tabular data → Spatial joins / overlays → Statistics & maps
```

## **Key Concepts**

### **Spatial Operations**

Core geoprocessing operations include:
- **Spatial join** — Transfer attributes from one layer to another based on spatial relationship (intersects, contains, nearest)
- **Buffer** — Create zones of a specified distance around features
- **Overlay (intersection, union, difference)** — Combine two polygon layers to produce new features
- **Dissolve** — Merge features based on a shared attribute
- **Clip** — Cut features to the boundary of another layer

### **Spatial Statistics**

Beyond simple overlay:
- **Point pattern analysis** — Kernel density estimation, nearest-neighbor analysis
- **Spatial autocorrelation** — Moran's I, Getis-Ord Gi* for hot spot detection
- **Zonal statistics** — Summarize raster values within vector zones (e.g., mean elevation per municipality)

## **Popular Tools**

- **QGIS** — Full-featured desktop GIS with built-in geoprocessing
- **PostGIS** — Spatial extension for PostgreSQL databases
- **GeoPandas** (Python) — Spatial operations on DataFrames
- **sf** (R) — Simple features for R
- **DuckDB Spatial** — Fast SQL-based spatial analytics
- **ArcGIS Pro** — Commercial GIS platform

## **Code Example: Spatial Analysis with Python**

<details>
<summary>**Click to expand Python geoprocessing workflow**</summary>

```python
#!/usr/bin/env python3
"""Spatial vector analysis: accessibility to public transport stops."""

import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt

# Step 1: Load data
municipalities = gpd.read_file("municipalities.gpkg")
transport_stops = gpd.read_file("public_transport_stops.gpkg")

# Ensure same CRS (project to metric CRS for distance calculations)
municipalities = municipalities.to_crs(epsg=25832)  # UTM 32N
transport_stops = transport_stops.to_crs(epsg=25832)

# Step 2: Create 500m buffer around each transport stop
buffers = transport_stops.copy()
buffers["geometry"] = transport_stops.buffer(500)

# Step 3: Dissolve all buffers into a single coverage area
coverage = buffers.dissolve()

# Step 4: Calculate area served per municipality
municipalities["total_area_m2"] = municipalities.area
intersection = gpd.overlay(municipalities, coverage, how="intersection")
intersection["covered_area_m2"] = intersection.area

# Aggregate covered area back to municipalities
covered_by_muni = intersection.groupby("municipality_name")["covered_area_m2"].sum().reset_index()
municipalities = municipalities.merge(covered_by_muni, on="municipality_name", how="left")
municipalities["covered_area_m2"] = municipalities["covered_area_m2"].fillna(0)
municipalities["coverage_pct"] = (municipalities["covered_area_m2"]
                                   / municipalities["total_area_m2"] * 100)

# Step 5: Count stops per municipality
stops_per_muni = gpd.sjoin(transport_stops, municipalities, predicate="within")
stop_counts = stops_per_muni.groupby("municipality_name").size().reset_index(name="n_stops")
municipalities = municipalities.merge(stop_counts, on="municipality_name", how="left")
municipalities["n_stops"] = municipalities["n_stops"].fillna(0).astype(int)

# Step 6: Visualize
fig, axes = plt.subplots(1, 2, figsize=(16, 8))
municipalities.plot(column="coverage_pct", cmap="RdYlGn", legend=True,
                    ax=axes[0], edgecolor="gray", linewidth=0.3)
axes[0].set_title("Public Transport Coverage (%)")
municipalities.plot(column="n_stops", cmap="YlOrRd", legend=True,
                    ax=axes[1], edgecolor="gray", linewidth=0.3)
axes[1].set_title("Number of Transport Stops")

plt.tight_layout()
plt.savefig("transport_accessibility.pdf", dpi=300)
plt.savefig("transport_accessibility.png", dpi=300)
print("Analysis complete!")
```

</details>

## **Code Example: Spatial Analysis with R**

<details>
<summary>**Click to expand R geoprocessing workflow**</summary>

```r
#!/usr/bin/env Rscript
# Spatial vector analysis: point-in-polygon and zonal statistics

library(sf)
library(dplyr)
library(ggplot2)

# Step 1: Load data
municipalities <- st_read("municipalities.gpkg") %>%
  st_transform(25832)  # UTM 32N

transport_stops <- st_read("public_transport_stops.gpkg") %>%
  st_transform(25832)

# Step 2: Count points in polygons
stops_in_muni <- st_join(transport_stops, municipalities, join = st_within)
stop_counts <- stops_in_muni %>%
  st_drop_geometry() %>%
  count(municipality_name, name = "n_stops")

municipalities <- municipalities %>%
  left_join(stop_counts, by = "municipality_name") %>%
  mutate(n_stops = replace_na(n_stops, 0))

# Step 3: Buffer and coverage analysis
buffers <- st_buffer(transport_stops, dist = 500) %>%
  st_union()

municipalities$total_area <- st_area(municipalities)
covered <- st_intersection(municipalities, buffers)
covered$covered_area <- st_area(covered)

coverage_df <- covered %>%
  st_drop_geometry() %>%
  group_by(municipality_name) %>%
  summarise(covered_area = sum(covered_area))

municipalities <- municipalities %>%
  left_join(coverage_df, by = "municipality_name") %>%
  mutate(
    covered_area = replace_na(as.numeric(covered_area), 0),
    coverage_pct = covered_area / as.numeric(total_area) * 100
  )

# Step 4: Visualize
ggplot(municipalities) +
  geom_sf(aes(fill = coverage_pct), color = "gray50", linewidth = 0.2) +
  scale_fill_viridis_c(name = "Coverage (%)") +
  theme_minimal() +
  ggtitle("Public Transport Coverage by Municipality")

ggsave("transport_coverage.pdf", width = 10, height = 8, dpi = 300)
ggsave("transport_coverage.png", width = 10, height = 8, dpi = 300)
```

</details>

## **Expected Outputs**

- **Attribute-enriched vector layers** — GeoPackage or Shapefile with computed spatial statistics
- **Thematic maps** — Choropleth maps showing spatial patterns (PDF + PNG)
- **Summary tables** — CSV with aggregated statistics per spatial unit

## **Computational Requirements**

| Task | CPU Cores | RAM | Storage | Time |
|------|-----------|-----|---------|------|
| Municipal-level analysis (single country) | 2-4 | 4-8 GB | `<1 GB` | Minutes |
| Large-scale overlay (millions of features) | 8-16 | 16-64 GB | 1-10 GB | 10-60 min |
| Database-backed analysis (PostGIS) | 4-8 | 16-32 GB | Varies | Varies |

## **Common Issues & Troubleshooting**

:::warning Common Problems

**Topology errors**
- Invalid geometries (self-intersections, unclosed polygons) cause overlay operations to fail
- Fix with `shapely.make_valid()` (Python) or `st_make_valid()` (R)

**CRS mismatch**
- Always reproject all layers to the same CRS before spatial operations
- Use metric CRS (e.g., UTM) for distance and area calculations, not WGS 84

**Memory issues with large datasets**
- Use spatial indexing (GeoPackage/PostGIS have built-in spatial indices)
- Process data in tiles or chunks
- Consider DuckDB Spatial or PostGIS for datasets too large for in-memory processing

:::

## **Key Considerations**

- **Always validate geometries** before running overlay operations
- **Use GeoPackage over Shapefile** for new projects (no size limits, supports NULL values, single file)
- **Document your CRS choices** — metric projections for analysis, geographic CRS for data sharing
- **Spatial joins can be many-to-many** — decide how to handle features that match multiple targets

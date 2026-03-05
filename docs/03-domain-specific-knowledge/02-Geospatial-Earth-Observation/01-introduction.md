---
title: Introduction to Geospatial & Earth Observation
slug: /domain-knowledge/geospatial/introduction
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **What is Geospatial & Earth Observation Data Science?**

Geospatial data science deals with data that has a spatial component — information tied to locations on the Earth's surface. Earth observation (EO) extends this by leveraging satellite, airborne, and ground-based sensors to systematically monitor the planet. Together, these fields combine geography, remote sensing, computer science, and statistics to extract knowledge from spatially referenced datasets.

## **The Need for Geospatial Data Management**

Geospatial and EO data present unique data management challenges:

- **Volume** — A single Sentinel-2 satellite scene is ~1 GB; global archives reach petabytes
- **Coordinate reference systems (CRS)** — Data must be consistently projected and georeferenced
- **Temporal dimension** — Many analyses require time series spanning years or decades
- **Heterogeneous sources** — Combining raster imagery, vector boundaries, point measurements, and tabular attributes
- **Resolution trade-offs** — Spatial, temporal, and spectral resolution vary across sensors and must be harmonized

## **Key Application Areas**

Geospatial and EO methods are applied across many disciplines:

- **Environmental monitoring** — Land use/land cover change, deforestation, glacier retreat
- **Climate science** — Temperature trends, precipitation patterns, carbon flux estimation
- **Urban planning** — Infrastructure mapping, population density, mobility analysis
- **Agriculture** — Crop monitoring, yield prediction, precision farming
- **Disaster management** — Flood mapping, wildfire detection, earthquake damage assessment
- **Ecology and biodiversity** — Habitat mapping, species distribution modeling

## **Common Tools and Software**

### **GIS Platforms**
- **QGIS** — Open-source desktop GIS for visualization and analysis
- **ArcGIS** — Commercial GIS platform (ESRI), widely used in industry and academia
- **Google Earth Engine** — Cloud-based platform for planetary-scale EO analysis

### **Programming Libraries**
- **GDAL/OGR** — The foundational library for reading/writing raster and vector geospatial formats
- **Rasterio** (Python) — Pythonic interface for raster data
- **GeoPandas** (Python) — Extends pandas with spatial operations
- **sf** (R) — Simple features for R, modern spatial data handling
- **terra** (R) — Raster and vector analysis
- **xarray** (Python) — Multi-dimensional labeled arrays, ideal for NetCDF/climate data

### **Remote Sensing**
- **SNAP** (ESA) — Sentinel Application Platform for satellite data processing
- **Orfeo ToolBox** — Open-source image processing for remote sensing
- **ENVI** — Commercial remote sensing analysis software

## **Infrastructure and Support**

Researchers working with geospatial data at TUM can access:

- **LRZ Linux Cluster** — For large-scale raster processing and modeling
- **Terrabyte** — LRZ platform specifically designed for Earth observation data science
- **Copernicus Data Space** — Free access to Sentinel satellite data
- **Google Earth Engine** — Cloud computing for global-scale analysis

## **Getting Started**

If you're new to geospatial data science:

1. **Explore our resources** — Check out the [Data Types](/domain-knowledge/geospatial/data-types) and [Use Cases](/domain-knowledge/geospatial/use-cases)
2. **Access infrastructure** — Apply for Terrabyte access at [LRZ](https://www.lrz.de/)
3. **Get support** — Contact [rdhub@mdsi.tum.de](mailto:rdhub@mdsi.tum.de) for data management guidance

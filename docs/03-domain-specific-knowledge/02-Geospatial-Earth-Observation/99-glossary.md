---
title: Glossary
slug: /domain-knowledge/geospatial/glossary
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Geospatial & Earth Observation Glossary**

This glossary provides definitions for common terms used in geospatial data science and Earth observation. Understanding these terms is essential for working with spatial data and communicating effectively with collaborators.

---

## **Spatial Data Fundamentals**

### **Raster Data**
Data represented as a regular grid of cells (pixels), where each cell holds a numerical value. Used for continuous phenomena like elevation, temperature, or satellite imagery.

### **Vector Data**
Data represented as discrete geometric features — points, lines, and polygons — with associated attribute tables. Used for discrete features like buildings, roads, and administrative boundaries.

### **Coordinate Reference System (CRS)**
A system that defines how coordinates (numbers) relate to positions on the Earth's surface. Includes a datum (model of the Earth's shape), a projection (how the 3D surface is flattened to 2D), and a coordinate system (units and axes).

### **EPSG Code**
A standardized numerical identifier for coordinate reference systems, maintained by the IOGP. Example: EPSG:4326 = WGS 84 (latitude/longitude).

### **Projection**
A mathematical transformation that converts 3D coordinates on the Earth's surface to 2D coordinates on a flat map. All projections introduce some distortion (area, shape, distance, or direction).

### **Datum**
A mathematical model of the Earth's shape used as a reference surface. Common datums: WGS 84 (global, GPS), ETRS89 (Europe).

### **Georeferencing**
The process of assigning real-world coordinates to data (an image, a map, a dataset) so it can be placed correctly on the Earth's surface.

---

## **Remote Sensing Terms**

### **Remote Sensing**
The acquisition of information about objects or phenomena from a distance, typically using sensors on satellites or aircraft.

### **Multispectral Imagery**
Imagery captured in multiple discrete wavelength bands (typically 4-12 bands), such as blue, green, red, and near-infrared.

### **Hyperspectral Imagery**
Imagery captured in many narrow, contiguous wavelength bands (typically 100-300+), providing a near-continuous spectrum for each pixel.

### **Spatial Resolution**
The size of a single pixel on the ground. A 10 m resolution means each pixel covers a 10 × 10 m area. Finer resolution = more detail but larger file sizes.

### **Temporal Resolution**
How frequently a sensor revisits the same location. Sentinel-2 has a 5-day revisit time at the equator.

### **Spectral Resolution**
The number and width of wavelength bands a sensor records. More and narrower bands = finer spectral detail.

### **Radiometric Resolution**
The sensitivity of a sensor to differences in energy, typically expressed in bits (e.g., 12-bit = 4096 levels).

### **NDVI (Normalized Difference Vegetation Index)**
A spectral index calculated as (NIR - Red) / (NIR + Red) that indicates vegetation health and density. Values range from -1 to +1, with healthy vegetation typically above 0.3.

### **Atmospheric Correction**
The process of removing atmospheric effects (scattering, absorption) from satellite imagery to obtain surface reflectance values.

### **Orthorectification**
The correction of geometric distortions in aerial or satellite imagery caused by terrain relief and sensor geometry, producing a geometrically accurate map-like product.

---

## **GIS and Analysis Terms**

### **GIS (Geographic Information System)**
A system for capturing, storing, analyzing, and visualizing spatial data and its associated attributes.

### **Spatial Join**
An operation that transfers attributes from one layer to another based on their spatial relationship (e.g., points falling within polygons).

### **Buffer**
A zone of specified distance created around a geographic feature. Used for proximity analysis (e.g., 500 m buffer around a school).

### **Overlay**
A set of operations (intersection, union, difference) that combine two spatial layers to produce new features.

### **Geocoding**
The process of converting addresses or place names into geographic coordinates (latitude/longitude).

### **Spatial Autocorrelation**
The degree to which values at nearby locations are similar. Measured by statistics like Moran's I.

### **Zonal Statistics**
The calculation of statistics (mean, sum, count, etc.) for raster cell values within vector polygon zones.

---

## **Data and Format Terms**

### **GeoTIFF**
A TIFF image format extended with embedded georeferencing metadata (CRS, extent, resolution).

### **Cloud Optimized GeoTIFF (COG)**
A GeoTIFF with internal tiling and overviews optimized for efficient partial reads over HTTP.

### **Shapefile**
A legacy multi-file vector format (.shp, .shx, .dbf, .prj) from ESRI. Being replaced by GeoPackage.

### **GeoPackage**
An open, SQLite-based container for vector and raster geospatial data. The modern replacement for Shapefiles.

### **NetCDF**
Network Common Data Form — a self-describing format for array-oriented scientific data, widely used in climate science.

### **STAC (SpatioTemporal Asset Catalog)**
A JSON-based specification for describing and cataloging geospatial assets for search and discovery.

### **WMS / WFS / WCS**
OGC web service standards for serving map images (WMS), vector features (WFS), and raster coverages (WCS) over the internet.

### **Point Cloud**
A collection of 3D points (X, Y, Z) representing the external surface of objects, typically acquired by LiDAR or photogrammetry.

### **DEM (Digital Elevation Model)**
A raster representation of terrain elevation. Related: DSM (Digital Surface Model) includes buildings and vegetation; DTM (Digital Terrain Model) represents bare earth.

---

**Need clarification on a term?** Contact us at [rdhub@mdsi.tum.de](mailto:rdhub@mdsi.tum.de) or explore our [Use Cases](/domain-knowledge/geospatial/use-cases) for practical examples.

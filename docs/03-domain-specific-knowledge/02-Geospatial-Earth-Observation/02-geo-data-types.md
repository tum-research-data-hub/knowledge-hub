---
title: Geospatial Data Types & Formats
slug: /domain-knowledge/geospatial/data-types
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Geospatial Data Types & File Formats**

Geospatial data comes in two fundamental categories — **raster** (grid-based) and **vector** (geometry-based) — plus several specialized formats for point clouds, time series, and tabular spatial data. Understanding these formats is essential for choosing the right tools and workflows.

---

## **Raster Data**

Raster data represents the world as a regular grid of cells (pixels), where each cell holds a value. This is the primary format for satellite imagery, elevation models, and continuous surfaces.

### **GeoTIFF (.tif, .tiff)**
The most widely used raster format. A standard TIFF image extended with georeferencing metadata (CRS, extent, resolution) embedded in the file header.
- **Use:** Satellite imagery, elevation models, land cover maps
- **Tools:** GDAL, Rasterio, QGIS, Google Earth Engine
- **Spec:** [OGC GeoTIFF Standard](https://www.ogc.org/standard/geotiff/)

### **Cloud Optimized GeoTIFF (COG)**
A GeoTIFF with internal tiling and overviews that allows efficient partial reads over HTTP, enabling cloud-native workflows without downloading entire files.
- **Use:** Cloud-hosted imagery, web map services
- **Tools:** GDAL, Rasterio, STAC catalogs
- **Spec:** [COG Specification](https://www.cogeo.org/)

### **NetCDF (.nc)**
Network Common Data Form — a self-describing, machine-independent format for array-oriented scientific data. Widely used in climate and atmospheric science.
- **Use:** Climate model output, reanalysis data, oceanographic measurements
- **Tools:** xarray (Python), ncdf4 (R), CDO, NCO
- **Spec:** [Unidata NetCDF](https://www.unidata.ucar.edu/software/netcdf/)

### **HDF5 (.h5, .hdf5) / HDF-EOS**
Hierarchical Data Format version 5 — a flexible format for large, complex datasets with internal groups and metadata. HDF-EOS is NASA's extension for Earth observation data.
- **Use:** Satellite sensor data (e.g., MODIS, Landsat Collection 2), multi-dimensional arrays
- **Tools:** h5py (Python), GDAL, HDFView
- **Spec:** [HDF Group](https://www.hdfgroup.org/solutions/hdf5/)

### **JPEG2000 (.jp2)**
A wavelet-based compressed image format used by the European Space Agency for Sentinel-2 imagery.
- **Use:** Sentinel-2 Level-1C and Level-2A products
- **Tools:** GDAL, SNAP, Rasterio
- **Spec:** [ISO/IEC 15444](https://www.iso.org/standard/78321.html)

### **Zarr (.zarr)**
A cloud-native, chunked, compressed array format designed for parallel read/write access. Increasingly used as a modern alternative to NetCDF/HDF5.
- **Use:** Analysis-ready climate data, large raster time series
- **Tools:** xarray + Zarr (Python), Dask
- **Spec:** [Zarr Specification](https://zarr.readthedocs.io/)

---

## **Vector Data**

Vector data represents geographic features as points, lines, and polygons with associated attribute tables.

### **GeoPackage (.gpkg)**
An open, SQLite-based format that can store multiple vector layers, raster tiles, and attribute data in a single file. The modern replacement for Shapefiles.
- **Use:** Administrative boundaries, infrastructure networks, point observations
- **Tools:** QGIS, GDAL/OGR, GeoPandas, sf (R)
- **Spec:** [OGC GeoPackage](https://www.geopackage.org/)

### **Shapefile (.shp + .shx + .dbf + .prj)**
The legacy vector format from ESRI, still widely used despite limitations (2 GB size limit, 10-character field names, no NULL values). Always consists of multiple associated files.
- **Use:** Legacy datasets, interoperability with ArcGIS workflows
- **Tools:** QGIS, GDAL/OGR, GeoPandas, sf (R)
- **Note:** Consider migrating to GeoPackage or GeoJSON for new projects

### **GeoJSON (.geojson)**
A lightweight, text-based vector format using JSON syntax. Human-readable and web-friendly, but not suitable for large datasets.
- **Use:** Web mapping, API responses, small to medium datasets
- **Tools:** Any JSON parser, GeoPandas, Leaflet, Mapbox
- **Spec:** [RFC 7946](https://tools.ietf.org/html/rfc7946)

### **KML / KMZ (.kml, .kmz)**
Keyhole Markup Language — an XML-based format developed for Google Earth. KMZ is a zipped KML with embedded resources.
- **Use:** Visualization in Google Earth, simple data sharing
- **Tools:** Google Earth, QGIS, GDAL

### **GeoParquet (.parquet)**
A columnar storage format combining Apache Parquet's efficiency with geospatial metadata. Designed for large-scale analytical queries.
- **Use:** Big data analytics, cloud-native spatial data processing
- **Tools:** GeoPandas, DuckDB, Apache Spark
- **Spec:** [GeoParquet](https://geoparquet.org/)

---

## **Point Cloud Data**

### **LAS / LAZ (.las, .laz)**
The standard format for LiDAR point cloud data. LAZ is the compressed version. Each point contains XYZ coordinates plus attributes like intensity, classification, and return number.
- **Use:** Terrain modeling, forestry, urban 3D mapping
- **Tools:** PDAL, CloudCompare, LAStools, lidR (R)
- **Spec:** [ASPRS LAS Specification](https://www.asprs.org/divisions-committees/lidar-division/laser-las-file-format-exchange-activities)

---

## **Metadata and Catalog Formats**

### **STAC (SpatioTemporal Asset Catalog)**
A JSON-based specification for describing geospatial assets (imagery, point clouds, etc.) to make them searchable and discoverable.
- **Use:** Satellite data catalogs, data discovery
- **Tools:** pystac (Python), rstac (R), STAC Browser
- **Spec:** [STAC Spec](https://stacspec.org/)

### **ISO 19115 / INSPIRE**
International metadata standards for geographic information, widely used in European spatial data infrastructures.
- **Use:** Official metadata records, INSPIRE-compliant data portals
- **Spec:** [ISO 19115](https://www.iso.org/standard/53798.html)

---

## **Format Selection Guide**

| Use Case | Recommended Format |
|----------|-------------------|
| Satellite imagery (local) | GeoTIFF |
| Satellite imagery (cloud) | Cloud Optimized GeoTIFF (COG) |
| Climate / atmospheric data | NetCDF or Zarr |
| Vector features (new projects) | GeoPackage |
| Vector features (web) | GeoJSON |
| Large-scale vector analytics | GeoParquet |
| LiDAR point clouds | LAS/LAZ |
| Data catalogs | STAC |

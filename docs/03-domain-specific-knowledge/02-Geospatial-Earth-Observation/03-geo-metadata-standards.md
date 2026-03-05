---
title: Geospatial Metadata & Standards
slug: /domain-knowledge/geospatial/metadata-standards
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Geospatial Metadata & Standards**

Proper metadata is critical in geospatial science — without knowing a dataset's coordinate reference system, temporal extent, or spatial resolution, the data is effectively unusable. Several well-established standards govern how geospatial metadata is created and shared.

---

## **Core Metadata Standards**

### **ISO 19115 / ISO 19139**
The international standard for geographic information metadata. ISO 19115 defines the schema; ISO 19139 provides the XML encoding.
- **Use:** Official metadata records, data portals, INSPIRE compliance
- **Key fields:** Title, abstract, spatial extent, temporal extent, CRS, lineage, contact information
- **Spec:** [ISO 19115](https://www.iso.org/standard/53798.html)

### **INSPIRE Metadata**
The European Union's spatial data infrastructure directive requires metadata conforming to ISO 19115 with additional INSPIRE-specific requirements.
- **Use:** European spatial data sharing and discovery
- **Spec:** [INSPIRE Technical Guidelines](https://inspire.ec.europa.eu/metadata/6541)

### **Dublin Core (with spatial extensions)**
A general-purpose metadata standard that can be extended with geospatial elements (coverage, spatial).
- **Use:** Simple metadata for data catalogs and repositories
- **Spec:** [Dublin Core](https://www.dublincore.org/)

### **STAC (SpatioTemporal Asset Catalog)**
A modern, JSON-based metadata specification designed for Earth observation data. Describes spatiotemporal assets with standardized fields and extensions.
- **Use:** Satellite data catalogs, cloud-native data discovery
- **Key fields:** Bounding box, datetime, assets (links to data files), properties
- **Spec:** [STAC Spec](https://stacspec.org/)

---

## **Coordinate Reference Systems (CRS)**

Understanding CRS is fundamental to all geospatial work:

- **WGS 84 (EPSG:4326)** — The most common geographic CRS (latitude/longitude), used by GPS
- **UTM Zones (EPSG:326xx)** — Projected CRS for metric measurements, divided into 60 zones globally
- **ETRS89 (EPSG:4258)** — European Terrestrial Reference System, standard for European data
- **Web Mercator (EPSG:3857)** — Used by web maps (Google Maps, OpenStreetMap), distorts area

**Best practice:** Always document the CRS (as an EPSG code) in your metadata. Use projected CRS (e.g., UTM) for area/distance calculations, and geographic CRS (e.g., WGS 84) for data exchange.

- **Lookup tool:** [EPSG.io](https://epsg.io/) — Search and convert CRS codes

---

## **Ontologies and Vocabularies**

### **CF Conventions (Climate and Forecast)**
Standard names and units for climate and forecast data, used with NetCDF files.
- **Use:** Climate model output, reanalysis datasets
- **Spec:** [CF Conventions](https://cfconventions.org/)

### **GEMET (General Multilingual Environmental Thesaurus)**
A controlled vocabulary for environmental topics, used in INSPIRE metadata.
- **Spec:** [GEMET](https://www.eionet.europa.eu/gemet/)

### **GCMD Keywords (NASA)**
Hierarchical keywords for Earth science data discovery.
- **Spec:** [GCMD Keywords](https://earthdata.nasa.gov/earth-observation-data/find-data/idn/gcmd-keywords)

---

## **Best Practices**

1. **Always include CRS information** — Either embedded in the file (GeoTIFF) or documented in metadata
2. **Document spatial and temporal extent** — Bounding box, date range, temporal resolution
3. **Specify data lineage** — Processing chain from raw data to analysis-ready products
4. **Use standard vocabularies** — CF Conventions for climate data, GCMD keywords for discovery
5. **Provide data quality information** — Accuracy, completeness, positional uncertainty

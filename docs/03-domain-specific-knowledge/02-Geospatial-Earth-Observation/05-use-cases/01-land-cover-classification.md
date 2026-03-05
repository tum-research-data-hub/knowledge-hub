---
title: Land Cover Classification from Satellite Imagery
slug: /domain-knowledge/geospatial/use-cases/land-cover-classification
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Use Case 1: Land Cover Classification from Satellite Imagery**

This analysis takes multispectral satellite imagery (e.g., Sentinel-2) and assigns each pixel or region to a land cover class such as forest, water, urban, or agricultural land. It is one of the most fundamental tasks in remote sensing and underpins applications from environmental monitoring to urban planning.

## **Workflow Overview**

```
Data acquisition → Preprocessing → Classification → Validation
Satellite scenes (JPEG2000/GeoTIFF) → Analysis-ready data → Classified map (GeoTIFF) → Accuracy report
```

## **Key Concepts**

### **Preprocessing**

Raw satellite data must be corrected before analysis:
- **Atmospheric correction** — Converts top-of-atmosphere (TOA) reflectance to surface reflectance (e.g., using Sen2Cor for Sentinel-2)
- **Cloud masking** — Identifies and removes cloud-contaminated pixels
- **Mosaicking and compositing** — Combines multiple scenes to create cloud-free coverage over larger areas or time periods

### **Classification Approaches**

- **Pixel-based classification** — Each pixel is classified independently based on its spectral values
- **Object-based classification (OBIA)** — Pixels are first grouped into meaningful segments, then classified
- **Supervised classification** — Uses labeled training data (e.g., Random Forest, SVM, deep learning)
- **Unsupervised classification** — Clusters pixels based on spectral similarity without training data (e.g., k-means, ISODATA)

### **Spectral Indices**

Derived indices enhance specific features:
- **NDVI** (Normalized Difference Vegetation Index) — Highlights vegetation
- **NDWI** (Normalized Difference Water Index) — Highlights water bodies
- **NDBI** (Normalized Difference Built-up Index) — Highlights urban areas

## **Popular Tools**

- **Google Earth Engine** — Cloud-based processing of global satellite archives
- **QGIS + Semi-Automatic Classification Plugin** — Desktop-based classification
- **scikit-learn** (Python) — Machine learning classifiers
- **Rasterio + GeoPandas** (Python) — Raster/vector processing
- **terra + sf** (R) — Raster and vector analysis in R
- **SNAP** (ESA) — Sentinel data preprocessing

## **Code Example: Sentinel-2 Classification with Python**

<details>
<summary>**Click to expand Python workflow**</summary>

```python
#!/usr/bin/env python3
"""Land cover classification from Sentinel-2 using Random Forest."""

import rasterio
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import geopandas as gpd
from rasterio.mask import mask

# Step 1: Load Sentinel-2 bands (B2=Blue, B3=Green, B4=Red, B8=NIR)
bands = {}
for band_name, band_file in [("B2", "B02.tif"), ("B3", "B03.tif"),
                               ("B4", "B04.tif"), ("B8", "B08.tif")]:
    with rasterio.open(band_file) as src:
        bands[band_name] = src.read(1).astype(float)
        profile = src.profile

# Step 2: Calculate spectral indices
ndvi = (bands["B8"] - bands["B4"]) / (bands["B8"] + bands["B4"] + 1e-10)
ndwi = (bands["B3"] - bands["B8"]) / (bands["B3"] + bands["B8"] + 1e-10)

# Step 3: Stack features into a multi-band array
feature_stack = np.stack([bands["B2"], bands["B3"], bands["B4"],
                          bands["B8"], ndvi, ndwi], axis=0)

# Step 4: Load training polygons (GeoPackage with 'class' column)
training_data = gpd.read_file("training_polygons.gpkg")

# Step 5: Extract pixel values for training areas
# (simplified — in practice, use rasterio.mask or rasterstats)
X_train = []  # feature vectors
y_train = []  # class labels

# Step 6: Train Random Forest classifier
rf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)
rf.fit(X_train, y_train)

# Step 7: Classify the full image
n_bands, rows, cols = feature_stack.shape
flat_features = feature_stack.reshape(n_bands, -1).T
prediction = rf.predict(flat_features).reshape(rows, cols)

# Step 8: Save classified map as GeoTIFF
profile.update(dtype=rasterio.uint8, count=1, nodata=0)
with rasterio.open("land_cover_map.tif", "w", **profile) as dst:
    dst.write(prediction.astype(np.uint8), 1)

print("Classification complete! Output: land_cover_map.tif")
```

</details>

## **Expected Outputs**

- **Classified map** — GeoTIFF with integer values representing land cover classes
- **Accuracy assessment** — Confusion matrix, overall accuracy, kappa coefficient, per-class F1 scores
- **Feature importance** — Ranking of which bands/indices contributed most to classification

## **Computational Requirements**

| Task | CPU Cores | RAM | Storage | Time |
|------|-----------|-----|---------|------|
| Single Sentinel-2 tile (100x100 km) | 4-8 | 8-16 GB | ~5 GB | 10-30 min |
| Regional mosaic (country-scale) | 16-32 | 32-64 GB | 50-200 GB | Hours |
| Global analysis (via GEE) | Cloud | Cloud | Cloud | Hours-days |

## **Common Issues & Troubleshooting**

:::warning Common Problems

**Low classification accuracy**
- Ensure training samples are representative and balanced across classes
- Add more spectral indices or texture features
- Try object-based classification instead of pixel-based

**Cloud contamination**
- Use cloud-masked composites (median or best-pixel) rather than single scenes
- Apply the Sentinel-2 Scene Classification Layer (SCL) for masking

**CRS mismatch**
- Ensure all input layers use the same CRS before analysis
- Sentinel-2 tiles use UTM zones — adjacent tiles may be in different UTM zones

:::

## **Key Considerations**

- **Training data quality** is the single most important factor for classification accuracy
- **Temporal compositing** (e.g., seasonal medians) often outperforms single-date classification
- **Always validate** with independent test data, not just training data
- Consider using **existing land cover products** (e.g., ESA WorldCover, Corine) as baseline references

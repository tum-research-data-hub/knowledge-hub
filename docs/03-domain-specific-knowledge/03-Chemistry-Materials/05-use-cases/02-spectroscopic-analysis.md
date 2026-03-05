---
title: Spectroscopic Data Analysis
slug: /domain-knowledge/chemistry-materials/use-cases/spectroscopic-analysis
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Use Case 2: Spectroscopic Data Analysis**

Spectroscopic techniques measure how matter interacts with electromagnetic radiation, producing spectra that reveal molecular structure, composition, and concentration. This use case covers the management and analysis of spectroscopic data from techniques like NMR, IR/Raman, UV-Vis, and mass spectrometry.

## **Workflow Overview**

```
Measurement → Data conversion → Processing → Interpretation
Raw vendor data → Open format (JCAMP-DX/mzML) → Processed spectra → Identified compounds / quantification
```

## **Key Concepts**

### **Common Spectroscopic Techniques**

| Technique | What it Measures | Key Information | Typical File Size |
|-----------|-----------------|-----------------|-------------------|
| **NMR** | Nuclear spin environments | Molecular structure, connectivity | 1-100 MB per experiment |
| **IR / Raman** | Molecular vibrations | Functional groups, bond types | 10 KB - 1 MB |
| **UV-Vis** | Electronic transitions | Concentration, conjugation | 1-100 KB |
| **Mass Spectrometry (MS)** | Mass-to-charge ratios | Molecular weight, fragmentation | 100 MB - 10 GB |
| **X-ray Diffraction (XRD)** | Crystal lattice spacings | Crystal structure, phase identification | 1-10 MB |

### **Data Processing Steps**

- **Baseline correction** — Remove background signal drift
- **Normalization** — Scale spectra for comparison
- **Peak picking** — Identify signal positions and intensities
- **Phase correction** (NMR) — Correct phase distortions in frequency-domain spectra
- **Calibration** — Map instrument response to known standards
- **Deconvolution** — Separate overlapping peaks

### **Data Format Challenges**

Each instrument manufacturer uses proprietary binary formats. Converting to open formats is critical for:
- Long-term data preservation
- Cross-instrument comparison
- Sharing and reproducibility

## **Popular Tools**

### **General**
- **Python** (SciPy, NumPy, matplotlib) — Custom processing and visualization
- **R** (hyperSpec, ChemoSpec) — Chemometric analysis

### **NMR**
- **TopSpin** (Bruker) — Acquisition and processing
- **MestReNova** — Processing and analysis
- **nmrglue** (Python) — Read/process NMR data programmatically

### **Mass Spectrometry**
- **ProteoWizard / msconvert** — Convert vendor formats to mzML
- **OpenMS** — LC-MS data processing
- **MZmine** — Untargeted metabolomics

### **IR/Raman**
- **OPUS** (Bruker) — IR data processing
- **SpectraGryph** — Free viewer for multiple spectral formats

## **Code Example: Processing IR Spectra with Python**

<details>
<summary>**Click to expand Python workflow**</summary>

```python
#!/usr/bin/env python3
"""Process and compare IR spectra from JCAMP-DX files."""

import numpy as np
from scipy.signal import savgol_filter, find_peaks
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from pathlib import Path

def read_jcamp(filepath: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a JCAMP-DX file and return wavenumber and absorbance arrays."""
    wavenumbers = []
    absorbance = []
    in_data = False

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("##XYDATA") or line.startswith("##XYPOINTS"):
                in_data = True
                continue
            if line.startswith("##END"):
                break
            if in_data and not line.startswith("##"):
                parts = line.split()
                if len(parts) >= 2:
                    wavenumbers.append(float(parts[0]))
                    absorbance.append(float(parts[1]))

    return np.array(wavenumbers), np.array(absorbance)


# Step 1: Load spectra
wn_sample, abs_sample = read_jcamp(Path("sample.jdx"))
wn_ref, abs_ref = read_jcamp(Path("reference.jdx"))

# Step 2: Baseline correction (simple rubber band method)
def rubber_band_baseline(wavenumbers, spectrum, n_points=50):
    """Simple convex hull baseline correction."""
    from scipy.spatial import ConvexHull
    points = np.column_stack([wavenumbers, spectrum])
    hull = ConvexHull(points)
    baseline_idx = sorted(set(hull.vertices))
    baseline_fn = interp1d(wavenumbers[baseline_idx], spectrum[baseline_idx],
                           kind="linear", fill_value="extrapolate")
    return spectrum - baseline_fn(wavenumbers)

abs_sample_corr = rubber_band_baseline(wn_sample, abs_sample)

# Step 3: Smooth with Savitzky-Golay filter
abs_smooth = savgol_filter(abs_sample_corr, window_length=15, polyorder=3)

# Step 4: Normalize (min-max)
abs_norm = (abs_smooth - abs_smooth.min()) / (abs_smooth.max() - abs_smooth.min())

# Step 5: Find peaks
peaks, properties = find_peaks(abs_norm, height=0.1, prominence=0.05)
peak_wavenumbers = wn_sample[peaks]

print(f"Found {len(peaks)} peaks at: {peak_wavenumbers}")

# Step 6: Plot
fig, axes = plt.subplots(2, 1, figsize=(12, 8))

axes[0].plot(wn_sample, abs_sample, label="Raw", alpha=0.5)
axes[0].plot(wn_sample, abs_sample_corr, label="Baseline-corrected")
axes[0].set_xlabel("Wavenumber (cm⁻¹)")
axes[0].set_ylabel("Absorbance")
axes[0].invert_xaxis()
axes[0].legend()
axes[0].set_title("Raw vs. Baseline-Corrected Spectrum")

axes[1].plot(wn_sample, abs_norm, label="Processed sample")
axes[1].plot(wn_sample[peaks], abs_norm[peaks], "rx", markersize=8, label="Detected peaks")
axes[1].set_xlabel("Wavenumber (cm⁻¹)")
axes[1].set_ylabel("Normalized Absorbance")
axes[1].invert_xaxis()
axes[1].legend()
axes[1].set_title("Processed Spectrum with Peak Detection")

plt.tight_layout()
plt.savefig("ir_analysis.pdf", dpi=300)
plt.savefig("ir_analysis.png", dpi=300)
print("Analysis complete!")
```

</details>

## **Expected Outputs**

- **Processed spectra** — Baseline-corrected, normalized spectra in open formats (JCAMP-DX, CSV)
- **Peak lists** — Tables of peak positions, intensities, and assignments
- **Comparison plots** — Overlay of sample vs. reference spectra (PDF + PNG)
- **Quantification results** — Concentrations from calibration curves
- **Structure assignments** — Molecular structure determined from spectral evidence

## **Computational Requirements**

| Task | CPU | RAM | Storage | Time |
|------|-----|-----|---------|------|
| Single spectrum processing | 1 | `<1 GB` | Minimal | Seconds |
| Batch processing (100s of spectra) | 2-4 | 2-4 GB | 1-10 GB | Minutes |
| LC-MS untargeted metabolomics | 4-8 | 8-32 GB | 10-100 GB | Hours |
| NMR structure elucidation (2D) | 2-4 | 4-8 GB | 1-5 GB | Minutes-hours |

## **Common Issues & Troubleshooting**

:::warning Common Problems

**Vendor format lock-in**
- Convert to open formats (JCAMP-DX for spectroscopy, mzML for MS) as early as possible
- Use ProteoWizard's `msconvert` for mass spectrometry data
- Document the conversion process and any parameters used

**Poor baseline**
- Automated baseline correction may fail for noisy spectra — try different algorithms
- Manual baseline correction may be needed for complex backgrounds
- Always visually inspect results

**Peak overlap**
- Use deconvolution (curve fitting) to separate overlapping signals
- Consider 2D techniques (e.g., 2D NMR) for complex mixtures

:::

## **Key Considerations**

- **Convert vendor formats early** — Raw proprietary files may become unreadable as software evolves
- **Store raw data alongside processed results** — Never discard the original measurement
- **Document instrument parameters** — Acquisition settings are essential for reproducibility
- **Use reference databases** — Compare spectra against SDBS, NIST, or MassBank for identification
- **Consider electronic lab notebooks** — Link spectra to sample preparation records using eLabFTW or Chemotion

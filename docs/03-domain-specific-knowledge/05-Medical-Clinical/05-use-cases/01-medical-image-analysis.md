---
title: Medical Image Analysis
slug: /domain-knowledge/medical-clinical/use-cases/medical-image-analysis
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Use Case 1: Medical Image Analysis**

Medical image analysis involves processing and interpreting imaging data (CT, MRI, X-ray, pathology slides) to support diagnosis, treatment planning, and research. This use case highlights the unique data management challenges — from DICOM handling and anonymization to deep learning-based segmentation.

## **Workflow Overview**

```
Acquisition → Anonymization → Preprocessing → Analysis → Reporting
DICOM (hospital PACS) → De-identified DICOM/NIfTI → Standardized images → Segmentation / classification → Results
```

## **Key Concepts**

### **DICOM and Data Flow**

Medical images originate from hospital Picture Archiving and Communication Systems (PACS) in DICOM format. Before research use:
1. **Export** from PACS (usually via institutional radiology department)
2. **Anonymize** — Remove all patient-identifiable information (name, ID, dates, institution)
3. **Convert** — Typically to NIfTI for neuroimaging or keep as anonymous DICOM for radiology
4. **Organize** — Structure files in a consistent directory hierarchy (e.g., BIDS for neuroimaging)

### **Anonymization**

DICOM files embed extensive patient metadata. De-identification must handle:
- **Direct identifiers** — Patient name, ID, birth date, address
- **Indirect identifiers** — Referring physician, institution name, study date
- **Burned-in annotations** — Text overlaid on the pixel data itself (requires pixel scrubbing)
- **DICOM UIDs** — Unique identifiers that can be used for re-identification

### **Common Analysis Tasks**

- **Segmentation** — Delineating anatomical structures or lesions (manual, semi-automatic, or AI-based)
- **Classification** — Predicting disease presence or severity from images
- **Registration** — Aligning images from different time points or modalities
- **Quantification** — Measuring volumes, areas, signal intensities
- **Radiomics** — Extracting quantitative features from images for statistical modeling

## **Popular Tools**

### **Viewing and Annotation**
- **3D Slicer** — Open-source platform for visualization, segmentation, and analysis
- **OHIF Viewer** — Web-based DICOM viewer
- **ITK-SNAP** — Interactive segmentation tool
- **QuPath** — Digital pathology analysis

### **Processing Libraries**
- **SimpleITK / ITK** — Image processing (filtering, registration, segmentation)
- **nibabel** (Python) — Read/write NIfTI and other neuroimaging formats
- **pydicom** (Python) — Read/write DICOM files
- **MONAI** (Python) — Deep learning framework for medical imaging built on PyTorch

### **Anonymization**
- **CTP (Clinical Trial Processor)** — RSNA's DICOM anonymization pipeline
- **deid** (Python) — Flexible DICOM de-identification library
- **DicomCleaner** — PixelMed's standalone anonymization tool

## **Code Example: DICOM Processing with Python**

<details>
<summary>**Click to expand Python DICOM workflow**</summary>

```python
#!/usr/bin/env python3
"""Process DICOM files: anonymize, convert to NIfTI, basic analysis."""

from pathlib import Path
import pydicom
from pydicom.uid import generate_uid
import SimpleITK as sitk
import numpy as np

# ============================================================
# Step 1: Anonymize DICOM files
# ============================================================

def anonymize_dicom(input_dir: Path, output_dir: Path) -> None:
    """Remove patient-identifiable information from DICOM files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    tags_to_remove = [
        "PatientName", "PatientID", "PatientBirthDate",
        "PatientAddress", "PatientTelephoneNumbers",
        "ReferringPhysicianName", "InstitutionName",
        "InstitutionAddress", "StudyDate", "StudyTime",
    ]
    
    for dcm_path in sorted(input_dir.glob("**/*.dcm")):
        ds = pydicom.dcmread(dcm_path)
        
        for tag_name in tags_to_remove:
            if hasattr(ds, tag_name):
                setattr(ds, tag_name, "ANONYMIZED")
        
        # Replace UIDs to prevent re-linking
        ds.StudyInstanceUID = generate_uid()
        ds.SeriesInstanceUID = generate_uid()
        ds.SOPInstanceUID = generate_uid()
        
        output_path = output_dir / dcm_path.name
        ds.save_as(output_path)
    
    print(f"Anonymized {len(list(input_dir.glob('**/*.dcm')))} files")

anonymize_dicom(Path("raw_dicom/"), Path("anon_dicom/"))

# ============================================================
# Step 2: Convert DICOM series to NIfTI
# ============================================================

reader = sitk.ImageSeriesReader()
dicom_names = reader.GetGDCMSeriesFileNames("anon_dicom/")
reader.SetFileNames(dicom_names)
image = reader.Execute()

sitk.WriteImage(image, "patient_scan.nii.gz")
print(f"Image size: {image.GetSize()}")
print(f"Spacing: {image.GetSpacing()} mm")

# ============================================================
# Step 3: Basic analysis — threshold segmentation
# ============================================================

# Simple threshold segmentation (e.g., bone from CT)
image_array = sitk.GetArrayFromImage(image)

# Bone threshold for CT (Hounsfield units > 300)
bone_mask = sitk.BinaryThreshold(image, lowerThreshold=300,
                                  upperThreshold=3000,
                                  insideValue=1, outsideValue=0)

# Calculate volume
spacing = image.GetSpacing()
voxel_volume_mm3 = spacing[0] * spacing[1] * spacing[2]
bone_voxels = sitk.GetArrayFromImage(bone_mask).sum()
bone_volume_cm3 = bone_voxels * voxel_volume_mm3 / 1000

print(f"Bone volume: {bone_volume_cm3:.1f} cm³")

# Save segmentation mask
sitk.WriteImage(bone_mask, "bone_segmentation.nii.gz")
print("Analysis complete!")
```

</details>

## **Code Example: Deep Learning Segmentation with MONAI**

<details>
<summary>**Click to expand MONAI segmentation pipeline**</summary>

```python
#!/usr/bin/env python3
"""Organ segmentation using MONAI and a pretrained U-Net."""

import torch
from monai.transforms import (
    Compose, LoadImaged, EnsureChannelFirstd,
    Spacingd, ScaleIntensityRanged, CropForegroundd,
)
from monai.networks.nets import UNet
from monai.inferers import sliding_window_inference
import nibabel as nib
import numpy as np

# Step 1: Define preprocessing transforms
transforms = Compose([
    LoadImaged(keys=["image"]),
    EnsureChannelFirstd(keys=["image"]),
    Spacingd(keys=["image"], pixdim=(1.5, 1.5, 2.0), mode="bilinear"),
    ScaleIntensityRanged(keys=["image"], a_min=-175, a_max=250,
                          b_min=0.0, b_max=1.0, clip=True),
])

# Step 2: Load and preprocess image
data = transforms({"image": "patient_scan.nii.gz"})
image_tensor = data["image"].unsqueeze(0)  # Add batch dimension

# Step 3: Load model
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = UNet(
    spatial_dims=3, in_channels=1, out_channels=14,
    channels=(16, 32, 64, 128, 256), strides=(2, 2, 2, 2),
    num_res_units=2,
).to(device)

# Load pretrained weights (example)
# model.load_state_dict(torch.load("pretrained_segmentation.pth"))
model.eval()

# Step 4: Run inference with sliding window
with torch.no_grad():
    output = sliding_window_inference(
        image_tensor.to(device), roi_size=(96, 96, 96),
        sw_batch_size=4, predictor=model,
    )
    prediction = torch.argmax(output, dim=1).cpu().numpy()[0]

# Step 5: Save segmentation result
seg_nifti = nib.Nifti1Image(prediction.astype(np.uint8),
                              nib.load("patient_scan.nii.gz").affine)
nib.save(seg_nifti, "organ_segmentation.nii.gz")

print(f"Segmentation saved with {len(np.unique(prediction))} classes")
print("Open in 3D Slicer for visualization")
```

</details>

## **Expected Outputs**

- **Anonymized images** — De-identified DICOM or NIfTI files safe for research use
- **Segmentation masks** — Label maps identifying anatomical structures or lesions
- **Quantitative measurements** — Volumes, diameters, signal intensities (CSV)
- **Classification results** — Predicted diagnoses with confidence scores
- **Radiomics features** — Extracted texture, shape, and intensity features (CSV)

## **Computational Requirements**

| Task | CPU Cores | GPU | RAM | Storage | Time |
|------|-----------|-----|-----|---------|------|
| DICOM anonymization (1000 files) | 2-4 | - | 4-8 GB | 1-5 GB | Minutes |
| Classical segmentation (single scan) | 4 | - | 8-16 GB | `<1 GB` | Minutes |
| Deep learning inference (single scan) | 4 | 1 (8+ GB VRAM) | 16-32 GB | `<1 GB` | 1-10 min |
| DL model training (100+ scans) | 8-16 | 1-4 (16+ GB VRAM) | 32-128 GB | 50-500 GB | Hours-days |

## **Common Issues & Troubleshooting**

:::warning Common Problems

**DICOM inconsistencies**
- Different scanners produce different DICOM tags — normalize before analysis
- Missing slices or inconsistent slice spacing — validate series completeness
- Use `dcmcheck` or `dciodvfy` tools to validate DICOM conformance

**Anonymization failures**
- Burned-in annotations on images require pixel-level de-identification
- Secondary capture images may contain patient info in unexpected locations
- Always validate anonymization with a DICOM tag viewer before sharing

**Memory errors during deep learning**
- Reduce patch/window size in sliding window inference
- Use mixed precision (FP16) training to halve memory usage
- For very large 3D images, process in overlapping patches

:::

## **Key Considerations**

- **Ethics approval first** — Ensure institutional ethics board approval before accessing patient data
- **Anonymize as early as possible** — De-identify data immediately upon export from clinical systems
- **Document your pipeline** — Every processing step must be reproducible
- **Separate raw and processed data** — Never modify original DICOM files
- **Consider federated approaches** — If data cannot leave the hospital, bring the analysis to the data

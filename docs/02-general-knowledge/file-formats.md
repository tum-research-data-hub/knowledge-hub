---
title: File formats
slug: /general-knowledge/file-formats
---

A file format specifies the way information is encoded in a file on disk. Depending on the type of data that are to be stored, different file formats are appropriate.

Choosing the most optimal file format for research data can lead to drastic improvements in terms of used disk space, used memory while reading and processing speed.

Some data are so domain-specific that new file formats have been invented but there are lots of file formats that serve a general purpose and lend themselves to a wide variety of usecases in a scientific context.

## Text-based file formats
Text-based file formats are, while not necessarily designed to be written by a human, human-readable as the information they contain is encoded as text that may be inspected or modified via a simple text editor.
They are optimal to store information that is meant to be read or changed by a human such as configuration files, small amounts of data that need to be stored in a most portable way and of course any data that are inherently related to text or speech.

| Name | File suffix | Notes | How to read |
| --- | --- | --- | --- |
| JSON | `.json` | nested data | text editor / json library in any programming language |
| yaml | `.yaml` / `.yml` | nested data | text editor / yaml library in any programming language |
| XML | `.xml` | nested data | text editor / xml library in any programming language |
| CSV | `.csv` | tabular data | text editor / python (pandas, polars) / csv library in any programming language |

## Binary file formats
While technically all files on a computer are binary, "binary files" refers to all file formats that are not text.
Even though it is possible to find a textual representation for all kinds of data, doing so introduces overhead both in disk space and needed processing power.
Hence, storing data in a binary format is mostly a matter of efficiency, although it can have an enormous impact.

Here we only list general purpose file formats for tabular data:

| Name | File suffix | Notes | How to read |
| --- | --- | --- | --- |
| [Parquet](https://parquet.apache.org/) | `.parquet` | columnar storage; supports boolean and numerical data, byte arrays, geospatial data, nested data and logical types | python (pandas, polars), Apache Arrow (C/C++, Julia, MATLAB, Python, R and many more) |
| [HDF5](https://www.hdfgroup.org/) | `.hdf5` / `.h5` | hierarchical format for large numerical arrays and datasets; widely used in scientific computing and deep learning; supports metadata and compression | python (h5py, pandas), MATLAB, R (rhdf5), C/C++, Julia |
| [Arrow IPC / Feather](https://arrow.apache.org/) | `.arrow` / `.feather` | columnar in-memory format; optimized for fast read/write with zero-copy access; less compression than Parquet but faster I/O | python (pandas, polars), Apache Arrow (C/C++, Julia, R and many more) |
| [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) | `.nc` / `.nc4` | array-oriented format popular in climate science, meteorology, and oceanography; supports metadata and unlimited dimensions | python (netCDF4, xarray), MATLAB, R (ncdf4), C/C++, Fortran, Java |
| [Zarr](https://zarr.dev/) | `.zarr` (directory) | chunked, compressed N-dimensional arrays; cloud-friendly (S3, GCS); can serve as an alternative to HDF5 for large-scale parallel I/O | python (zarr, xarray), Julia |
| [NumPy binary](https://numpy.org/doc/stable/reference/generated/numpy.save.html) | `.npy` / `.npz` | simple format for serializing NumPy arrays; `.npz` bundles multiple arrays in a zip archive | python (numpy) |

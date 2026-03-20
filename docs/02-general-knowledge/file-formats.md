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
| XML | `.xml` | nested data |  |
| CSV | `.csv` | tabular data | text editor / python (pandas, polars) / csv library in any programming language |

## Binary file formats
While technically all files on a computer are binary, "binary files" refers to all file formats that are not text.
Even though it is possible to find a textual representation for all kinds of data, doing so introduces overhead both in disk space and needed processing power.
Hence, storing data in a binary format is mostly a matter of efficiency, although it can have an enormous impact.

Here we only list general purpose file formats for tabular data:

| Name | File suffix | Notes | How to read |
| --- | --- | --- | --- |
| [parquet](https://parquet.apache.org/) | `.parquet` | supports boolean and numerical data, byte arrays, geospatial data, nested data and logical types  | python (pandas, polars), apache arrow (implementations available in C/C++, Julia, Matlab, Python, R and many more) |
| HDF5 |  |  |  |

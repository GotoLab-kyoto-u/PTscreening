# PTscreening

## Overview
This script aggregates data exported from UNIFI and extracts relevant information to calculate the conversion rate for each combination of peptide, prenyl donor, and prenyltransferase.
The pipeline automates peak parsing, reference-based assignment, and downstream conversion rate calculation from high-throughput screening datasets.

## Directory Structure
```text
main/
├─ src/
│  ├─ extract_data.py
│  └─ analyze_data.py
├─ example/
│  ├─ UNIFI_raw_data/
│  ├─ reference_data/
│  ├─ rename.csv
│  ├─ protein_pt_list.csv
│  ├─ pt_name_list.csv
│  ├─ expected_output/
│  └─ README.md
├─ LICENSE
├─ README.md
└─ requirements.txt
```


## Input Data Specification
PTscreening requires the following input files:

### 1. UNIFI raw data (.xls)
- Exported from Waters UNIFI software  
⚠ Note: The internal structure of `.xls` files depends on UNIFI export settings.
Consistent export format must be used across all samples.

### 2. reference_data/
Reference table used for peak annotation and assignment.  
This file defines expected retention times for unmodified peptides and is used to match LC-MS peaks from UNIFI raw data.
#### Required columns:
- Protein name  
  Identifier of the peptide
- Ref RT (min)  
  Reference retention time in minutes used for peak matching

### 3. CSV file listing peptide and prenyltransferase names used in the analysis
In this documentation, it is referred to as `protein_pt_list.csv` as an example.
#### Required columns:
- Protein name
  peptide names
- PT_Name
  prenyltransferase names

### 4. CSV file used for renaming sheet names (optional)
A user-defined CSV file used to map raw UNIFI sheet names to standardized identifiers for downstream analysis.  
This file is optional and can have any filename.  
In this documentation, it is referred to as `rename.csv` as an example.
#### Required columns:
- original name
- standardized name
#### Behavior:
- If "standardized name" is empty, the corresponding sample is ignored during processing.
- Such samples are not parsed and are excluded from all downstream analyses.
- This mechanism is intentionally used to allow selective exclusion of samples without modifying raw UNIFI data.

### 5. CSV file used for renaming prenyltransferase names (optional)
Mapping table between prenyltransferase IDs and prenyltransferase names used for annotation consistency.  
In this documentation, it is referred to as `pt_name_list.csv` as an example.
#### Required columns:
- PT_No
  prenyltransferase IDs
- PT_name
  prenyltransferase names

## Tested Environment
- Python 3.14.2
- pandas 2.3.3
- pyexcel 0.7.4

## Example Data
- `UNIFI_raw_data/`: Example raw data exported from UNIFI.
- `reference_data/`: Reference files used for peak assignment.
- `rename.csv`: CSV file used for renaming sheet names.
- `protein_pt_list.csv`: CSV file listing peptide and prenyltransferase names used in the analysis.
- `pt_name_list.csv`: CSV file used for renaming prenyltransferase names.

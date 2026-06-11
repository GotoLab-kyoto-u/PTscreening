# PTscreening Example Dataset

This directory provides a minimal example dataset for testing PTscreening.  
It demonstrates how to run the pipeline using example UNIFI-exported LC-MS data.  
⚠ This dataset is only for testing the software and is not intended for publication.

---

## Directory Structure

```text
example/
├─ UNIFI_raw_data/        # Example raw UNIFI .xls files
├─ reference_data/        # Reference table used for peak annotation and assignment
├─ rename.csv             # (optional) CSV file used for renaming sheet names
├─ protein_pt_list.csv    # CSV file listing peptide and prenyltransferase names used in the analysis
├─ pt_name_list.csv       # (optional) CSV file used for renaming prenyltransferase names
├─ expected_output/       # Expected output results
└─ README.md
```

## Step 1: Install dependencies
Run the following command from the root directory:

```bash
pip install -r requirements.txt
```

## Step 2: Run data extraction
Execute the following command:

```bash
python src/extract_data.py \
  -i example/UNIFI_raw_data/* \
  -r example/reference_data/* \
  -o example/extracted_data \
  -s example/rename.csv
```

This step generates intermediate processed data.  
All outputs are written to the directory specified by `-o / --output`.
- `all_unknownPT_data.csv`
  - Contains all detected modified and unmodified peptides observed in LC-MS data from reaction mixtures containing unknown prenyltransferases (identifiers starting with "PT").
- `all_knownPT_data.csv`
  - Contains all detected modified and unmodified peptides observed in LC-MS data from reaction mixtures containing known prenyltransferases (suffix "*F").
- `sample_sheet_list.csv`
  - Summary table of processed sheets

## Step 3: Run analysis
Execute the following command:

```bash
python src/analyze_data.py \
  -i example/extracted_data/all_*knownPT_data.csv \
  -o example/result \
  -l example/protein_pt_list.csv \
  -s example/sample_sheet_list.csv \
  -r example/pt_name_list.csv
```

The script generates processed conversion analysis results from extracted UNIFI LC-MS data.  
All outputs are written to the directory specified by `-o / --output`.
- Main processed dataset
  - `final_df.csv`
    - Contains the final aggregated results for each peptide × prenyltransferase × prenyl donor condition.
- Pivot tables (summary matrices)  
  Matrix-formatted CSV files are generated, with all metrics aggregated at the peptide × prenyltransferase × prenyl donor level.
  - `pct_modification.csv`
    - Average conversion rate (%modification) per condition.
    - %modification represents the proportion of LC-MS signal corresponding to prenylated product peaks relative to the sum of substrate and product peak intensities.
  - `pct_biproduct.csv`
    - Average fraction of minor product peaks (%biproduct) per condition.
    - %biproduct represents the fraction of all non-dominant product peaks, defined as all prenylated product signals other than the peak with the highest intensity.
  - `No_prenylation.csv`
    - Number of prenylation events per condition.
  - `No_replicate.csv`
    - Number of replicates per condition.

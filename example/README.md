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
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract data from UNIFI data.

This script aggregates data exported from UNIFI
and extracts relevant information.
Warnings are issued if inappropriate files or sheets are read,
or if duplicate data is detected.

Usage: python3 ./src/extract_data.py [-h] [-i [INPUT ...]]
       [-s SHEET] [-r [REF ...]] [-o OUTPUT] [--separate]

options:
    -h, --help            show this help message and exit
    -i, --input [INPUT ...]
                          Path to the UNIFI-exported .xls files
    -s, --sheet SHEET     Path to the CSV file used for renaming sheet names
    -r, --ref [REF ...]   Path to the reference .xlsx files
    -o, --output OUTPUT   Path to the output directory for the processed data
    --separate            Decide whether to output the results separately
                          for each peptide library

Output:
    Processed CSV files are written to the specified output directory.
    If the --separate flag is used,
    results are output separately for each peptide library.

"""

import os
import glob
import pandas as pd
import re
import argparse
import warnings
import pyexcel
import copy
from collections import Counter

warnings.filterwarnings(
    "ignore",
    message="Cannot parse header or footer so it will be ignored",
    category=UserWarning,
    module="openpyxl.worksheet.header_footer"
)

Cxprenylation_dict = {"D": "C5prenylation", "G": "C10prenylation",
                      "F": "C15prenylation"}


def read_xls(xls_path):
    """Read an XLS file and convert each sheet to a pandas DataFrame.

    Args:
        xls_path (str): Path to the .xls file.

    Returns:
        dict: A dictionary where keys are sheet names (str)
              and values are pandas DataFrames.
              Each DataFrame contains the sheet's data,
              with the first row used as column names.

    """

    book = pyexcel.get_book(file_name=xls_path)

    dfs = {}
    for sheet_name in book.sheet_names():
        if sheet_name and sheet_name.strip():
            data = book[sheet_name].get_array()
            dfs[sheet_name] = pd.DataFrame(data[1:], columns=data[0])

    return dfs


def edit_sheetname(sheet_name_original):
    """Organize the sheet name.

    Args:
        sheet_name_original (str): Original sheet name.

    Returns:
        list: List of edited sheet name.

    """

    if sheet_name_original in sheet_name_dict:
        sheet_name_dict_tmp = copy.deepcopy(sheet_name_dict)
        sheet_name_edit_l = sheet_name_dict_tmp[sheet_name_original]
    else:
        sheet_name_edit_l = [sheet_name_original]
    for i in range(len(sheet_name_edit_l)):
        sheet_name_edit = sheet_name_edit_l[i]
        pattern = r"^Sample_\d{6}_[^_]+_(PT\d+|[^_]{3}F)_[DGF](?:_.+)*$"

        if not re.match(pattern, sheet_name_edit):
            parts = sheet_name_edit.split("_")
            # Classify each element
            date = next((p for p in parts if re.fullmatch(r"\d{6}", p)), None)
            if date is None:
                matches = re.findall(r"\d{6}", os.path.basename(data_file))
                date = matches[-1] if matches else None

            pt = next((p for p in parts if re.fullmatch(r"PT\d+|[^_]{3}F", p)),
                      None)
            if pt is None:
                matches = re.findall(r"PT\d+|[^_]{3}F", os.path.basename(data_file))
                pt = matches[-1] if matches else None

            dgf = next((p for p in parts if p in ["D", "G", "F"]), None)
            if dgf is None:
                matches = re.findall(r"_[DGF]", os.path.basename(data_file))
                dgf = matches[-1].replace("_", "") if matches else None

            peptide_library = next((p for p in parts if p in reference_file),
                                   None)
            if peptide_library is None:
                peptide_library = next(
                    (p for p in os.path.basename(data_file).split("_")
                     if p in reference_file), None
                )

            # Check if all elements are present
            if not (date and pt and dgf and peptide_library):
                print(
                    f"!!! Invalid File Name ( {sheet_name_original} / "
                    f"{sheet_name_edit} ) !!!"
                )
                sheet_name_edit_l[i] = None
                continue
            rest = [p for p in parts[1:]
                    if p not in {date, pt, dgf, peptide_library}]

            # Reorder and reconstruct the sheet name
            new_name = f"Sample_{date}_{peptide_library}_{pt}_{dgf}"
            if rest:
                new_name += "_" + "_".join(rest)

            sheet_name_edit_l[i] = new_name

    return sheet_name_edit_l


def extract_data(data_file):
    """Extracts relevant data from UNIFI-exported .xls or .xlsx files.

    Args:
        data_file (str): Path to the input .xls or .xlsx file.

    Returns:
        tuple:
            bool: Success flag.
                  True if data extraction was successful, False otherwise.
            pandas.DataFrame or None:
                Extracted data combined from all valid sheets.
            list of str or None: List of processed sample sheet names.
            set of str or None:
                Set of peptide library names extracted from the file.

        Notes:
            - If the file format is unsupported or contains no valid data,
              the function returns (False, None, None, None).
    """

    ext = data_file.split(".")[-1]
    if ext == "xlsx":
        data_file_dict = pd.read_excel(data_file, sheet_name=None)
    elif ext == "xls":
        data_file_dict = read_xls(data_file)
    else:
        print("!!! Unsupported file format. Only .xls or .xlsx "
              "are supported !!!")
        return False, None, None, None

    # List to store rows that meet the specified conditions
    data_rows = []
    # List to store peptide library names
    peptide_library_l_tmp = []
    # List to store sample sheet names
    sample_sheet_l_tmp = []

    # Process each sheet in the data-containing files
    for sheet_name, data_sheet in data_file_dict.items():
        # Organize the sample sheet name
        sheet_name_edit_l = edit_sheetname(sheet_name)
        sheet_name_edit_l = [x for x in sheet_name_edit_l if x is not None]
        # print(sheet_name_edit_l)
        if sheet_name_edit_l == []:
            continue

        for sheet_name_edit in sheet_name_edit_l:
            # Add the sheet name to a list to check for duplicates
            sample_sheet_l_tmp.append(sheet_name_edit)

            # Add the peptide library name
            peptide_library_l_tmp.append(sheet_name_edit.split("_")[2])

            # Extract entries where the protein name matches peptides
            # listed in the reference file
            # and the label is ":1&" (i.e., excluding -NH2 and -H2O)
            merged_data_sheet = pd.merge(
                data_sheet,
                reference_file[sheet_name_edit.split("_")[2]][
                    ["Protein name", "Ref RT (min)"]],
                on="Protein name", how="left").copy()
            data_sheet_extract = merged_data_sheet[
                (merged_data_sheet["Protein name"].isin(reference_file[
                    sheet_name_edit.split("_")[2]]["Protein name"]))
                & (merged_data_sheet["Label"].isin([":1&"]))].copy()
            data_sheet_extract["Sheet_Name"] = sheet_name_edit
            # Calculate Expected mass from Observed mass and  Mass error (ppm)
            if "Expected mass (Da)" not in data_sheet_extract.columns:
                data_sheet_extract["Expected mass (Da)"] = (
                    data_sheet_extract["Observed mass (Da)"]
                    / (1 + data_sheet_extract["Mass error (ppm)"] * 1e-6)
                )
            if not data_sheet_extract.empty:
                data_rows.append(data_sheet_extract)

    # Combine modification_data and not_modification_data
    # into a single DataFrame
    if data_rows:
        data = pd.concat(data_rows)
    else:
        print("This file do not contain data.")
        return False, None, None, None

    data["PT_Name"] = (
        data["Sheet_Name"]
        .str.extract(r"(PT[0-9]*)", expand=False)
        .fillna(
            data["Sheet_Name"].str.extract(r"_([^_F]{3}F)", expand=False)
        )
    )
    data["Modifiers"] = data["Modifiers"].astype(str).fillna("")
    data["Cxprenylation"] = data["Modifiers"].str.extract(r"(C.*prenylation)")
    data["No_prenylation"] = data["Modifiers"].str.extract(r"\(([^)]+)\)")

    return True, data, sample_sheet_l_tmp, set(peptide_library_l_tmp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", default=["./data/UNIFI_raw_data/*"],
        help="Path to the UNIFI-exported .xls files", nargs="*"
    )
    parser.add_argument(
        "-s", "--sheet", default=None,
        help="Path to the CSV file used for renaming sheet names"
    )
    parser.add_argument(
        "-r", "--ref", default=["./data/reference_data/*.xlsx"],
        help="Path to the reference .xlsx files", nargs="*"
    )
    parser.add_argument(
        "-o", "--output", default="./data/extracted_data/",
        help="Path to the output directory for the processed data"
    )
    parser.add_argument(
        "--separate", action="store_true",
        help="Decide whether to output the results separately "
        "for each peptide library")
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    data_files = []
    for x in args.input:
        data_files += glob.glob(x)
    # print(data_files)

    sheet_name_dict = {}
    if args.sheet is not None:
        with open(args.sheet, "r") as f:
            for line in f:
                data = line.rstrip().split(",")
                sheet_name_dict.setdefault(data[0], []).append(data[1])
        # print(sheet_name_dict)

    ref_files = []
    for x in args.ref:
        ref_files += glob.glob(x)
    # print(ref_files)
    reference_file = {
        x.split("/")[-1].split("_reference")[0]:
        pd.read_excel(x) for x in ref_files
    }

    sample_sheet_l = []        # List to check for duplicates
    peptide_library_l = set()  # Set to store peptide library name
    data_l = []                # List to store data

    for data_file in data_files:
        print(data_file)
        (
            ok,
            new_data,
            new_sample_sheet_l,
            new_peptide_library_l
        ) = extract_data(data_file)
        if ok:
            data_l += [new_data]
            sample_sheet_l += new_sample_sheet_l
            peptide_library_l = peptide_library_l.union(new_peptide_library_l)

    # Combine all data
    data_all = pd.concat(data_l, ignore_index=True)

    # Extract unkown PTs
    data_unknown = data_all[
        data_all["PT_Name"].str.contains("PT", na=False)
        ].reset_index(drop=True)
    data_unknown.to_csv(os.path.join(args.output, "all_unknownPT_data.csv"),
                        index=False)
    # Extract kown PTs
    data_known = data_all[
        data_all["PT_Name"].str.fullmatch(r"([^_F]{3}F)", na=False)
        ].reset_index(drop=True)
    data_known.to_csv(os.path.join(args.output, "all_knownPT_data.csv"),
                      index=False)

    if args.separate:
        for peptide_library_name in peptide_library_l:
            data_unknown_extract = data_unknown[
                data_unknown["Sheet_Name"].str.split("_").str[2]
                == peptide_library_name].reset_index(drop=True)
            data_unknown_extract.to_csv(
                os.path.join(
                    args.output, f"{peptide_library_name}_unknownPT_data.csv"
                ), index=False
            )
            data_known_extract = data_known[
                data_known["Sheet_Name"].str.split("_").str[2]
                == peptide_library_name].reset_index(drop=True)
            data_known_extract.to_csv(
                os.path.join(
                    args.output, f"{peptide_library_name}_knownPT_data.csv"
                ), index=False
            )
    # Check for duplicate sheet names
    if len(sample_sheet_l) != len(set(sample_sheet_l)):
        counts = Counter(sample_sheet_l)
        # Extract elements that are duplicated
        duplicates = [item for item, count in counts.items() if count > 1]
        print("The duplicated sheets are listed below for your review.")
        print(duplicates)

    with open(args.output+"/sample_sheet_list.csv", "w") as f:
        f.write("Sheet_Name,Protein_name_list,PT_Name,Cxprenylation\n")
        for x in set(sample_sheet_l):
            f.write(
                f"{x},{x.split('_')[2]},{x.split('_')[3]},"
                f"{Cxprenylation_dict[x.split('_')[4]]}\n"
            )
    print("\nCompleted!")

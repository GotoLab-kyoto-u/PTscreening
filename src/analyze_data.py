#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analyze data from UNIFI data.

This script processes extracted UNIFI data to obtain relevant
peptide and prenyltransferase information,
and calculates the conversion rate based on specified criteria.
It aggregates data from CSV files of the extracted data,
performs filtering, calculates conversion rates,
and generates output CSV files.

Usage: python3 ./src/analyze_data.py [-h] [-i [INPUT ...]] [-o OUTPUT]
       [-l LIST] [-s SAMPLESHEET] [-r RENAME] [-e]

options:
    -h, --help            show this help message and exit
    -i, --input [INPUT ...]
                          Path to the CSV files of the extracted data
    -o, --output OUTPUT   Path to the output directory for the processed data
    -l, --list LIST       CSV file path for prenyltransferase and peptide names
    -s, --samplesheet SAMPLESHEET
                          Path to sample_sheet_list.csv
    -r, --rename RENAME   CSV file path containing the correspondence table of
                          prenyltransferase names
    -e, --edit            Decide whether to make a simple csv file

Output:
    Processed CSV files are written to the specified output directory.

"""


import os
import glob
import argparse
import pandas as pd
from pandas.api.types import CategoricalDtype


# Change the following values to change the threshold
#  RT threshold for each Cx prenylation type
thresholds = {"C5prenylation": 1.5, "C10prenylation": 2.6,
              "C15prenylation": 3.5}
#  Threshold for mass error (ppm)
Mass_error_TH = 7
#  Minimum response value to consider a peak as non-noise
Response_TH = 10 ** 6
#  Conversion rate threshold below which the modification is likely to be noise
modification_TH = 0.05
#  RT tolerance (± RT_TH1 (min)) for selecting peaks
#  that match the observed RT (min) in the reference file
RT_TH1 = 0.2
#  RT tolerance (± RT_TH2 (min)) for determining
#  whether two peaks are considered identical
RT_TH2 = 0.1
#  RT tolerance (± RT_TH2 (min)) for determining
#  whether two peaks are considered identical
RT_TH2 = 0.1
#  The maximum number of prenyl groups that can be introduced
No_prenylation_TH = 2


# Change the following values for editing table
#  Dictionary for converting "Cxprenylation" to "prenyl donor"
Cxprenylation_dict = {"C5prenylation": "DMAPP", "C10prenylation": "GPP",
                      "C15prenylation": "FPP"}


def modify_data(combined_data, output_path):
    # Split data into modified and non-modified based on
    #  whether Cxprenylation has a value
    if "Cxprenylation" not in combined_data.columns:
        combined_data["Cxprenylation"] = None

    modification_data = combined_data[
        combined_data["Cxprenylation"].notna()
        ].copy()
    not_modification_data = combined_data[
        combined_data["Cxprenylation"].isna()
        ].copy()
    modification_data["No_prenylation"] = modification_data[
        "No_prenylation"
        ].fillna(1)

    # Filter the data
    #  For modification_data:
    #   Set threshold for "Mass error (ppm)"
    #   Set threshold for "No_prenylation"
    modification_data = modification_data[
        ((modification_data["Mass error (ppm)"] >= - Mass_error_TH)
         & (modification_data["Mass error (ppm)"] <= Mass_error_TH))
        & (modification_data["No_prenylation"] <= No_prenylation_TH)
    ]
    #   Set RT thresholds for each Cxprenylation and extract rows
    #   with RT values exceeding the thresholds
    modification_data = modification_data[
        modification_data.apply(
            lambda row: (
                row["Observed RT (min)"] - row["Ref RT (min)"]
            ) > thresholds[row["Cxprenylation"]], axis=1)
    ]

    #  For not_modification_data:
    #   For each name, extract rows whose time falls within ±RT_TH1
    #   of the "Observed RT (min)" corresponding to
    #   the same "Protein name" in the reference file
    not_modification_data = not_modification_data[
        (not_modification_data["Observed RT (min)"]
         >= not_modification_data["Ref RT (min)"] - RT_TH1) &
        (not_modification_data["Observed RT (min)"]
         <= not_modification_data["Ref RT (min)"] + RT_TH1)
    ]

    # Extract the entry with the largest Response
    # within each "Protein name" and "Sheet_Name" group
    max_response_rows = not_modification_data.groupby(
        ["Protein name", "Sheet_Name"])["Response"].idxmax()
    not_modification_data = not_modification_data.loc[max_response_rows]

    # Calculate %modification
    #  Rename the "Response" column
    not_modification_data = not_modification_data.rename(
        columns={"Response": "Response_unmodified"}
    )
    not_modification_data.to_csv(
        os.path.join(output_path, "not_modification_data.csv"), index=False
    )

    # For rows with the same "Protein name" and "Sheet_Name",
    # calculate the sum of "Response" and store it
    # in the "response_sum" column (combine major and minor products)
    modification_data["response_sum"] = modification_data.groupby(
        ["Protein name", "Sheet_Name"])["Response"].transform("sum")

    # Merge modification_data and filtered_not_modification_data
    # using “Protein name” and “Sheet_Name” as keys
    merged_data = pd.merge(modification_data,
                           not_modification_data[["Protein name", "Sheet_Name",
                                                  "Response_unmodified"]],
                           on=["Protein name", "Sheet_Name"], how="left")
    merged_data["Response_unmodified"] = merged_data[
        "Response_unmodified"].fillna(0)

    return merged_data


# Remove peaks that are clearly noise
def remove_noise(merged_data):
    # Extract the following:
    # - Cases where Response_unmodified is not detected,
    #   Response is smaller than Response_TH, and Charge is an integer
    # - Cases where Response_unmodified is detected,
    #   the pre–noise-removal conversion rate is smaller than modification_TH,
    #   and Charge is an integer
    target_rows = merged_data[
        (((merged_data["Response_unmodified"] != 0)
          & (merged_data["Response"] < (
              merged_data["response_sum"] + merged_data["Response_unmodified"]
          )*modification_TH))
         |
         ((merged_data["Response_unmodified"] == 0)
          & (merged_data["Response"] < Response_TH)))
        & (merged_data["Charge"] == merged_data["Charge"].astype(int))
        ]
    noise_idx = []
    for idx, row in target_rows.iterrows():
        # From experiments using different enzymes/prenyl donors,
        # extract peaks with the same mass
        # whose retention time differs by no more than RT_TH2
        ref_rows = merged_data[
            (merged_data["Protein name"] == row["Protein name"])
            & (merged_data["Cxprenylation"] == row["Cxprenylation"])
            & (merged_data["No_prenylation"] == row["No_prenylation"])
            & (merged_data["Observed RT (min)"]
               >= row["Observed RT (min)"] - RT_TH2)
            & (merged_data["Observed RT (min)"]
               <= row["Observed RT (min)"] + RT_TH2)
            ]
        # Do not consider samples measured multiple times
        ref_rows = ref_rows[~(ref_rows["PT_Name"] == row["PT_Name"])]
        # If no corresponding peak is found, treat it as not noise
        if len(ref_rows) == 0:
            continue
        # If a nearby peak exists with a larger Response, treat it as not noise
        # - Cases where Response_unmodified is not detected and
        #   Response ≥ Response_TH
        elif len(
            ref_rows[
                (ref_rows["Response_unmodified"] != 0)
                & (ref_rows["Response"]
                   >= (ref_rows["response_sum"] +
                       ref_rows["Response_unmodified"]) * modification_TH)
                ]
        ) > 0:
            continue
        # - Cases where Response_unmodified is detected and
        #   the pre–noise-removal conversion rate ≥ modification_TH
        elif len(
            ref_rows[
                (ref_rows["Response_unmodified"] == 0)
                & (ref_rows["Response"] > Response_TH)
                ]
        ) > 0:
            continue
        # If none of the above conditions are met, classify the peak as noise
        else:
            noise_idx.append(idx)

    merged_data_filter = merged_data.drop(noise_idx).reset_index(drop=True)
    merged_data_filter["response_sum"] = merged_data_filter.groupby(
        ["Protein name", "Sheet_Name"])["Response"].transform("sum")

    return merged_data_filter


# Calculate the %modification
def calculate_pct_modification(merged_data_filter,
                               samplesheet_path, output_path):
    # Calculate the %modification
    # (the proportion of each product’s "Response"
    # relative to the substrate + all products)
    merged_data_filter["%modification"] = (
        merged_data_filter["Response"] * 100
        / (merged_data_filter["response_sum"] +
           merged_data_filter["Response_unmodified"])
    )

    merged_data_filter.to_csv(
        os.path.join(output_path, "merged_data.csv"), index=False
    )

    # Calculate the proportion of byproducts
    # Calculate %modification_sum, the total %modification of all products
    merged_data_filter["%modification_sum"] = merged_data_filter.groupby(
        ["Protein name", "Sheet_Name"])["%modification"].transform("sum")

    # Extract the row corresponding to the major product
    main_products = merged_data_filter.loc[merged_data_filter.groupby(
        ["Protein name", "Sheet_Name"])["%modification"].idxmax()]

    # Add a column for the number of replicates
    sample_sheet_list = pd.read_csv(samplesheet_path)
    sample_sheet_list["No_replicate"] = (
        sample_sheet_list.groupby(
            ["Protein_name_list", "PT_Name", "Cxprenylation"]
        )["Sheet_Name"].transform("count")
    )
    main_products = pd.merge(main_products, sample_sheet_list,
                             on=["Sheet_Name", "PT_Name", "Cxprenylation"])

    # Calculate avg_%modification_sum and avg_%modification
    # (means across replicates)
    main_products["sum_%modification_sum"] = main_products.groupby(
        ["Protein name", "PT_Name", "Cxprenylation", "No_prenylation"]
    )["%modification_sum"].transform("sum")
    main_products["sum_%modification"] = main_products.groupby(
        ["Protein name", "PT_Name", "Cxprenylation", "No_prenylation"]
    )["%modification"].transform("sum")
    main_products["avg_%modification_sum"] = (
        main_products["sum_%modification_sum"] / main_products["No_replicate"])
    main_products["avg_%modification"] = (
        main_products["sum_%modification"] / main_products["No_replicate"])

    # Calculate %biproduct
    main_products["avg_%biproduct"] = 100 - 100 * (
        main_products["avg_%modification"] /
        main_products["avg_%modification_sum"])
    final_df = main_products.loc[
        main_products.groupby(
            ["Protein name", "PT_Name", "Cxprenylation"]
        )["%modification_sum"].idxmax()]

    final_df.to_csv(os.path.join(output_path, "final_df.csv"), index=False)

    return final_df


# Create a pivot table (%modification)
def make_modification_pivottable(final_df, protein_names, pt_names,
                                 cxprenylations, PT_name_dict, output_path):
    modification_table = pd.pivot_table(
        final_df,
        values="avg_%modification_sum",
        columns="Protein name",
        index=["PT_Name", "Cxprenylation"],
        fill_value=None
    )
    modification_table = modification_table.reindex(
        columns=protein_names,
        index=pd.MultiIndex.from_product([pt_names, cxprenylations])
    )
    modification_table.index = [
        f"{PT_name_dict[pt]}_{cx}" for pt, cx in modification_table.index
    ]
    modification_table.to_csv(
        os.path.join(output_path, "pct_modification.csv")
    )


# Create a pivot table (%biproduct)
def make_biproduct_pivottable(final_df, protein_names, pt_names,
                              cxprenylations, PT_name_dict, output_path):
    biproduct_table = pd.pivot_table(
        final_df,
        values="avg_%biproduct",
        columns="Protein name",
        index=["PT_Name", "Cxprenylation"],
        fill_value=None
    )
    biproduct_table = biproduct_table.reindex(
        columns=protein_names,
        index=pd.MultiIndex.from_product([pt_names, cxprenylations])
    )
    biproduct_table.index = [f"{PT_name_dict[pt]}_{cx}"
                             for pt, cx in biproduct_table.index]
    biproduct_table.to_csv(os.path.join(output_path, "pct_biproduct.csv"))


# Create a pivot table (No_prenylation)
def make_Noprenylation_pivottable(final_df, protein_names, pt_names,
                                  cxprenylations, PT_name_dict, output_path):
    No_table = pd.pivot_table(
        final_df,
        values="No_prenylation",
        columns="Protein name",
        index=["PT_Name", "Cxprenylation"],
        fill_value=None
    )
    No_table = No_table.reindex(
        columns=protein_names,
        index=pd.MultiIndex.from_product([pt_names, cxprenylations])
    )
    No_table.index = [f"{PT_name_dict[pt]}_{cx}" for pt, cx in No_table.index]
    No_table.to_csv(os.path.join(output_path, "No_prenylation.csv"))


# Create a pivot table (No_replicate)
def make_Noreplicate_pivottable(final_df, protein_names, pt_names,
                                cxprenylations, PT_name_dict, output_path):
    replicate_table = pd.pivot_table(
        final_df,
        values="No_replicate",
        columns="Protein name",
        index=["PT_Name", "Cxprenylation"],
        fill_value=None
    )
    replicate_table = replicate_table.reindex(
        columns=protein_names,
        index=pd.MultiIndex.from_product([pt_names, cxprenylations])
    )
    replicate_table.index = [f"{PT_name_dict[pt]}_{cx}"
                             for pt, cx in replicate_table.index]
    replicate_table.to_csv(os.path.join(output_path, "No_replicate.csv"))


# Make a simple csv file
def edit_final_df(final_df, merged_data_filter, PT_name_dict,
                  Cxprenylation_dict, output_path):
    final_df_simple = final_df[["Sheet_Name", "PT_Name", "Cxprenylation",
                                "Protein name", "avg_%modification_sum",
                                "avg_%biproduct", "No_replicate"]].copy()
    merged_data_filter_simple = merged_data_filter[
        ["Sheet_Name", "PT_Name", "Cxprenylation",
         "Protein name", "Observed RT (min)", "Observed mass (Da)",
         "Expected mass (Da)", "Mass error (ppm)", "%modification"]
        ].copy()
    merged_data_filter_simple = pd.merge(merged_data_filter_simple,
                                         final_df_simple,
                                         on=["Sheet_Name", "PT_Name",
                                             "Cxprenylation", "Protein name"],
                                         how="inner")

    # Convert PT_Name to a categorical type and set its order
    PT_order = list(PT_name_dict.keys())
    cat_type = CategoricalDtype(categories=PT_order, ordered=True)
    merged_data_filter_simple["PT_Name"] = merged_data_filter_simple[
        "PT_Name"].astype(cat_type)

    # Convert Cxprenylation to a categorical type and set its order
    Cxprenylation_order = list(Cxprenylation_dict.keys())
    cat_type2 = CategoricalDtype(categories=Cxprenylation_order, ordered=True)
    merged_data_filter_simple["Cxprenylation"] = merged_data_filter_simple[
        "Cxprenylation"].astype(cat_type2)

    # Sort the table
    merged_data_filter_simple = merged_data_filter_simple.sort_values(
        by=["PT_Name", "Cxprenylation", "Protein name", "%modification"],
        ascending=[True, True, True, False]
    ).reset_index(drop=True)

    # Rename the values in the columns "PT_Name", "Cxprenylation",
    # and "Protein name"
    merged_data_filter_simple["PT_Name"] = merged_data_filter_simple[
        "PT_Name"
        ].map(PT_name_dict)
    merged_data_filter_simple["Cxprenylation"] = merged_data_filter_simple[
        "Cxprenylation"
        ].map(Cxprenylation_dict)
    merged_data_filter_simple["Protein name"] = (
        merged_data_filter_simple["Protein name"]
        .str.replace(r"^PTsub_", "", regex=True)
        .str.replace("_", "", regex=True)
    )
    merged_data_filter_simple = merged_data_filter_simple.drop(
        columns=["Sheet_Name"]
    )
    # Rename the column names
    merged_data_filter_simple = merged_data_filter_simple.rename(
        columns={"PT_Name": "enzyme name",
                 "Cxprenylation": "prenyl donor",
                 "Protein name": "peptide name",
                 "Observed RT (min)": "observed retention time (min)",
                 "Observed mass (Da)": "observed mass (Da)",
                 "Expected mass (Da)": "expected mass (Da)",
                 "Mass error (ppm)": "mass error (ppm)",
                 "avg_%modification_sum": "LC-MS conversion (%)"}
    )
    merged_data_filter_simple.to_csv(
        os.path.join(output_path, "final_df_simple.csv"), index=False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input",
        default=["./data/extracted_data/all_*knownPT_data.csv"], nargs="*",
        help="Path to the CSV files of the extracted data"
    )
    parser.add_argument(
        "-o", "--output", default="./data/results/",
        help="Path to the output directory for the processed data"
    )
    parser.add_argument(
        "-l", "--list", default="./data/protein_pt_list.csv",
        help="CSV file path for prenyltransferase and peptide names"
    )
    parser.add_argument(
        "-s", "--samplesheet",
        default="./data/extracted_data/sample_sheet_list.csv",
        help="Path to sample_sheet_list.csv"
    )
    parser.add_argument(
        "-r", "--rename", default="",
        help=(
            "CSV file path containing the correspondence table"
            " of prenyltransferase names"
        )
    )
    parser.add_argument(
        "-e", "--edit", action="store_true",
        help="Decide whether to make a simple csv file"
    )
    args = parser.parse_args()

    # Combine all data
    csv_files = []
    for x in args.input:
        csv_files += glob.glob(x)
    print(csv_files)
    csv_data = []
    for csv_file in csv_files:
        data = pd.read_csv(csv_file)
        csv_data.append(data)
    combined_data = pd.concat(csv_data, ignore_index=True)

    # Make output directory
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    # Make list of "Protein name", "PT_Name", and "Cxprenylation"
    protein_pt_list = pd.read_csv(args.list)
    protein_names = protein_pt_list["Protein name"].dropna().unique().tolist()
    pt_names = protein_pt_list["PT_Name"].dropna().unique().tolist()
    cxprenylations = list(Cxprenylation_dict.keys())

    # Make dictionary to rename the values in the columns "PT_Name"
    if args.rename == "":
        PT_name_dict = {x: x for x in pt_names}
    else:
        PT_name_list = pd.read_csv(args.rename)
        PT_name_dict = PT_name_list.set_index("PT_No")["PT_name"].to_dict()
    print(PT_name_dict)

    # Analyze data
    merged_data = modify_data(combined_data, args.output)
    merged_data_filter = remove_noise(merged_data)
    final_df = calculate_pct_modification(
        merged_data_filter, args.samplesheet, args.output
    )

    make_modification_pivottable(
        final_df, protein_names, pt_names, cxprenylations,
        PT_name_dict, args.output
    )
    make_biproduct_pivottable(
        final_df, protein_names, pt_names, cxprenylations,
        PT_name_dict, args.output
    )
    make_Noprenylation_pivottable(
        final_df, protein_names, pt_names, cxprenylations,
        PT_name_dict, args.output
    )
    make_Noreplicate_pivottable(
        final_df, protein_names, pt_names, cxprenylations,
        PT_name_dict, args.output
    )

    # Make a simple csv file
    if args.edit:
        edit_final_df(
            final_df, merged_data_filter, PT_name_dict,
            Cxprenylation_dict, args.output
        )

    print("Completed!")

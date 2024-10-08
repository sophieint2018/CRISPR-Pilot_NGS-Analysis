import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob

overall_reads_data = {}
overall_percentages_data = {}
"""
Includes the dfs to export to excel with the filtered data.
overall_data = {
    "Atm" = {
        "A1" = [WT, InFrame, FS]
        "A2" = ...
        .
        .
        .
        "A7" =

    }
    "Cdkn2a" = ...
    .
    .
    .
    "Trp53" = ...
}
"""

gRNAs = ["Trp53"]
important_data = ["TargetSequence", "Reads", "Type", "Pct", "IndelLength"]
threshold = 0.2
num_samples = 6

sample_names = []


def convert_excel_to_pd(path: str, sheet_name: str):
    """
    Converts excel documents into a panda dataframe and returns the frame.
    """
    xls = pd.ExcelFile(path)
    df = pd.read_excel(xls, sheet_name)
    df.columns = df.columns.map(str)
    return df


def filter_data(df: pd):
    """
    This method will take in the dataframe of data from Genewiz and filter
    out columns of data that are not of interest (based on the global
    variable called important_data) and then exclude data of the following
    characteristics:
    - reads with a pct below 0.2
    - reads that were classified as Base Changes by Genewiz

    Additionally, if the Type of read is an Insertion, Deletion, or both, 
    the method will fill in a new column called "Classification" that
    defines the indel as InFrame or FS (Frame Shift). This data is made
    into a dataframe that is returned by this method. 
    """
    # Removes unwanted columns from the dataframe
    df_filtered = df[important_data]

    # Remove rows with Pcts lower than the decided threshold
    df_pct_filtered = df_filtered[df_filtered['Pct'] >= threshold]

    # Isolate the "Type" column
    df_column = df_pct_filtered["Type"]
    # Creates an identical column with True if "Base Change" is not in the Type Column
    df_column = ~df_column.str.contains("Base Change")
    # Remove columns that are labeled as Base Changes
    df_basechange_filtered = df_pct_filtered[df_column]

    # Make a new df
    final_df = df_basechange_filtered.copy()

    # Find the row numbers of the Insertion/Deletion reads
    rows_to_check = []
    mask = final_df["Type"] == "Insertion"
    rows_to_check.extend(final_df.index[mask].tolist())
    mask = final_df["Type"] == "Deletion"
    rows_to_check.extend(final_df.index[mask].tolist())
    mask = final_df["Type"] == "Insertion and Deletion"
    rows_to_check.extend(final_df.index[mask].tolist())
    rows_to_check.sort()

    # Add the empty Classification column
    final_df.loc[:, 'Classification'] = np.nan

    # Add classification (InFrame, FS) to the "Classification" Column
    for index in rows_to_check:
        indel_length = final_df.loc[index, "IndelLength"]
        if (indel_length % 3 == 0):
            final_df.loc[index, "Classification"] = "InFrame"
        else:
            final_df.loc[index, "Classification"] = "FS"

    # Adding WT to Classification
    mask = final_df["Type"] == "WT"
    rows_to_check = (final_df.index[mask].tolist())
    for index in rows_to_check:
        final_df.loc[index, "Classification"] = "WT"

    # Check for NaN values
    has_nan = final_df.isna().any().any()
    if (has_nan):
        print("The Type of reads classified by Genewiz is more than Insertion, Deletion, and WT")
        print("Modify the script to account for those values")
        print(final_df)
        return

    return (final_df)


def analyze_df(df: pd):
    """
    This method takes in data from the filtered dataframe and analyzes the data for reads
    and percentage of data that is WT, InFrame, and FS. These are returned as a tuple:
    ( [WT #reads, InFrame #reads, FS #reads],
      [WT percentage, InFrame percentage, FS percentage])
    """
    total_reads = df["Reads"].sum()

    # Add up reads for WT
    WT_reads = 0
    mask = df["Classification"] == "WT"
    rows_to_check = (df.index[mask].tolist())
    for index in rows_to_check:
        WT_reads += df.loc[index, "Reads"]

    # Add up reads for FS
    FS_reads = 0
    mask = df["Classification"] == "FS"
    rows_to_check = (df.index[mask].tolist())
    for index in rows_to_check:
        FS_reads += df.loc[index, "Reads"]

    # Add up reads for InFrame
    InFrame_reads = 0
    mask = df["Classification"] == "InFrame"
    rows_to_check = (df.index[mask].tolist())
    for index in rows_to_check:
        InFrame_reads += df.loc[index, "Reads"]

    # Check that the total number of reads is equivalent to the number of WT, InFrame, and FS reads
    if (total_reads != (WT_reads + InFrame_reads + FS_reads)):
        print(
            "The total number of reads does not equal the sum of WT, InFrame, and FS reads.")
        print("Something was likely classified incorrectly or another error occurred")
        return

    # Calculate final analytics (total reads and percentages)
    read_data = [WT_reads, InFrame_reads, FS_reads]
    percentage_data = [(WT_reads/total_reads)*100,
                       (InFrame_reads/total_reads)*100, (FS_reads/total_reads)*100]

    data = [read_data, percentage_data]
    return data


def create_dfs(df: pd):
    """
    This method takes the data stored in the dataframe and outputs a
    summary spreadsheet that can go into the CRISPRCutting.py script
    for graph analysis. 
    """

    print("Creating the data table for the dataframe")

    # Creating lists for each sample of all the data appended together (gRNA alphabetical order)
    final_data_dict = {}
    final_data_dict = {key: [] for key in sample_names}
    for gRNA in gRNAs:
        dict_with_samples = df[gRNA]
        for sample in sample_names:
            existing_list = final_data_dict[sample]
            existing_list.extend(dict_with_samples[sample])
            final_data_dict[sample] = existing_list

    # Create a list for the column names of the dataframe
    column_names = []
    for gRNA in gRNAs:
        column_names.append("{} WT".format(gRNA))
        column_names.append("{} InFrame".format(gRNA))
        column_names.append("{} FS".format(gRNA))

    length = len(column_names)

    # Check lengths and create final list to import to dataframe
    final_data_list = []
    for sample in sample_names:
        sample_data = final_data_dict[sample]
        if (len(sample_data) == length):
            final_data_list.append(sample_data)
        else:
            print("Error. The amount of data for each sample is unequal.")
            print(sample_data)
            return

    print("Creating the df")
    final_df = pd.DataFrame(final_data_list, columns=column_names)

    # Adding in a first column that labels the rows of the first column
    # with Genewiz sample number
    final_df.insert(0, "Genewiz Sample", sample_names)

    # Adding a %Mutant column for each gRNA
    for gRNA in gRNAs:
        final_df["{} Mutant".format(gRNA)] = final_df["{} InFrame".format(
            gRNA)] + final_df["{} FS".format(gRNA)]
    print(final_df)
    return final_df


def main(current_path: str, data_path: str, sheet_names: list):
    global overall_reads_data
    global overall_percentages_data

    # Path for the output data
    filter_excel_sheet_path = current_path + "/Analyzed Data"

    # Find the excel files to analyze, analyze for each excel sheet and
    # page in the excel sheet, export final data
    paths_to_data_files = glob.glob(os.path.join(data_path, "*.xlsx"))
    print(paths_to_data_files)
    for file in paths_to_data_files:
        gRNA = file.split("/")[-1].split(" ")[0]
        print("Analyzing {} data:".format(gRNA))

        gRNA_reads_dic = {}
        gRNA_percentages_dic = {}

        for sheet in sheet_names:
            print("{} sheet".format(sheet))
            sample = sheet.split("-")[-1]

            # Filter the excel sheet down
            df = convert_excel_to_pd(file, sheet)
            df = filter_data(df)
            export_excel_name = filter_excel_sheet_path + \
                "/{} {}.xlsx".format(gRNA, sample)
            df.to_excel(export_excel_name, index=False)

            # Get overall data values from the filtered dataframe
            excel_data = analyze_df(df)
            gRNA_reads_dic[sample] = excel_data[0]
            gRNA_percentages_dic[sample] = excel_data[1]
        overall_reads_data[gRNA] = gRNA_reads_dic
        overall_percentages_data[gRNA] = gRNA_percentages_dic

    # Create dfs and export
    print("Creating reads dataframe")
    read_df_to_export = create_dfs(overall_reads_data)
    print("Exporting reads dataframe to excel")
    read_df_to_export.to_excel(
        current_path + "/CRISPRCutting Results/CRISPRCutting_Reads.xlsx", index=False)

    print("Creating percentage dataframe")
    percentage_df_to_export = create_dfs(overall_percentages_data)
    print("Exporting percentage dataframe to excel")
    percentage_df_to_export.to_excel(
        current_path + "/CRISPRCutting Results/CRISPRCutting_Percentages.xlsx", index=False)


if __name__ == "__main__":
    current_path = os.getcwd()
    path_to_data = current_path + "/Raw Data"

    # Create sheet names
    sheet_names = []
    for num in range(1, num_samples + 1):
        sheet_names.append("Plate1-A{}".format(num, num))
        sample_names.append("A{}".format(num))

    main(current_path, path_to_data, sheet_names)

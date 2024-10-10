import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

df = pd.DataFrame()
gRNAs = ["Trp53"]
sample_dict = {}
color_palette = ["#b30000", "#4421af", "#0d88e6",
                 "#00b7c7", "#5ad45a", "#8be04e", "#ebdc78"]


def convert_excel_to_pd(file_name: str):
    """
    Converts excel documents into a panda dataframe.
    """
    global df
    df = pd.read_excel(file_name)
    df.columns = df.columns.map(str)
    # print(df)


def calculate_mutant_wt_percentage():
    """
    Method will take the global list of gRNAs and look for
    the %Mutant and %FS values for each gRNA and populate
    new columns of data for "%In-frame" and "%WT" for each
    sample.

    It will also check that the sum equals %100 before returning

    e.g.
    Samples     Atm %Mutant     Atm %FS    Atm %In-frame    Atm %WT
    """

    # Calculates the In-Frame and WT columns
    for gRNA in gRNAs:
        df["{} %InFrame".format(gRNA)] = df["{} %Mutant".format(gRNA)
                                            ] - df["{} %FS".format(gRNA)]
        df["{} %WT".format(gRNA)] = 100 - df["{} %Mutant".format(gRNA)
                                             ]

    # Double checks that the column math was done correctly
    check = True
    for gRNA in gRNAs:
        if not all(df["{} %InFrame".format(gRNA)] + df["{} %FS".format(gRNA)] + df["{} %WT".format(gRNA)] == 100):
            check = False

    if check is False:
        raise Exception(
            "PD Column Math Failed. Double check the code in the calculate_mutant_wt_percentage() method")

    # print(df)


def create_dict_for_each_sample_data():
    """
    This method will iterate across the dataframe and create a master dictionary of data that
    will be used to create the stacked bar graph.

    sample_dict = {
        A1 (203) = {
            "WT" = [Atm, Cdkn2a, Ppm1d, Pten, Trp53]
            "InFrame" = [Atm, Cdkn2a, Ppm1d, Pten, Trp53]
            "FS" = [Atm, Cdkn2a, Ppm1d, Pten, Trp53]
        }
    }

    """
    # Adds keys to the sample dictionary based on the data in the Genewiz Sample Column
    global sample_dict
    sample_dict = {key: None for key in df["Genewiz Sample"]}
    for key in sample_dict:

        # Gets the row number for the sample
        row_number = df.index.get_loc(df[df["Genewiz Sample"] == key].index[0])

        # Creates empty arrays for each data type for the stacked bar graph
        WT = np.zeros(len(gRNAs))
        InFrame = np.zeros(len(gRNAs))
        FS = np.zeros(len(gRNAs))

        # Iterates across all possible gRNAs and adds the WT, InFrame, FS data to the arrays
        for gRNA in gRNAs:
            WT[gRNAs.index(gRNA)] = df.loc[row_number, "{} WT".format(gRNA)]
            InFrame[[gRNAs.index(gRNA)]] = df.loc[row_number,
                                                  "{} InFrame".format(gRNA)]
            FS[gRNAs.index(gRNA)] = df.loc[row_number, "{} FS".format(gRNA)]

        # Create a temporary dictionary that holds all the data together
        temp_dict = {}
        temp_dict["WT"] = WT
        temp_dict["InFrame"] = InFrame
        temp_dict["FS"] = FS
        # temp_dict["Genotype"] = df.loc[row_number, "Genotype"]

        # Set each sample key's value to the temporary dictionary
        sample_dict[key] = temp_dict


def create_stacked_bar_graph():
    """
    This method generates a separate stacked bar graph of all the data in sample_dict.

    Figures are saved to a subfolder titled "Figures"

    """
    WT = []
    InFrame = []
    FS = []
    for key in sample_dict.keys():
        sample = sample_dict[key]
        # genotype = sample["Genotype"]
        WT.extend(sample["WT"])
        InFrame.extend(sample["InFrame"])
        FS.extend(sample["FS"])

    WT_array = np.array(WT)
    InFrame_array = np.array(InFrame)
    FS_array = np.array(FS)
    sample_label = list(range(1, len(WT) + 1))
    sample_label_array = np.array(sample_label)

    # Create the graph
    plt.figure(figsize=(10, 6))
    plt.bar(sample_label_array, WT_array, label="WT", color=color_palette[0])
    plt.bar(sample_label_array, InFrame_array, bottom=WT_array,
            label="In-Frame", color=color_palette[1])
    plt.bar(sample_label_array, FS_array, bottom=WT_array + InFrame_array,
            label="FS", color=color_palette[2])
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xlabel("Sample Number")
    plt.ylabel("Percent of Reads")

    plt.ylim(0, 100 + 5)

    plt.tight_layout()
    plt.subplots_adjust(right=0.8)
    plt.savefig("Figures/CRISPRCuttingAnalysis.png".format(key), dpi=600)

    plt.show()


def main(file_name: str):
    convert_excel_to_pd(file_name)
    # calculate_mutant_wt_percentage()
    create_dict_for_each_sample_data()
    create_stacked_bar_graph()


if __name__ == "__main__":
    current_path = os.getcwd()
    path_to_data = current_path + "/CRISPRCutting Results/CRISPRCutting_Percentages.xlsx"
    main(path_to_data)

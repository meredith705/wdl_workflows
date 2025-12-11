import argparse
import pandas as pd
import numpy as np
import pysam
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import datetime
import json
from pathlib import Path

labelsize=10
ticksize=9

sns.set(
    style="white",  # sets the background to white
    rc={
        "xtick.labelsize": labelsize,
        "ytick.labelsize": labelsize,
        "axes.labelsize": labelsize,
        "axes.titlesize": ticksize,
        "figure.facecolor": "white"  # ensures figure background is white
    }
)

def plot_data(df):
    """ FN: hifiasm only
        FP: shasta only. 

        Scatter plots:
      Panel 1: TP-base vs TP-comp
      Panel 2: FP vs FN
      Panel 3: precision vs recall
    """
    figsize=(8, 12)

    sns.set(style="whitegrid")

    fig, axs = plt.subplots(3, 1, figsize=figsize)

    # ----------------------
    # Panel 1: TP-base vs TP-comp
    # ----------------------
    sns.scatterplot(
        data=df,
        x="TP-base",
        y="TP-comp",
        ax=axs[0],
        s=60,
        edgecolor="k"
    )
    axs[0].set_title("TP-base vs TP-comp")
    axs[0].set_xlabel("TP-base")
    axs[0].set_ylabel("TP-comp")

    # label points
    for _, r in df.iterrows():
        axs[0].text(r["TP-base"], r["TP-comp"], r["sample"], fontsize=9)

    # ----------------------
    # Panel 2: FP vs FN
    # ----------------------
    sns.scatterplot(
        data=df,
        x="FP",
        y="FN",
        ax=axs[1],
        s=60,
        edgecolor="k"
    )
    axs[1].set_title("FP vs FN")
    axs[1].set_xlabel("Shasta only SVs")
    axs[1].set_ylabel("Hifiasm only SVs")

    for _, r in df.iterrows():
        axs[1].text(r["FP"], r["FN"], r["sample"], fontsize=9)

    # ----------------------
    # Panel 3: precision vs recall
    # ----------------------
    sns.scatterplot(
        data=df,
        x="precision",
        y="recall",
        ax=axs[2],
        s=60,
        edgecolor="k"
    )
    axs[2].set_title("Precision vs Recall")
    axs[2].set_xlabel("precision")
    axs[2].set_ylabel("recall")

    for _, r in df.iterrows():
        axs[2].text(r["precision"], r["recall"], r["sample"], fontsize=9)

    plt.tight_layout()
    plt.savefig(f"RUSH_truvariBench_summary.png",dpi=300, facecolor='white', transparent=False)




def load_jsons(dirs):


    # fix this to plot from dirs 
    # dirs = ["path/to/dir1", "path/to/dir2"]  # directories with summary.json

    metrics_list = []

    for d in dirs:
        json_path = Path(d) / "summary.json"
        if not json_path.exists():
            continue

        with open(json_path) as f:
            data = json.load(f)

        # print(json.dumps(data, indent=4))

        # get TP-base, TP-comp, FP, FN, precision, recall, f1
        metrics = {}

        # Extract values
        metrics["TP-base"] = data.get("TP-base", 0)
        metrics["TP-comp"] = data.get("TP-comp", 0)
        metrics["FP"] = data.get("FP", 0)
        metrics["FN"] = data.get("FN", 0)

        # Precision, recall, f1 may be under 'stats' or 'metrics'
        metrics["precision"] = data.get("precision", 0)
        metrics["recall"] = data.get("recall", 0)
        metrics["f1"] = data.get("f1", 0)

        # Store sample/assembler info
        metrics["sample"] = "RUSH_"+Path(d).name.split("_")[1]

        metrics_list.append(metrics)

    # Combine into a DataFrame
    df_metrics = pd.DataFrame(metrics_list)

    print(df_metrics.head())

    plot_data(df_metrics)


if __name__ == "__main__":

    # Make an argument parser
    parser = argparse.ArgumentParser(description="Plot a summary.json from Truvari Bench from a list of directories of output.")

    # add argument for the input vcf file
    parser.add_argument(
        "-d","--dirs",
        nargs="+",
        help="List of directories containing summary.json files"
    )

    # parser.add_argument(
    #     "-f","--hifiasm_count_tsv",
    #     type=str,
    #     required=True,
    #     help="Path to the input file to be plotted."
    # )


    # parser.add_argument(
    #     "-o","--output_directory",
    #     type=str,
    #     required=True,
    #     help="Path to the input qtl tsv to be analyzed."
    # )

    args = parser.parse_args()

    load_jsons(args.dirs)



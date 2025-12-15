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
import gzip

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

def plot_total_bp(df):
    sns.set(style="whitegrid")

    plt.figure(figsize=(10,6))
    sns.barplot(
        data=df,
        x="sample",
        y="total_bp",
        hue="assembler",
        palette="Set2",
        edgecolor="k"
    )

    plt.ylabel("Total bp (INS + DEL)")
    plt.title("Total inserted + deleted bp per assembler")
    plt.tight_layout()
    plt.show()

def open_vcf(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def sum_ins_del(vcf_path):
    total_ins = 0
    total_del = 0

    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            ref = fields[3]
            alts = fields[4].split(",")

            for alt in alts:
                # skip symbolic calls like <INS>, <DEL>, <DUP>, etc
                if alt.startswith("<") and alt.endswith(">"):
                    continue

                # compute length delta
                d = len(alt) - len(ref)

                if d > 0:
                    total_ins += d
                elif d < 0:
                    total_del += abs(d)

    return total_ins, total_del

def count_bps_fn_fp(dirs):

    indel_len_list = []

    for d in dirs:

        sample = "RUSH_"+Path(d).name.split("_")[1]

        shastaOnly = Path(d) / "fp.vcf.gz"
        if not shastaOnly.exists():
            continue

        shasta_ins, shasta_del = sum_ins_del(shastaOnly)

        hifiasmOnly = Path(d) / "fn.vcf.gz"
        if not hifiasmOnly.exists():
            continue

        hifiasm_ins, hifiasm_del = sum_ins_del(hifiasmOnly)

        indel_len_list.append({
            "sample": sample,
            "assembler": "shasta",
            "ins_bp": shasta_ins,
            "del_bp": shasta_del,
            "total_bp": shasta_ins + shasta_del
        })

        indel_len_list.append({
            "sample": sample,
            "assembler": "hifiasm",
            "ins_bp": hifiasm_ins,
            "del_bp": hifiasm_del,
            "total_bp": hifiasm_ins + hifiasm_del
        })

        df = pd.DataFrame(indel_len_list)

        plot_total_bp(df)



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
    df_metrics.to_csv("RUSH_HifiasmShasta_truvariBench_summary.csv")

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

    args = parser.parse_args()

    load_jsons(args.dirs)


    count_bps_fn_fp(args.dirs)



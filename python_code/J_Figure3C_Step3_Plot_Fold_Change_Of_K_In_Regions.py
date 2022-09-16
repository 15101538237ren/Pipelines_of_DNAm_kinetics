# -*- coding:utf-8 -*-
import os, argparse, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Plot Figure 3c: Fold Change of K in Different Region Context')

parser.add_argument('-k_fp', type=str, default="../data/K-500bp-CpGd_HMM_DNase_NOS_Rmvd_Reversible_Sites.bed", help= 'The file path where the bed file of Ks is stored')

parser.add_argument('-region_dir', type=str, default="../data/K_region_intersect", help= 'The directory where the bed file with K and region intersected')

parser.add_argument('-region_type', type=str, default="TFBS", choices=["TFBS", "Histone_Modification", "Genomic_Regions", "CGI"], help='The region type, choices are: Genomic_Regions, Histone_Modification, TFBS, CGI')

parser.add_argument('-min_meth', type=float, default=0.5,
                    help='The minimal WGBS methylation level for CpG filtering, default: 0')

parser.add_argument('-max_meth', type=float, default=1.01,
                    help='The maximum WGBS methylation level for CpG filtering, default: 1.01')

parser.add_argument('-meth_category', type=str, default="IM", choices=["IM", "HM", "GM"],
                    help='The category of K values')

parser.add_argument('-k_category', type=str, default="slow", choices=["slow", "medium", "fast", "all"],
                    help='The category of K values')

parser.add_argument('-fig_path', type=str, default="../figures/J_Figures/fig3c_fold_change.svg", help='The path of figure3c output')


file_ordered_names_dict = {
            "Genomic_Regions": ["Genome","DNase", "Enhancer", "CGI", "Promoter",
                                "TFBR", "Exons", "Introns", "Intergenic", "5UTR", "3UTR",
                                "H2AZ", "LINE", "SINE", "LTR"],#
            "Histone_Modification": ["H3k4me1Not3_no_k27me3_no_k27ac", "H3k4me1Not3_k27me3", "H3k4me1Not3_k27ac",
                                                    "H3k4me3Not1_no_k27me3_no_k27ac", "H3k4me3Not1_k27me3", "H3k4me3Not1_k27ac",
                                                    "H3K9Me3", "H3K36me3"],
            "TFBS" : ["DNase", "CTCF", "EP300", "EZH2", "H2AZ", "tfbs_cluster_v3"]
            }
file_labels_dict = {
            "Genomic_Regions": ["Genome", "DNase", "Enhancer", "CGI", "Promoter",
                                "TFBR", "Exons", "Introns", "Intergenic", "5UTR", "3UTR",
                                "H2AZ", "LINE", "SINE", "LTR"],
            "Histone_Modification": ["H3k4me1 Solo", "H3k4me1 + H3k27me3", "H3k4me1 + H3k27ac",
                                        "H3k4me3 Solo", "H3k4me3 + H3k27me3", "H3k4me3 + H3k27ac", "H3K9Me3", "H3K36me3"],
            "TFBS": ["DNase", "CTCF", "EP300", "EZH2", "H2AZ", "General\nTFBSs"]
                }

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def get_meth_indexs(df, min_meth, max_meth):
    meth_inds = [mid for mid, meth in enumerate(df[:, 3]) if
                 (meth != "." and ((float(meth) > min_meth) and float(meth) < max_meth))]
    return df[meth_inds, :]

def get_regional_percentage(regional_k_fp, meth_category, min_meth, max_meth, min_k, max_k):
    k_sub_category_df = pd.read_csv(regional_k_fp, sep="\t", header=None).values
    k_sub_category_df = get_meth_indexs(k_sub_category_df, min_meth, max_meth)
    k_number = k_sub_category_df.shape[0]
    if meth_category == "IM":
        min_meth = 0.5
        max_meth = 0.8
        k_sub_category_df = get_meth_indexs(k_sub_category_df, min_meth, max_meth)
    elif meth_category == "HM":
        min_meth = 0.8
        max_meth = 1.01
        k_sub_category_df = get_meth_indexs(k_sub_category_df, min_meth, max_meth)
    ks_in_k_sub_category_df = k_sub_category_df[:, 4].astype(float)
    k_sub_category_df = k_sub_category_df[np.logical_and(ks_in_k_sub_category_df > min_k, ks_in_k_sub_category_df < max_k), :]
    sub_percentage_of_k = k_sub_category_df.shape[0] / k_number
    return sub_percentage_of_k

def plot_fold_change(k_fp, region_dir, region_type, k_category, meth_category, min_meth, max_meth, fig_path):
    file_ordered_names = file_ordered_names_dict[region_type]
    file_labels = file_labels_dict[region_type]
    n_files = len(file_labels)
    ind = np.arange(n_files)
    percentages = np.zeros((n_files, 1))
    fold_changes = np.zeros((n_files, 1))

    print("Processing Ks")
    k_sub_category_df = pd.read_csv(k_fp, sep="\t", header=None).values
    k_sub_category_df = get_meth_indexs(k_sub_category_df, min_meth, max_meth)

    ks_in_k_sub_category_df = k_sub_category_df[:, 4].astype(float)
    slow_k_pct = np.percentile(ks_in_k_sub_category_df, 33)
    mid_k_pct = np.percentile(ks_in_k_sub_category_df, 67)
    min_k = 0
    max_k = 50
    if k_category == "slow":
        max_k = slow_k_pct
    elif k_category == "medium":
        min_k = slow_k_pct
        max_k = mid_k_pct
    elif k_category == "fast":
        min_k = mid_k_pct
    else:
        pass
    sub_percentage_of_k = get_regional_percentage(k_fp, meth_category, min_meth, max_meth, min_k, max_k)
    for cid, file_label in enumerate(file_labels):
        print("Processing %s" % file_label)
        regional_k_fp = os.path.join(region_dir, region_type, "%s.bed" % file_ordered_names[cid])
        percentages[cid] = get_regional_percentage(regional_k_fp, meth_category, min_meth, max_meth, min_k, max_k)
        fold_changes[cid] = percentages[cid]/sub_percentage_of_k

    plt.rc('font', weight='bold')
    plt.rc('xtick', labelsize=10)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    plt.subplots_adjust(wspace=0.4, hspace=0.5, bottom=0.2)
    ax.bar(ind, tuple(list(fold_changes[:, 0])), 0.38, color="blue")
    for cid, file_label in enumerate(file_labels):
        fold_change = fold_changes[cid, 0]
        ax.text(cid - 0.1, fold_change + 0.3, "%.1f" % fold_change, color="black", weight='bold', fontsize=12)
    ax.set_ylabel('Fold Enrichment Over All ks', weight='bold', fontsize=14)
    ax.set_xticks(ind)
    ax.set_xticklabels(file_labels)
    ax.set_ylim([0, fold_changes.max() + 1])
    # ax.tick_params(axis='x', rotation=90)
    mkdir(os.path.dirname(fig_path))
    plt.savefig(fig_path, dpi=300)

if __name__ == "__main__":
    args = parser.parse_args()
    k_fp = args.k_fp
    region_dir = args.region_dir
    region_type = args.region_type
    k_category = args.k_category
    meth_category = args.meth_category
    min_meth = args.min_meth
    max_meth = args.max_meth
    fig_path = args.fig_path
    plot_fold_change(k_fp, region_dir, region_type, k_category, meth_category, min_meth, max_meth, fig_path)
# -*- coding:utf-8 -*-
import os, argparse, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Plot Figure 3c: Chrom HMM Distribution')

parser.add_argument('-chromHMM_fp', type=str, default="../data/HmmH1hescHMM.bed", help= 'The file path where the annotated ChromHMM bed file for H1-ESC is stored')

parser.add_argument('-k_fp', type=str, default="../data/K-500bp-CpGd_HMM_DNase_NOS_Rmvd_Reversible_Sites.bed", help= 'The file path where the bed file of Ks is stored')

parser.add_argument('-k_in_tfbs_cluster_fp', type=str, default="../data/TFBS_clusters/wgEncodeRegTfbsClusteredV3_K.bed", help='The file path where K intersected with wgEncodeRegTfbsClusteredV3.bed stored')

parser.add_argument('-k_in_tfbs_jsd_over_95_pct_fp', type=str, default="../data/TFBS_clusters/General_TFs_JSD_Over_95_pct_K.bed", help='The file path where K intersected with General_TFs_JSD_Over_95_pct.bed stored')

parser.add_argument('-min_meth', type=float, default=0.5,
                    help='The minimal WGBS methylation level for CpG filtering, default: 0')

parser.add_argument('-max_meth', type=float, default=1.01,
                    help='The maximum WGBS methylation level for CpG filtering, default: 1.01')

parser.add_argument('-k_category', type=str, default="slow", choices=["slow", "medium", "fast", "all"],
                    help='The category of K values')

parser.add_argument('-fig_path', type=str, default="../figures/J_Figures/fig3c.svg", help='The path of figure3c output')


def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def get_meth_indexs(df, min_meth, max_meth):
    meth_inds = [mid for mid, meth in enumerate(df[:, 3]) if
                 (meth != "." and ((float(meth) > min_meth) and float(meth) < max_meth))]
    return df[meth_inds, :]
def plot_chrmoHMM(chromHMM_fp, k_fp, k_category, min_meth, max_meth, fig_path, only_enhancer_in_fold_change = True):
    chrmoHMM_labels = "Active Promoter", "Weak Promoter", "Poised Promoter", "Strong Enhancer", "Strong Enhancer", "Weak Enhancer", "Weak Enhancer", "Insulator", "Txn Transition", "Txn Elongation", "Weak Txn", "Repressed", "Heterochrom/lo", "Repetitive/CNV", "Repetitive/CNV"
    cm = plt.get_cmap('gist_rainbow')
    category_labels = ["H1\n(ESC)", "All\nrates", "%s\nrates" % k_category.capitalize(), "%s rates\nin TFBS" % k_category.capitalize(), "%s rates\nin TFBS\n>95%% JSD" % k_category.capitalize()]
    ind = np.arange(len(category_labels))
    chrHMMs = np.arange(len(chrmoHMM_labels))
    class_data = np.zeros((len(chrHMMs), len(ind)))

    chromHMM_df = pd.read_csv(chromHMM_fp, sep="\t", header=None).values
    k_df = pd.read_csv(k_fp, sep="\t", header=None).values
    k_sub_category_df = pd.read_csv(k_fp, sep="\t", header=None).values
    k_in_tfbs_cluster_df = pd.read_csv(k_in_tfbs_cluster_fp, sep="\t", header=None).values
    k_in_tfbs_jsd_over_95_pct_df = pd.read_csv(k_in_tfbs_jsd_over_95_pct_fp, sep="\t", header=None).values

    k_df = get_meth_indexs(k_df, min_meth, max_meth)
    k_sub_category_df = get_meth_indexs(k_sub_category_df, min_meth, max_meth)
    k_in_tfbs_cluster_df = get_meth_indexs(k_in_tfbs_cluster_df, min_meth, max_meth)
    k_in_tfbs_jsd_over_95_pct_df = get_meth_indexs(k_in_tfbs_jsd_over_95_pct_df, min_meth, max_meth)

    ks_in_k_sub_category_df = k_sub_category_df[:, 4].astype(float)
    slow_k_pct = np.percentile(ks_in_k_sub_category_df, 33)
    mid_k_pct = np.percentile(ks_in_k_sub_category_df, 67)
    if k_category == "slow":
        k_sub_category_df = k_sub_category_df[ks_in_k_sub_category_df < slow_k_pct, :]
        k_in_tfbs_cluster_df = k_in_tfbs_cluster_df[k_in_tfbs_cluster_df[:, 4].astype(float) < slow_k_pct, :]
        k_in_tfbs_jsd_over_95_pct_df = k_in_tfbs_jsd_over_95_pct_df[k_in_tfbs_jsd_over_95_pct_df[:, 4].astype(float) < slow_k_pct, :]
    elif k_category == "medium":
        k_sub_category_df = k_sub_category_df[np.logical_and(ks_in_k_sub_category_df > slow_k_pct, ks_in_k_sub_category_df < mid_k_pct), :]
        k_in_tfbs_cluster_df = k_in_tfbs_cluster_df[np.logical_and(k_in_tfbs_cluster_df[:, 4].astype(float) > slow_k_pct, k_in_tfbs_cluster_df[:, 4].astype(float) < mid_k_pct), :]
        k_in_tfbs_jsd_over_95_pct_df = k_in_tfbs_jsd_over_95_pct_df[np.logical_and(k_in_tfbs_jsd_over_95_pct_df[:, 4].astype(float) > slow_k_pct, k_in_tfbs_jsd_over_95_pct_df[:, 4].astype(float) < mid_k_pct), :]
    elif k_category == "fast":
        k_sub_category_df = k_sub_category_df[ks_in_k_sub_category_df > mid_k_pct, :]
        k_in_tfbs_cluster_df = k_in_tfbs_cluster_df[k_in_tfbs_cluster_df[:, 4].astype(float) > mid_k_pct, :]
        k_in_tfbs_jsd_over_95_pct_df = k_in_tfbs_jsd_over_95_pct_df[k_in_tfbs_jsd_over_95_pct_df[:, 4].astype(float) > mid_k_pct, :]
    else:
        pass
    number_of_ks = []
    for cid, class_label in enumerate(category_labels):
        print("Processing %s" % class_label)
        if cid == 0:
            categories = chromHMM_df
            classes = categories[:, -1].astype(int)
        else:
            if cid == 1:
                categories = k_df
            elif cid == 2:
                categories = k_sub_category_df
            elif cid == 3:
                categories = k_in_tfbs_cluster_df
            else:
                categories = k_in_tfbs_jsd_over_95_pct_df
            classes = categories[:, -3]
            classes = classes[classes != "."].astype(int)
        number_of_ks.append(classes.shape[0])
        if cid == 0:
            dist = categories[:, -2].astype(int) - categories[:, -3].astype(int)
            unique = np.unique(classes)
            counts = []
            for ui in unique:
                counts.append(np.sum(dist[classes == ui]))
            counts = np.array(counts)
        else:
            unique, counts = np.unique(classes, return_counts=True)
        unique_dc = dict(zip(unique.astype(int), counts / float(counts.sum())))
        for cls in chrHMMs:
            if cls + 1 in unique_dc.keys():
                class_data[cls][cid] = unique_dc[cls + 1]
            else:
                class_data[cls][cid] = 0.

    plt.rc('font', weight='bold')
    plt.rc('xtick', labelsize=10)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels

    fig, axs = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.4, hspace=0.5, bottom=0.2)
    ax = axs[0]
    pls = []
    for cls in chrHMMs:
        if cls == 0:
            pl = ax.bar(ind, tuple(list(class_data[cls, :])), 0.35)
        else:
            sum_prev = list(np.sum(class_data[0: cls, :], axis=0))
            pl = ax.bar(ind, tuple(list(class_data[cls, :])), 0.35, bottom=tuple(sum_prev))
        for item in pl:
            item.set_color(cm(1. * cls / len(chrmoHMM_labels)))
        pls.append(pl)
    ax.set_ylabel('Percentage', weight='bold', fontsize=16)
    ax.set_xticks(ind)
    ax.set_xticklabels(category_labels)
    ax.set_yticks(np.arange(0.2, 1.1, 0.2))
    ax.set_yticklabels(["%d" % int(item * 100) for item in np.arange(0.2, 1.1, 0.2)])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(tuple(pls), tuple(chrmoHMM_labels), loc='center left', bbox_to_anchor=(1, 0.5))

    column_vec = class_data[:, 0]
    fold_change_data = class_data[:, :] / column_vec[:, None]
    ax = axs[1]
    pls = []
    legend_label = chrmoHMM_labels if not only_enhancer_in_fold_change else ["Strong Enhancer", "Strong Enhancer", "Weak Enhancer", "Weak Enhancer"]
    cls_included = []

    for cls in chrHMMs:
        if not only_enhancer_in_fold_change or (only_enhancer_in_fold_change and "Enhancer" in chrmoHMM_labels[cls]):
            if cls == 0 or not cls_included:
                pl = ax.bar(ind, tuple(list(fold_change_data[cls, :])), 0.35)
            else:
                sum_prev = list(np.sum(fold_change_data[np.array(cls_included), :], axis=0))
                pl = ax.bar(ind, tuple(list(fold_change_data[cls, :])), 0.35, bottom=tuple(sum_prev))
            for item in pl:
                item.set_color(cm(1. * cls / len(chrmoHMM_labels)))
            pls.append(pl)
            cls_included.append(cls)
    cls_sum = np.sum(fold_change_data[np.array(cls_included), :], axis=0)
    for cid, class_label in enumerate(category_labels):
        if cid != 0:
            ax.text(cid - 0.3, cls_sum[cid] + 2, "#CpG:\n%d" % number_of_ks[cid], fontsize=8, color='black', weight='bold')
    ax.set_ylabel('Fold Change', weight='bold', fontsize=16)
    ax.set_xticks(ind)
    ax.set_xticklabels(category_labels)
    ax.set_ylim(0, np.max(cls_sum) + 10 )
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(tuple(pls), tuple(legend_label), loc='center left', bbox_to_anchor=(1, 0.5))

    mkdir(os.path.dirname(fig_path))
    plt.savefig(fig_path, dpi=300)


if __name__ == "__main__":
    args = parser.parse_args()
    chromHMM_fp = args.chromHMM_fp
    k_fp = args.k_fp
    k_in_tfbs_cluster_fp = args.k_in_tfbs_cluster_fp
    k_in_tfbs_jsd_over_95_pct_fp = args.k_in_tfbs_jsd_over_95_pct_fp
    k_category = args.k_category
    min_meth = args.min_meth
    max_meth = args.max_meth
    fig_path = args.fig_path
    plot_chrmoHMM(chromHMM_fp, k_fp, k_category, min_meth, max_meth, fig_path, only_enhancer_in_fold_change=True)
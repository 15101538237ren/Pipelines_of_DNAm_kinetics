# -*- coding:utf-8 -*-
import os, argparse, sys
import pandas as pd
import numpy as np
from scipy.spatial.distance import jensenshannon
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Plot the Histogram iof Ks in TFBSs vs in Genome')

parser.add_argument('-k_fp', type=str, default="../data/K-500bp-CpGd_HMM_DNase_NOS_Rmvd_Reversible_Sites.bed",
                    help='The file path where the bed file of Ks is stored')

parser.add_argument('-tfbs_k_fp', type=str, default="../data/TFBS_clusters/wgEncodeRegTfbsClusteredV3_Ks/POU5F1.bed",
                    help='The file path where K intersected with TFs stored')

parser.add_argument('-min_meth', type=float, default=0,
                    help='The minimal WGBS methylation level for CpG filtering, default: 0')

parser.add_argument('-max_meth', type=float, default=1.01,
                    help='The maximum WGBS methylation level for CpG filtering, default: 1.01')

parser.add_argument('-fig_path', type=str, default="../figures/J_Figures/fig3h_POU5F1.svg", help='The path of figure3a output')

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def plot_k_hist_tfbs_vs_genome(k_fp, tfbs_k_fp, fig_path, min_meth = 0, max_meth = 1.0 + 1e-6, NBIN = 30):
    print("Reading %s" % k_fp)
    df_genome = pd.read_csv(k_fp, sep='\t', header=None).values
    meth_inds = [mid for mid, meth in enumerate(df_genome[:, 3]) if (meth != "." and ((float(meth) > min_meth) and float(meth) < max_meth))]
    ks_genome = df_genome[meth_inds, 4].astype(float)

    df_in_tfbs = pd.read_csv(tfbs_k_fp, sep='\t', header=None).values
    print("Reading %s" % tfbs_k_fp)
    meth_inds = [mid for mid, meth in enumerate(df_in_tfbs[:, 3]) if (meth != "." and ((float(meth) > min_meth) and float(meth) < max_meth))]
    k_in_tfbs = df_in_tfbs[meth_inds, 4].astype(float)

    plt.rc('font', weight='bold')
    plt.rc('xtick', labelsize=12)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=12)  # fontsize of the tick labels
    tfbs_k_dir = os.path.dirname(tfbs_k_fp)
    tfbs_names = []
    for file_name in os.listdir(tfbs_k_dir):
        if file_name.endswith(".bed"):
            tfbs_name = os.path.basename(file_name).split(".")[0]
            tfbs_names.append(tfbs_name)
    tfbs_dict = {}
    for tid, tfbs in enumerate(sorted(tfbs_names)):
        tfbs_dict[tfbs] = tid

    tfbs = os.path.basename(tfbs_k_fp).split(".")[0]
    _, ax = plt.subplots(1, 1, figsize=(5, 5))
    plt.subplots_adjust(left=0.3, bottom=0.2)
    cm = plt.get_cmap('gist_rainbow')
    max_k = 10.0
    bins = [i * (max_k/NBIN) for i in range(NBIN + 1)]
    xhist, _, _ = ax.hist(ks_genome, bins=bins, density=True, color='black', edgecolor='black', alpha=0.5,
                          linewidth=0.5, label="CpGs in the Genome")
    yhist, _, _ = ax.hist(k_in_tfbs, bins=bins, density=True, color=cm(1. * tfbs_dict[tfbs] / 161), edgecolor='black',
                          alpha=0.5, linewidth=0.5, label="CpGs in %s binding sites" % tfbs)
    jsd = jensenshannon(xhist, yhist)
    ax.legend(loc='upper right')
    ax.text(5, .5, "JSD=%.2f" % jsd, color="black", weight='bold', fontsize=12)
    ax.set_xlim([0, max_k])
    ax.set_ylim([0, 0.95])
    ax.set_xlabel('Remethylation Rate($hr^{-1}$)', weight='bold', fontsize=14)
    ax.set_ylabel('Probability', weight='bold', fontsize=14)
    mkdir(os.path.dirname(fig_path))
    plt.savefig(fig_path)


if __name__ == "__main__":
    args = parser.parse_args()
    if len(sys.argv) < 3:
        parser.print_help()
    else:
        k_fp = args.k_fp
        tfbs_k_fp = args.tfbs_k_fp
        min_meth = args.min_meth
        max_meth = args.max_meth
        fig_path = args.fig_path
        plot_k_hist_tfbs_vs_genome(k_fp, tfbs_k_fp, fig_path, min_meth, max_meth)

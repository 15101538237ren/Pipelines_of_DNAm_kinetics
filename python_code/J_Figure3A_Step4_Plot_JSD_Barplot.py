# -*- coding:utf-8 -*-
import os, argparse, sys, math
import pandas as pd
import numpy as np
from scipy.spatial.distance import jensenshannon
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Plot the Barplot of JSDs for the K in TFBSs')

parser.add_argument('-k_fp', type=str, default="../data/K-500bp-CpGd_HMM_DNase_NOS_Rmvd_Reversible_Sites.bed",
                    help='The file path where the bed file of Ks is stored')

parser.add_argument('-tfbs_k_dir', type=str, default="../data/TFBS_clusters/wgEncodeRegTfbsClusteredV3_Ks",
                    help='The directory where K intersected with TFs stored')

parser.add_argument('-min_meth', type=float, default=0,
                    help='The minimal WGBS methylation level for CpG filtering, default: 0')

parser.add_argument('-max_meth', type=float, default=1.01,
                    help='The maximum WGBS methylation level for CpG filtering, default: 1.01')

parser.add_argument('-percentile', type=int, default=95,
                    help='The percentile(int) for filtering high jsd tfbs, default: 95')

parser.add_argument('-with_tfbs', type=int, default=0,
                    help='Whether add tfbs names on x-axis ticks, default: 0(False)')

parser.add_argument('-fig_path', type=str, default="../figures/J_Figures/fig3h.svg", help='The path of figure3a output')

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def calc_jsd_for_tfbs(k_fp, tfbs_k_dir, min_meth = 0, max_meth = 1.0 + 1e-6, NBIN = 30):
    print("Reading %s" % k_fp)
    df_genome = pd.read_csv(k_fp, sep='\t', header=None).values
    meth_inds = [mid for mid, meth in enumerate(df_genome[:, 3]) if (meth != "." and ((float(meth) > min_meth) and float(meth) < max_meth))]
    ks_genome = df_genome[meth_inds, 4].astype(float)
    max_k = 10.0
    bins = [i * (max_k / NBIN) for i in range(NBIN + 1)]
    hist_of_genome_k, _ = np.histogram(ks_genome, bins=bins, density=True)
    dict_of_jsd = {}

    for file_name in os.listdir(tfbs_k_dir):
        if file_name.endswith(".bed"):
            tfbs_name = os.path.basename(file_name).split(".")[0]
            input_fp = os.path.join(tfbs_k_dir, file_name)
            df_in_tfbs = pd.read_csv(input_fp, sep='\t', header=None).values
            print("Reading %s" % file_name)
            meth_inds = [mid for mid, meth in enumerate(df_in_tfbs[:, 3]) if (meth != "." and ((float(meth) > min_meth) and float(meth) < max_meth))]
            k_in_tfbs = df_in_tfbs[meth_inds, 4].astype(float)
            hist_of_k_in_tfbs, _ = np.histogram(k_in_tfbs, bins=bins, density=True)
            jsd = jensenshannon(hist_of_genome_k, hist_of_k_in_tfbs)
            dict_of_jsd[tfbs_name] = jsd

    jsd_tsv_fp = os.path.join(tfbs_k_dir, "jsd.tsv")
    with open(jsd_tsv_fp, "w") as jsd_f:
        for tfbs in sorted(dict_of_jsd.keys()):
            ltw = "%s\t%.3f\n" % (tfbs, dict_of_jsd[tfbs] )
            jsd_f.write(ltw)
    print("Saved %s successful!" % jsd_tsv_fp)

def barplot_for_jsd(tfbs_k_dir, fig_path, percentile= 95, with_tfbs=False):
    jsd_tsv_fp = os.path.join(tfbs_k_dir, "jsd.tsv")
    df_jsd = pd.read_csv(jsd_tsv_fp, sep='\t', header=None).values
    tfbs_names = df_jsd[:, 0].astype(str)
    jsds = df_jsd[:, 1].astype(float)
    jsd_pct = np.percentile(jsds, percentile)
    print("%d%% Percentile of JSD is: %.2f" % (percentile, float(jsd_pct)))
    n_tfbs = len(tfbs_names)

    plt.rc('font', weight='bold')
    plt.rc('xtick', labelsize=6)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)  # fontsize of the tick labels

    _, ax = plt.subplots(1, 1, figsize=(12, 5))
    ind = np.arange(n_tfbs)
    ax.bar(ind, jsds, color="blue", edgecolor='black', linewidth=0.5)
    ax.set_xlim([-0.5, n_tfbs])
    if with_tfbs:
        ax.set_xticks(ind)
        ax.set_xticklabels(tfbs_names, rotation=90)
    else:
        ax.set_xticks([])
    ax.axhline(y=jsd_pct, color='red', linestyle='--')
    ax.text(n_tfbs * 0.85, jsd_pct + 0.05, "%d%% JSD\n Percentile" % percentile, color="red", fontsize=12)
    for tid, tfbs in enumerate(tfbs_names):
        jsd = jsds[tid]
        if jsd > jsd_pct or tfbs in ["CTCF", "GATA1"]:
            ax.text(tid - 3.5, jsd + 0.012, tfbs, color="black", weight='bold', fontsize=12)
            if tfbs in ["CTCF", "GATA1"]:
                ax.scatter(tid, jsd + 0.008, s=10, marker='*', color='black')
    ylim = 0.6
    ax.set_yticks([0.1 * i for i in range(int(math.floor(ylim/0.1)) + 1)])
    ax.set_ylim([0, ylim])
    ax.set_xlabel('Binding regions for specific transcription factors', weight='bold', fontsize=16)
    ax.set_ylabel('Jensen-Shannon Divergence', weight='bold', fontsize=16)
    mkdir(os.path.dirname(fig_path))
    plt.savefig(fig_path)

if __name__ == "__main__":
    args = parser.parse_args()
    if len(sys.argv) < 3:
        parser.print_help()
    else:
        k_fp = args.k_fp
        tfbs_k_dir = args.tfbs_k_dir
        min_meth = args.min_meth
        max_meth = args.max_meth
        fig_path = args.fig_path
        percentile = args.percentile
        with_tfbs = args.with_tfbs
        #calc_jsd_for_tfbs(k_fp, tfbs_k_dir, min_meth, max_meth)
        barplot_for_jsd(tfbs_k_dir, fig_path, percentile, with_tfbs)
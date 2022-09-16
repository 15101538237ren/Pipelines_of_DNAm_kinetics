# -*- coding:utf-8 -*-
import os, argparse, math, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import variation

parser = argparse.ArgumentParser(description='Plot the Variation of Expression vs K near Gene TSS')

parser.add_argument('-k_gene_fp', type=str, default="../data/hg19_genes/GENE_Ks.bed",
                    help='The file path where K distance to Gene TSS stored')

parser.add_argument('-gene_fp', type=str, default="../data/hg19_genes/GENE.bed",
                    help='The file path where K distance to Gene TSS stored')

parser.add_argument('-h1_exprs_matrix_fp', type=str, default="../data/single-cell-data/hESC_expression_matrix.csv",
                    help='The csv file where single cell expression matrix at hESC stage are stored')

parser.add_argument('-max_dist', type=int, default=15000,
                    help='The max distance to TSS for CgP filtering')

parser.add_argument('-min_meth', type=float, default=0,
                    help='The min WGBS methylation level for CpG filtering, default: 0')

parser.add_argument('-max_meth', type=float, default=1.01,
                    help='The max WGBS methylation level for CpG filtering, default: 1.01')

parser.add_argument('-n_category', type=int, default=5,
                    help='The number of categories for Gene Expression Variation')

parser.add_argument('-fig_path', type=str, default="../figures/J_Figures/fig3f.svg", help='The path of figure3f output')

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def plot_k_near_tss_categorized_by_exprs_coef_var(k_gene_fp, h1_exprs_matrix_fp, max_dist, fig_path, n_category = 5, minCpG=20, window_sz = 2000, step_sz = 200, min_samples=3):
    print("Start Execute Plotting the Variation of Expression vs K near Gene TSS")
    k_dist_df = pd.read_csv(k_gene_fp, sep="\t", header=None).values
    gene_df = pd.read_csv(args.gene_fp, sep="\t", header=None).values.astype(str)
    genes = gene_df[:, -1]
    gene_dict = {gene: [gene_df[gid, 0], gene_df[gid, 1], gene_df[gid, 2], gene_df[gid, 3]]  for gid, gene in enumerate(genes) }
    meth_inds = [mid for mid, meth in enumerate(k_dist_df[:, 3]) if
                 (meth != "." and ((float(meth) > min_meth) and float(meth) < max_meth))]
    k_dist_df = k_dist_df[meth_inds, :]
    dists = k_dist_df[:, -1].astype(int)
    dist_ind = np.logical_and(dists >= -max_dist, dists <= max_dist)
    k_dist_df = k_dist_df[dist_ind, :]
    genes, counts = np.unique(k_dist_df[:, -2], return_counts=True)
    filtered_genes = np.array([genes[cid] for cid, count in enumerate(counts) if count > minCpG])
    print("Preprocessed Ks with distance bed file")

    hesc_exprs_matrix = pd.read_csv(h1_exprs_matrix_fp, sep=",", header=0, index_col=0)
    hesc_exprs = hesc_exprs_matrix.values

    genes_in_exprs_mat = hesc_exprs_matrix.index.to_numpy()
    overlapped_genes = np.in1d(genes_in_exprs_mat, filtered_genes).nonzero()[0].tolist()
    valid_genes_indices = np.where(np.sum(hesc_exprs_matrix > 0, axis=1) >= min_samples)[0].tolist()
    intersected_gene_idxs = intersection(overlapped_genes, valid_genes_indices)
    cv_of_intersect_genes = np.array(
        [variation(hesc_exprs[idx, hesc_exprs[idx, :] > 0]) for idx in intersected_gene_idxs])
    intersected_gene_names = genes_in_exprs_mat[intersected_gene_idxs]
    sorted_indices_of_expr_cv_arr = np.argsort(cv_of_intersect_genes)[::-1]
    n_genes = len(intersected_gene_idxs)
    category_names = [str(i + 1) for i in
                      range(n_category)]  # "Highest Exprs Var Genes", "Higher Var", "Mid Var", "Lower Var", "Lowest Var"
    each_proportion = int(math.floor(float(n_genes) / n_category))
    print("Preprocessed Expression Matrix")

    plt.rc('font', weight='bold')
    plt.rc('xtick', labelsize=12)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=12)  # fontsize of the tick labels

    fig, ax = plt.subplots(1, 1, figsize=(5,  3))
    plt.subplots_adjust(hspace=0.5, bottom=0.2)
    NSTEPS = int((2 * max_dist - window_sz) / step_sz)

    cm = plt.get_cmap('gist_rainbow')
    genes_in_category = []
    for i in range(n_category):
        print("Preprocessed Category %d / %d" % (i + 1, n_category))
        gene_names_in_proprotion = []
        for j in range(each_proportion):
            idx = i * each_proportion + j
            gene_idx = sorted_indices_of_expr_cv_arr[idx]
            gene_name = intersected_gene_names[gene_idx]
            gene_names_in_proprotion.append(gene_name)
            variation_level = cv_of_intersect_genes[gene_idx]
            genes_in_category.append(gene_dict[gene_name] + [gene_name, variation_level, i + 1])
        gene_names_in_last_cols = k_dist_df[:, -2].astype(str)
        gene_names_in_proprotion = np.array(gene_names_in_proprotion).astype(str)
        k_gene_prop_inds = np.in1d(gene_names_in_last_cols, gene_names_in_proprotion)

        k_and_dist = k_dist_df[k_gene_prop_inds, :][:, [4, 6]].astype("float")
        d_of_each_bin = []
        avg_of_each_bin = []
        for di in range(NSTEPS):
            d1 = -max_dist + step_sz * di
            d2 = -max_dist + window_sz + step_sz * di
            d_of_each_bin.append((d1 + d2) / 2.)
            avg_of_each_bin.append(
                np.mean(k_and_dist[np.logical_and(k_and_dist[:, 1] > d1, k_and_dist[:, 1] < d2), 0]))
        avg_of_each_bin = np.array(avg_of_each_bin)
        ax.plot(d_of_each_bin, avg_of_each_bin, color=cm(1. * i / n_category), label=category_names[i])

    export_df = pd.DataFrame(genes_in_category, columns=["Chr","Start","End","Strand","Gene", "Variation", "Category"])
    export_df.to_csv("../data/Genes_with_scExprVariation_Categories.csv", index=False)
    ax.set_xticks([-max_dist, 0, max_dist])
    ax.set_xlim([-max_dist, max_dist])
    ax.set_yticks([3, 3.5])
    ax.set_ylim([2.95, 3.5])
    ax.set_xlabel("Distance to gene body region (bp)", weight='bold', fontsize=12)
    ax.set_ylabel('Remethylation Rate($hr^{-1}$)', weight='bold', fontsize=12)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(fig_path, dpi=300)
    print("Plot %s successful!" % fig_path)

if __name__ == "__main__":
    args = parser.parse_args()
    if len(sys.argv) < 3:
        parser.print_help()
    else:
        k_gene_fp = args.k_gene_fp
        h1_exprs_matrix_fp = args.h1_exprs_matrix_fp
        n_category = args.n_category
        max_dist = args.max_dist
        min_meth = args.min_meth
        max_meth = args.max_meth
        fig_path = args.fig_path
        plot_k_near_tss_categorized_by_exprs_coef_var(k_gene_fp, h1_exprs_matrix_fp, max_dist, fig_path, n_category=n_category)

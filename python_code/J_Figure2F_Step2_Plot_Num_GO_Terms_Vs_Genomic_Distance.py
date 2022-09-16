# -*- coding:utf-8 -*-
import os, argparse, sys, collections
import pandas as pd
from pandas.errors import ParserError
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import interpolate

parser = argparse.ArgumentParser(description='Step 2 for fig2f: Plot # of GO Terms Vs Genomic Distance')

parser.add_argument('-great_results_dir', type=str, default="../data/GREAT_results/IMSK", help='The directory where GREAT exported results(.tsv) for subcategorized Ks stored')

parser.add_argument('-fig_path', type=str, default="../figures/J_Figures/fig2f.svg", help='The path of figure2f output')


def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def plot_distance_effect_on_go_term(great_results_dir, fig_path, min_dist = -2, max_dist = 28):
    GO_CLASS_NAMES = ["GO Biological Process"]
    COL_NAMES = {'CATEGORY': '# Ontology', 'GO_DESCRIPTION': ' Term Name ', 'GO_ID': ' Term ID ',
                 'Q_Value': '  Hyper FDR Q-Val'}
    base_distance = 5
    plt.rc('font', weight='bold')
    plt.rc('xtick', labelsize=14)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    distances = [item + base_distance for item in [i for i in range(min_dist, max_dist)]]
    dist_and_n = np.zeros((len(distances), 2))
    dist_and_n_dev_go_terms = np.zeros((len(distances), 2))
    dist_and_number_of_go_term_dict = {}
    dist_and_number_of_developmental_go_term_dict = {}
    for file_name in os.listdir(great_results_dir):
        if file_name.endswith(".tsv"):
            input_fp = os.path.join(great_results_dir, file_name)
            dist = file_name.strip(".tsv").split("_")[-1]
            try:
                df = pd.read_csv(input_fp, sep='\t', header=1, index_col=False)
                df_categories = list(df[COL_NAMES['CATEGORY']].values)
                go_term_indexs = [idx for idx, item in enumerate(df_categories) if item in GO_CLASS_NAMES]
                num_of_go_term = len(go_term_indexs)
                df_rows = df.iloc[np.array(go_term_indexs), :]
                df_biological_process = list(df_rows.loc[:, COL_NAMES['GO_DESCRIPTION']].values)
                developmetental_go_terms = []
                for item in df_biological_process:
                    for key_word in ["development", "differentiation", "morphogenesis", "pattern specification"]:
                        if key_word in item:
                            developmetental_go_terms.append(item)
                num_of_developmental_go_term = len(developmetental_go_terms)
            except ParserError:
                num_of_go_term = 0
                num_of_developmental_go_term = 0
            dist_and_number_of_go_term_dict[base_distance + int(dist)] = num_of_go_term
            dist_and_number_of_developmental_go_term_dict[base_distance + int(dist)] = num_of_developmental_go_term
    sorted_dict = collections.OrderedDict(sorted(dist_and_number_of_go_term_dict.items()))
    sorted_dev_dict = collections.OrderedDict(sorted(dist_and_number_of_developmental_go_term_dict.items()))
    for kid, dist in enumerate(distances):
        dist_and_n[kid, 0] = dist
        dist_and_n_dev_go_terms[kid, 0] = dist
        if dist in sorted_dict.keys():
            dist_and_n[kid, 1] = sorted_dict[dist]
            dist_and_n_dev_go_terms[kid, 1] = sorted_dev_dict[dist]
        else:
            try:
                dist_and_n[kid, 1] = int(round((sorted_dict[dist -1] + sorted_dict[dist + 1] )/ 2.))
                dist_and_n_dev_go_terms[kid, 1] = int(round((sorted_dev_dict[dist - 1] + sorted_dev_dict[dist + 1]) / 2.))
            except KeyError as e:
                print(e)

    x1 = dist_and_n[:, 0]
    y1 = dist_and_n[:, 1]
    df = pd.DataFrame(dist_and_n, columns=['Distance', 'N_Term'])

    df_fp = os.path.join(great_results_dir, "%s_dist_vs_n_go_term.csv" % great_results_dir.split("/")[-1])
    df.to_csv(df_fp, sep=",", header=True, index=False, )
    tck1 = interpolate.splrep(x1, y1, s=0)
    xnew = np.linspace(3, 39, num=100, endpoint=True)
    ynew = interpolate.splev(xnew, tck1, der=0)
    ax.plot(xnew, ynew, "-", color="red", linewidth=3, label="All Biological Processes")
    ax.plot(x1, y1, '.', markersize=4, color="red")

    x2 = dist_and_n_dev_go_terms[:, 0]
    y2 = dist_and_n_dev_go_terms[:, 1]
    tck2 = interpolate.splrep(x2, y2, s=0)
    ynew2 = interpolate.splev(xnew, tck2, der=0)
    ax.plot(xnew, ynew2, "-", color="blue", linewidth=3, label="Developmental-associated")
    ax.plot(x2, y2, '.', markersize=4, color="blue")
    ax.legend()
    ax.set_xlabel('Max distance to TSS in both directions (Kb)', weight='bold', fontsize=14)
    ax.set_ylabel("# of GO Terms", weight='bold', fontsize=14)
    ax.set_xlim([0, 25])
    ax.set_ylim([-5, 100])
    plt.savefig(fig_path, dpi=300)

if __name__ == "__main__":
    args = parser.parse_args()
    great_results_dir = args.great_results_dir
    fig_path = args.fig_path
    plot_distance_effect_on_go_term(great_results_dir, fig_path, min_dist=-2, max_dist=24)
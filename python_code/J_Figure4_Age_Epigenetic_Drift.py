# -*- coding:utf-8 -*-
import os, argparse, sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Convert reversible sites from mat file to bed file')
parser.add_argument('-region_intersect_dir', type=str, default="../data/CD4Tell_Age_WGBS/Region_Intersection", help='The directory where regional intersected K with Young and Old WGBS values')
parser.add_argument('-fig_path', type=str, default="../figures/J_Figures/fig4.pdf", help='The path of figure output')

FILE_ORDERED_NAMES = {
            "Genomic_Regions":
                ["Promoter", "Enhancer", "CGI", "LINE", "SINE", "LTR"],
            "Histone_Modification":
                ["H3k4me1Not3_no_k27me3_no_k27ac", "H3k4me1Not3_k27me3", "H3k4me1Not3_k27ac", "H3k4me3Not1_no_k27me3_no_k27ac", "H3k4me3Not1_k27me3", "H3k4me3Not1_k27ac", "H3K9Me3", "H3K36me3"],
            "TFBS":
                ["DNase", "CTCF", "EP300", "EZH2", "tfbs_cluster_v3", "general_jsd_over_threshold"]
            }

FILE_LABELS = {
            "Genomic_Regions":
                ["Promoter", "Enhancer", "CGI",
                 "LINE", "SINE", "LTR"],
            "Histone_Modification":
                ["k4me1", "k4me1 + k27me3", "k4me1 + k27ac",
                 "k4me3", "k4me3 + k27me3", "k4me3 + k27ac",
                 "K9Me3", "K36me3"],
            "TFBS":
                ["DNase", "CTCF", "EP300", "EZH2",
                 "TFBS", "JSD > 95%\npercentile TFBS"]
            }

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def plot_region_specific_methy_diff_between_young_and_aged(region_intersect_dir, fig_path):
    regions = ["Genomic_Regions", "Histone_Modification", "TFBS"]
    hypos = ["Hypo", "Hyper"]
    N_ROW, N_COL = 3, 2
    cm = plt.get_cmap('gist_rainbow')
    fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * 8, N_ROW * 4))
    sns.set(color_codes=True)
    # fig.tight_layout()
    # plt.subplots_adjust(left=0.125, right=1.2, bottom=0.3, top=0.9, wspace=0.3, hspace=0.2)
    for hid, hypo in enumerate(hypos):
        for rgid, region in enumerate(regions):
            row = rgid
            col = hid
            ax = axs[col] if N_ROW == 1 else axs[row][col]
            for rlid, region_label in enumerate(FILE_ORDERED_NAMES[region]):
                print("Reading %s" % region_label)
                file_fp = os.path.join(region_intersect_dir, region, "%s.bed" % region_label)
                data_df = pd.read_csv(file_fp, sep="\t", header=None).values
                methy_diff = data_df[:, -1].astype("float")
                if hid == 0:
                    final_methy_diff = methy_diff[methy_diff < 0]
                    ks = data_df[methy_diff < 0, 4].astype("float")
                else:
                    final_methy_diff = methy_diff[methy_diff >= 0]
                    ks = data_df[methy_diff >= 0, 4].astype("float")

                sns.regplot(x=ks, y=final_methy_diff, ax=ax, lowess=True, scatter=False, line_kws={'color': cm(1. * rlid / len(FILE_ORDERED_NAMES[region]))}, label=FILE_LABELS[region][rlid])#n_boot=500, + ": #" + str(ks[ks > 5].shape[0])
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            #ax.set_xlim([0, 10])
            if hid:
                #ax.set_ylim([0, 0.3])
                ax.set_title("Hyper-methy in Aged")
            else:
                ax.set_title("Hypo-methy in Aged")
                #ax.set_ylim([-0.3, 0])
                ax.set_ylim(ax.get_ylim()[::-1])
            ax.set_ylabel("Methylation_{Young}")
    mkdir(os.path.dirname(fig_path)) # Create directory for figure output
    plt.savefig(fig_path, dpi=300)
if __name__ == "__main__":
    args = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
    else:
        region_intersect_dir = args.region_intersect_dir
        fig_path = args.fig_path
        plot_region_specific_methy_diff_between_young_and_aged(region_intersect_dir, fig_path)
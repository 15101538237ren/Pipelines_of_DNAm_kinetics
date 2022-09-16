# -*- coding:utf-8 -*-
import os, argparse, sys
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='1. Extracting H1-ESC Specific TFBS and 2. Split it by TFs')
parser.add_argument('-tfbs_cluster_fp', type=str, default="../data/TFBS_clusters/wgEncodeRegTfbsClusteredV3.bed", help='The file path where wgEncodeRegTfbsClusteredV3.bed stored')

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def get_h1_cell_ids(tfbs_cluster_fp):
    tfbs_cluster_dir = os.path.dirname(tfbs_cluster_fp)
    tfbs_cells_fp = os.path.join(tfbs_cluster_dir, "wgEncodeRegTfbsCellsV3.tab")
    df = pd.read_csv(tfbs_cells_fp, sep='\t', header=None).values
    abrvs = df[:, 1].astype(str)
    h1_ids = df[abrvs == "1", 0].astype(str)
    return h1_ids

def extract_h1_specific_tfbs(tfbs_cluster_fp):
    print("Start Extracting H1-ESC Specific TFBS")
    h1_ids = set(get_h1_cell_ids(tfbs_cluster_fp))
    df = pd.read_csv(tfbs_cluster_fp, sep='\t', header=None).values
    cell_ids = df[:, -2].astype(str)
    h1_row_ids = []
    tfbs_cluster_dir = os.path.dirname(tfbs_cluster_fp)
    h1_specific_tfbs_fp = os.path.join(tfbs_cluster_dir, "h1_specific_tfbs.bed")
    for row_id, cell_id_row in enumerate(cell_ids):
        if row_id and (row_id + 1) % 500000 == 0:
            print("%.1f%% for Extracting H1-ESC Specific TFBS" % ((row_id + 1.0) * 100.0 / cell_ids.shape[0]))
        cids = cell_id_row.strip("\n").split(",")
        for cid in cids:
            if cid in h1_ids:
                h1_row_ids.append(row_id)
                break
    h1_df = df[h1_row_ids, 0:4]
    np.savetxt(h1_specific_tfbs_fp, h1_df, fmt="%s\t%d\t%d\t%s\n", delimiter="", newline="")
    print("Saved %s Successful!" % h1_specific_tfbs_fp)

def split_bed_file_by_tfbs_name(tfbs_fp, tf_col = 3):
    tfbs_cluster_dir = os.path.dirname(tfbs_fp)
    tfbs_fn = os.path.basename(tfbs_fp).split(".")[0]
    tfbs_dir = os.path.join(tfbs_cluster_dir, tfbs_fn)
    mkdir(tfbs_dir)
    try:
        df = pd.read_csv(tfbs_fp, sep='\t', header=None).values
        tfbs_vals = df[:, tf_col]
        tfbs_names = np.unique(tfbs_vals)
        for tid, tfbs in enumerate(tfbs_names):
            print("%.1f%% for Split %s by TF: %s" % ((tid + 1.0) * 100.0 / len(tfbs_names), tfbs_fn, tfbs))
            df_lines = df[tfbs_vals == tfbs, 0: (tf_col + 1)]
            out_fp = os.path.join(tfbs_dir, "%s.bed" % tfbs)
            np.savetxt(out_fp, df_lines[:], fmt="%s\t%d\t%d\t%s\n", delimiter="", newline="")
    except Exception as e:
        print(e)

if __name__ == "__main__":
    args = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
    else:
        tfbs_cluster_fp = args.tfbs_cluster_fp
        extract_h1_specific_tfbs(tfbs_cluster_fp)
        tfbs_fp = os.path.join(os.path.dirname(tfbs_cluster_fp), "h1_specific_tfbs.bed")
        split_bed_file_by_tfbs_name(tfbs_fp)
        tfbs_fp = tfbs_cluster_fp
        split_bed_file_by_tfbs_name(tfbs_fp)
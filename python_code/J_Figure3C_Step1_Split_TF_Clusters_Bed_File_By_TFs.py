# -*- coding:utf-8 -*-
import os, argparse, sys
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Split TF Cluster File By TFs')

parser.add_argument('-tfbs_cluster_fp', type=str, default="../data/TFBS_clusters/wgEncodeRegTfbsClusteredV3.bed", help='The file path where wgEncodeRegTfbsClusteredV3.bed stored')

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def split_bed_file_by_tfbs_name(tfbs_cluster_fp, tf_col = 3):
    tfbs_cluster_dir = os.path.dirname(tfbs_cluster_fp)
    tfbs_fn = os.path.basename(tfbs_cluster_fp).split(".")[0]
    tfbs_dir = os.path.join(tfbs_cluster_dir, "%s_TFs" % tfbs_fn)
    mkdir(tfbs_dir)
    try:
        df = pd.read_csv(tfbs_cluster_fp, sep='\t', header=None).values
        tfbs_vals = df[:, tf_col]
        tfbs_names = np.unique(tfbs_vals)
        tf_lines = { tfbs : [] for tfbs in tfbs_names}
        with open(tfbs_cluster_fp, "r") as tfbs_f:
            lines = tfbs_f.readlines()
            n_line = len(lines)
            for lid, line in enumerate(lines):
                if lid and (lid + 1) % 500000 == 0:
                    print("%.1f%% for Split TF Cluster By TFs" % (
                            (lid + 1.0) * 100.0 / n_line))
                tf = line.split("\t")[3]
                tf_lines[tf].append(line)
        for tid, tfbs in enumerate(tfbs_names):
            out_fp = os.path.join(tfbs_dir, "%s.bed" % tfbs)
            with open(out_fp, "w") as out_f:
                print("Writing for %s" % tfbs)
                for line in tf_lines[tfbs]:
                    out_f.write(line)
    except Exception as e:
        print(e)

if __name__ == "__main__":
    args = parser.parse_args()
    tfbs_cluster_fp = args.tfbs_cluster_fp
    split_bed_file_by_tfbs_name(tfbs_cluster_fp)
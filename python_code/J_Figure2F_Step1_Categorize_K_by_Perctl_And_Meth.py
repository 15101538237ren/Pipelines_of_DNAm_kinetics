# -*- coding:utf-8 -*-
import os, argparse, sys, enum, random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Step 1 for fig2f: Categorize K by Percentile and WGBS Methylation Level')

parser.add_argument('-k_fp', type=str, default="../data/K-500bp-CpGd_HMM_DNase_NOS_Rmvd_Reversible_Sites.bed", help= 'The file path where the bed file of Ks is stored')

parser.add_argument('-min_meth', type=float, default=0.5,
                    help='The minimal WGBS methylation level for CpG filtering, default: 0')

parser.add_argument('-max_meth', type=float, default=1.01,
                    help='The maximum WGBS methylation level for CpG filtering, default: 1.01')

parser.add_argument('-m_category', type=str, default="IM", choices=["GM", "IM", "HM"],
                    help='The category of methylation level for CpGs filtering, choices are: IM (Intermediate-methylated), HM (highly-methylated), GM (greater-than-half-methylated)')

parser.add_argument('-k_category', type=str, default="slow", choices=["slow", "medium", "fast", "all"],
                    help='The category of K values')

# Enum for size units
class SIZE_UNIT(enum.Enum):
   BYTES = 1
   KB = 2
   MB = 3
   GB = 4
def convert_unit(size_in_bytes, unit):
   """ Convert the size from bytes to other units like KB, MB or GB"""
   if unit == SIZE_UNIT.KB:
       return size_in_bytes/1024
   elif unit == SIZE_UNIT.MB:
       return size_in_bytes/(1024*1024)
   elif unit == SIZE_UNIT.GB:
       return size_in_bytes/(1024*1024*1024)
   else:
       return size_in_bytes

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def get_meth_indexs(df, m_category):
    if m_category == "GM":
        max_meth = 1.01
        min_meth = 0.5
    elif m_category == "IM":
        max_meth = 0.8
        min_meth = 0.5
    else:
        max_meth = 1.01
        min_meth = 0.8
    meth_inds = [mid for mid, meth in enumerate(df[:, 3]) if
                 (meth != "." and ((float(meth) > min_meth) and float(meth) < max_meth))]
    return df[meth_inds, :], min_meth, max_meth

def categorize_k_and_export(k_fp, m_category, k_category):
    k_df = pd.read_csv(k_fp, sep="\t", header=None).values
    k_df, min_meth, max_meth = get_meth_indexs(k_df, m_category)
    ks = k_df[:, 4].astype(float)
    slow_k_pct = np.percentile(ks, 33)
    mid_k_pct = np.percentile(ks, 67)
    k_label = "AK"
    if k_category == "slow":
        k_df = k_df[ks < slow_k_pct, :]
        k_label = "SK"
    elif k_category == "medium":
        k_df = k_df[np.logical_and(ks > slow_k_pct, ks < mid_k_pct), :]
        k_label = "MK"
    elif k_category == "fast":
        k_df = k_df[ks > mid_k_pct, :]
        k_label = "FK"
    else:
        pass
    categoried_k_dir = os.path.join(os.path.dirname(k_fp), "Categorized_k")
    mkdir(categoried_k_dir)
    export_fp = os.path.join(categoried_k_dir, "%s%s.bed" % (m_category, k_label))
    np.savetxt(export_fp, k_df[:, :3], fmt="%s\t%d\t%d\n", delimiter="", newline="")
    print("Exported k with %s kinetics and %.2f < WGBS methylation < %.2f ! at %s" % (k_category.capitalize(), min_meth, max_meth, export_fp))
    return export_fp

def down_sampling_to_target_size(input_fp, target_file_sz_in_MB = 4.7, seed = 41):
    byte_size = os.path.getsize(input_fp)
    mb_size = convert_unit(byte_size, SIZE_UNIT.MB)
    sampling_ratio = target_file_sz_in_MB/mb_size
    k_df = pd.read_csv(input_fp, sep="\t", header=None).values
    n_row = k_df.shape[0]
    target_n_row = int(n_row * sampling_ratio)
    row_ids = range(n_row)
    random.seed(seed)
    sampled_row_ids = random.sample(row_ids, target_n_row)
    k_df_sampled = k_df[sampled_row_ids, :]
    file_name = os.path.basename(input_fp).split(".")[0]
    export_fp = os.path.join(os.path.dirname(input_fp), "%s_sampled.bed" % file_name)
    np.savetxt(export_fp, k_df_sampled[:, :], fmt="%s\t%d\t%d\n", delimiter="", newline="")
    print("Exported %s with sampling ratio %.2f!" % (export_fp, sampling_ratio))
    return export_fp

if __name__ == "__main__":
    args = parser.parse_args()
    k_fp = args.k_fp
    k_category = args.k_category
    m_category = args.m_category
    categorized_fp = categorize_k_and_export(k_fp, m_category, k_category)
    if m_category != "IM":
        down_sampling_to_target_size(categorized_fp, target_file_sz_in_MB=4.7)
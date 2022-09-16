# -*- coding:utf-8 -*-
import os, argparse, sys, re
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Preprocess Genomic Annotations')
parser.add_argument('-annotation_dir', type=str, default="../data/Genomic_Context/RAW", help='The data directory where annotations bedfiles are stored')
parser.add_argument('-out_dir', type=str, default="../data/Genomic_Context/Genomic_Regions", help='The output directory')

chromosomes = ["chr%d" % i for i in range(1, 23)]
def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def process_non_repeats(region, args):
    in_fp = os.path.join(args.annotation_dir, "%s.bed" % region)
    out_fp = os.path.join(args.out_dir, "%s.bed" % region)
    n_skip = 1 if region == "CGI" else None
    data_df = pd.read_csv(in_fp, sep="\t", header=None, skiprows=n_skip).values
    data_df = data_df[:, 0:3].astype(str)
    data_df = data_df[np.isin(data_df[:, 0], chromosomes), :]
    np.savetxt(out_fp, data_df, delimiter="\t", fmt='%s\t%s\t%s')
    print("Saved %s sucessful!" % out_fp)

def process_Repeats(region, args):
    in_fp = os.path.join(args.annotation_dir, "%s.bed" % region)
    data_df = pd.read_csv(in_fp, sep="\t", header=None, skiprows=1).values
    data_df = data_df[:, [0, 1, 2, 4]].astype(str)
    data_df = data_df[np.isin(data_df[:, 0], chromosomes), :]

    repeat_classes = ["LINE", "SINE", "LTR"]
    for repeat_class in repeat_classes:
        out_fp = os.path.join(args.out_dir, "%s.bed" % repeat_class)
        repeat_df = data_df[data_df[:, -1] == repeat_class, 0:3]
        np.savetxt(out_fp, repeat_df, delimiter="\t", fmt='%s\t%s\t%s')
        print("Saved %s sucessful!" % out_fp)

def preprocess_region(region, args):
    if region != "Repeats":
        process_non_repeats(region, args)
    else:
        process_Repeats(region, args)

if __name__ == "__main__":
    args = parser.parse_args()
    mkdir(args.out_dir)

    regions = ["Repeats"]#"CGI", "Exons", "Introns", "5UTR", "3UTR", "Repeats"
    for region in regions:
        preprocess_region(region, args)
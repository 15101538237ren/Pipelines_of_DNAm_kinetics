# -*- coding:utf-8 -*-
import os, argparse
import scipy.io as sio
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Convert from bed file to mat file or vice versa')
parser.add_argument('-bed_dir', type=str, default="../data/K_region_intersect/Genomic_Regions", help='The bed file directory')
parser.add_argument('-mat_dir', type=str, default="../data/K_region_intersect/Genomic_Regions_mat", help='The mat file directory')

chr_offset_sz = 1000000000
chromosomes = ["chr%d" % j for j in range(1, 23)] + ["chrX", "chrY"]
chrom_start_loci = {chrm : (i + 1) * chr_offset_sz for i, chrm in enumerate(chromosomes)}

Genomic_Regions = ["Genome"]#"CGI", "Enhancer", "SINE", "DNase", "Promoter",
                                # "TFBR", "Exons", "Introns", "Intergenic", "5UTR", "3UTR",
                                # "H2AZ", "LINE", "LTR"]

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def convert_from_bed_to_mat(bed_fp, with_region_index=False):
    bed_df = pd.read_csv(bed_fp, sep="\t", header=None).values
    for ir, row in enumerate(bed_df):
        chrom = str(row[0])
        bed_df[ir, 1 : 3] += chrom_start_loci[chrom]

    off_set = 1 if with_region_index else 0
    bed_df = bed_df[bed_df[:, 8 + off_set] != ".", :]
    bed_df = bed_df[bed_df[:, 6 + off_set] != ".", :]
    bed_df[bed_df[:, 3] == ".", 3] = -1

    fittedSites = bed_df[:, 1].astype(int)
    inferedRates = bed_df[:, 4].astype(float)
    inferredMethyFrac = bed_df[:, 5].astype(float)

    WGBS = bed_df[:, 3].astype(float)
    CpGd = bed_df[:, 6 + off_set].astype(int)
    DNase = bed_df[:, 8 + off_set].astype(float)
    base_name = os.path.basename(bed_fp).split(".")[0]
    suffix = "_with_Region_Index" if with_region_index else ""
    mat_fp = os.path.join(args.mat_dir, "kineticRates_%s%s.mat" % (base_name, suffix))
    data_dict = {'fittedSites': np.transpose(fittedSites), 'inferedRates': inferedRates, 'inferredMethyFrac': inferredMethyFrac, 'WGBS': WGBS, "CpGd": CpGd, "DNase": DNase}
    if with_region_index:
        RegionIndex = bed_df[:, 6].astype(int)
        data_dict["RegionIndex"] = RegionIndex
    sio.savemat(mat_fp, data_dict)
    print("Saved at %s successful!" % mat_fp)

if __name__ == "__main__":
    args = parser.parse_args()
    bed_dir = args.bed_dir
    mkdir(args.mat_dir)
    for region in Genomic_Regions:
        bed_fp = os.path.join(bed_dir, "%s.bed" % region)
        convert_from_bed_to_mat(bed_fp, with_region_index=False)
    # convert_from_bed_to_mat(bed_fp)
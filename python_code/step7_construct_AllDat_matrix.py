# -*- coding:utf-8 -*-
import os, argparse, sys, re
import numpy as np
import pandas as pd
import scipy.io as sio

parser = argparse.ArgumentParser(description='Step 7: Construct AllDat Matrix in .mat format for each chromosome')
parser.add_argument('-data_dir', type=str, default="../data/", help='The data directory where merged 0h, 1h, 4h, 16h bed files for each chromosomes stored')
parser.add_argument('-out_dir', type=str, default="../data/AllDat", help='The output directory to store the AllDat_chrx.mat files')

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def construct_allDat_matrix(data_dir, out_dir, sep="\t"):
    TIME_POINTS = ['0h', '1h', '4h', '16h']
    CHROMOSOMES = [str(i) for i in range(1, 23)]
    for chromosome in CHROMOSOMES:
        data_dict = {} # the dict to temperatorily store reads information
        for t_idx, time_point in enumerate(TIME_POINTS):
            repli_name = time_point + '_' + 'rep0'
            input_fp = os.path.join(data_dir, repli_name, 'chr' + chromosome + '.bed')
            with open(input_fp, "r") as input_file:
                print("Processing Chr %s, Time %s" % (chromosome, time_point))
                lines = input_file.readlines()
                for lid, line in enumerate(lines):
                    if lid and (lid + 1) % 100000 == 0:
                        print("%.1f%% for Chr %s, Time %s" % ((lid + 1.0) * 100.0 / len(lines), chromosome, time_point))
                    line_contents = re.split(sep, line.strip("\n"))
                    try:
                        chr_i, start, end, methy_reads, unmethy_reads = line_contents
                        methy_reads = int(methy_reads)
                        unmethy_reads = int(unmethy_reads)
                        site = int(start)
                        if site not in data_dict.keys():
                            data_dict[site] = np.zeros((len(TIME_POINTS), 2))
                        data_dict[site][t_idx][0] = methy_reads
                        data_dict[site][t_idx][1] = unmethy_reads
                    except ValueError as e:
                        print("ValueError in Chr %s, Time %s, Line %d :%s" % (chromosome, time_point, lid + 1, line))
        sorted_site_data = sorted(data_dict.items())
        sites = []
        mat_data = []
        for site, reads in sorted_site_data:
            sites.append(site)
            mat_data.append(reads)
        mat_data = np.array(mat_data)
        MAT_DICT = {'AllDat': mat_data, 'sites': sites}
        mkdir(out_dir)
        sio.savemat(os.path.join(out_dir, "AllDat_chr%s.mat" % chromosome), MAT_DICT)
        print("**!!FINISHED Chr %s PROCESSING!!**" % chromosome)

def construct_allDat_bed_file(data_dir, out_dir, sep="\t"):

    TIME_POINTS = ['0h', '1h', '4h', '16h']

    CHROMOSOMES = [str(i) for i in range(1, 23)]

    ltws = []

    for chromosome in CHROMOSOMES:
        data_dict = {} # the dict to temperatorily store reads information
        for t_idx, time_point in enumerate(TIME_POINTS):
            repli_name = time_point + '_' + 'rep0'
            input_fp = os.path.join(data_dir, repli_name, 'chr' + chromosome + '.bed')
            with open(input_fp, "r") as input_file:
                print("Processing Chr %s, Time %s" % (chromosome, time_point))
                lines = input_file.readlines()
                for lid, line in enumerate(lines):
                    if lid and (lid + 1) % 100000 == 0:
                        print("%.1f%% for Chr %s, Time %s" % ((lid + 1.0) * 100.0 / len(lines), chromosome, time_point))
                    line_contents = re.split(sep, line.strip("\n"))
                    try:
                        chr_i, start, end, methy_reads, unmethy_reads = line_contents
                        methy_reads = int(methy_reads)
                        unmethy_reads = int(unmethy_reads)
                        site = int(start)
                        if site not in data_dict.keys():
                            data_dict[site] = np.zeros((len(TIME_POINTS), 2))
                        data_dict[site][t_idx][0] = methy_reads
                        data_dict[site][t_idx][1] = unmethy_reads
                    except ValueError as e:
                        print("ValueError in Chr %s, Time %s, Line %d :%s" % (chromosome, time_point, lid + 1, line))

        for site, dct in data_dict.items():
            vals = ["chr%s" % chromosome, str(site), str(site + 1)]
            for meth in range(2):
                for t_idx, time_point in enumerate(TIME_POINTS):
                    vals.append(str(dct[t_idx][meth]))
            ltw = "\t".join(vals) + "\n"
            ltws.append(ltw)

    mkdir(out_dir)
    outfp = os.path.join(out_dir, "AllDat.bed")
    with open(outfp, "w") as alldata_f:
        alldata_f.writelines(ltws)
    print("**!!FINISHED WRITING into %s !!**" % outfp)

def get_chrom_start_locis():
    chr_offset_sz = 1000000000
    chromosomes = ["chr%d" % j for j in range(1, 23)] + ["chrX", "chrY"]
    chrom_start_loci = {chrm: (i + 1) * chr_offset_sz for i, chrm in enumerate(chromosomes)}
    return chrom_start_loci

def convert_allDat_bed_to_mat(data_fp, out_fp):
    print("**!!START %s PROCESSING!!**" % data_fp)
    data_df = pd.read_csv(data_fp, sep="\t", header=None).values
    TIME_POINTS = ['0h', '1h', '4h', '16h']
    chrom_start_loci_dict = get_chrom_start_locis()
    data_dict = {}# the dict to temperatorily store reads information
    for item in data_df:
        chrom, start, end = item[0], item[1], item[2]
        site = chrom_start_loci_dict[chrom] + int(start)
        offset = 0
        if site not in data_dict.keys():
            data_dict[site] = np.zeros((len(TIME_POINTS), 2))
        for meth in range(2):
            for t_idx, time_point in enumerate(TIME_POINTS):
                data_dict[site][t_idx][meth] = int(item[3+offset])
                offset += 1
    sorted_site_data = sorted(data_dict.items())
    sites = []
    mat_data = []
    for site, reads in sorted_site_data:
        sites.append(site)
        mat_data.append(reads)
    mat_data = np.array(mat_data)
    MAT_DICT = {'AllDat': mat_data, 'sites': sites}
    sio.savemat(out_fp, MAT_DICT)
    print("**!!FINISHED %s PROCESSING!!**" % data_fp)
if __name__ == "__main__":
    args = parser.parse_args()
    # if len(sys.argv) <= 1:
    #     parser.print_help()
    # else:
    # data_dir = args.data_dir
    # out_dir = args.out_dir
    # construct_allDat_bed_file(data_dir, out_dir)
    regions = ["CGI", "Promoter",
                                "TFBR", "Exons", "Introns", "Intergenic", "5UTR", "3UTR",
                                "H2AZ", "LINE", "SINE", "LTR", "DNase", "Enhancer"]
    input_dir = "../data/AllDat/AllDat_Region_Intersect"
    output_dir = "../data/AllDat/AllDat_Region_Intersect_mat"
    mkdir(output_dir)
    for region in regions:
        input_fp = os.path.join(input_dir, "AllDat_%s.bed" % (region))
        output_fp = os.path.join(output_dir, "AllDat_%s.mat" % (region))
        convert_allDat_bed_to_mat(input_fp, output_fp)
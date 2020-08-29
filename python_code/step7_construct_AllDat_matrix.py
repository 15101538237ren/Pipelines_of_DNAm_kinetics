# -*- coding:utf-8 -*-
import os, argparse, sys, re
import numpy as np
import scipy.io as sio

parser = argparse.ArgumentParser(description='Step 7: Construct AllDat Matrix in .mat format for each chromosome')
parser.add_argument('-data_dir', type=str, help='The data directory where merged 0h, 1h, 4h, 16h bed files for each chromosomes stored')
parser.add_argument('-out_dir', type=str, help='The output directory to store the AllDat_chrx.mat files')

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

if __name__ == "__main__":
    args = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
    else:
        data_dir = args.data_dir
        out_dir = args.out_dir
        construct_allDat_matrix(data_dir, out_dir)
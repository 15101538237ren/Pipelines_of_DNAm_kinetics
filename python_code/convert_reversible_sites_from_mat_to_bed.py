# -*- coding:utf-8 -*-
import os, argparse, sys
import scipy.io as sio

parser = argparse.ArgumentParser(description='Convert reversible sites from mat file to bed file')
parser.add_argument('-sites_dir', type=str, default="../data/ReversibleSites", help='The data directory where reversible sites stored')

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def convert_reversible_sites_from_mat_to_bed(sites_dir):
    CHROMOSOMES = [str(i) for i in range(1, 23)]
    ltws = []
    for chromosome in CHROMOSOMES:
        print("Started chr%s in %s" % (chromosome, sites_dir))
        data_fp = os.path.join(sites_dir, "ReversibleSites_chr%s.mat" % chromosome)
        mat_data = sio.loadmat(data_fp)
        reversibleSites = mat_data["reversibleSites"][:, 0]
        for sid, site in enumerate(reversibleSites):
            ltw = "chr%s\t%d\t%d\n" % (chromosome, site, site + 1)
            ltws.append(ltw)
    k_file_name = "K_reversible.bed"
    output_fp = os.path.join(sites_dir, k_file_name)
    with open(output_fp, "w") as out_f:
        ltw = "".join(ltws)
        out_f.write(ltw)
        print("Write %s Successful!" % k_file_name)

if __name__ == "__main__":
    args = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
    else:
        sites_dir = args.sites_dir
        convert_reversible_sites_from_mat_to_bed(sites_dir)
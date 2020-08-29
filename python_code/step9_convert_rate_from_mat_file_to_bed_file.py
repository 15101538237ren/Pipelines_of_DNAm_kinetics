# -*- coding:utf-8 -*-
import os, argparse, sys
import scipy.io as sio

parser = argparse.ArgumentParser(description='Step 9: Convert rate from mat file to bed file')
parser.add_argument('-rate_dir', type=str, default="../data/Rates", help='The data directory where Rates and Rates_RandomT stored')

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def convert_rates_mat_into_bed(rate_dir):
    CHROMOSOMES = [str(i) for i in range(1, 23)]
    ltws = []
    for chromosome in CHROMOSOMES:
        print("Started chr%s in %s" % (chromosome, rate_dir))
        data_fp = os.path.join(rate_dir, "%s_chr%s.mat" % ("Rates", chromosome))
        mat_data = sio.loadmat(data_fp)
        inferedRates = mat_data["inferedRates"][:, 0]
        indexs_retained = inferedRates > 1e-6
        fittedSites = mat_data["fittedSites"][indexs_retained]
        inferredMethyFrac = mat_data["inferredMethyFrac"][indexs_retained, 0]
        inferedRates = inferedRates[indexs_retained]
        for sid, site in enumerate(fittedSites):
            ltw = "chr%s\t%d\t%d\t%.4f\t%.4f\n" % (chromosome, site, site + 1, inferedRates[sid], inferredMethyFrac[sid])
            ltws.append(ltw)
    output_fp = os.path.join(rate_dir, "K.bed")
    with open(output_fp, "w") as out_f:
        ltw = "".join(ltws)
        out_f.write(ltw)
        print("Write K.bed Successful!")

if __name__ == "__main__":
    args = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
    else:
        rate_dir = args.rate_dir
        convert_rates_mat_into_bed(rate_dir)
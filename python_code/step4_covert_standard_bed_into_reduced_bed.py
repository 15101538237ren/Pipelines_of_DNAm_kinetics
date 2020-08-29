# -*- coding:utf-8 -*-
import os, argparse, sys

parser = argparse.ArgumentParser(description='Step 4: Covert standard bed format into reduced bed format')
parser.add_argument('-input', type=str, help='The input standard bed file path')
parser.add_argument('-output', type=str, help='The output reduced bed file path')

def covert_standard_bed_into_reduced_bed(input_fp, output_fp, sep="\t"):
    ltws = []
    with open(input_fp, "r") as input_f:
        lines = input_f.readlines()
        print("Started conversion")
        for lid, line in enumerate(lines):
            if lid and (lid + 1) % 100000 == 0:
                print("Parsed %.1f %% of total lines" % ((lid + 1.0) * 10.0 / len(lines)) )
            chri, start_idx, end_idx, val, _ = line.split(sep)
            start_idx = int(start_idx) + 1
            meth, unmeth = val.strip("\'").split("/")
            ltw = "%s\t%d\t%s\t%s\t%s" % (chri, start_idx, end_idx, meth, unmeth)
            ltws.append(ltw)
    ltw = "\n".join(ltws) + "\n"
    with open(output_fp, "w") as output_f:
        output_f.write(ltw)
        print("Converted %s Successful" % input_fp)
if __name__ == "__main__":
    args = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
    else:
        input_fp = args.input
        output_fp = args.output
        covert_standard_bed_into_reduced_bed(input_fp, output_fp)
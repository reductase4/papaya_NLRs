#!/user/bin/python3

#Author: Jiang Qian
#Data: 2022.3.14
#usage: python /data01/jiangqian/RGAugury/RGAugury_pipeline/pfam_out_filter.py -i Panicum_virgatum.pfam.local.search.out -e 0.001 -o Panicum_virgatum.pfam.local.search.out.filter

import argparse

parser = argparse.ArgumentParser(description='filter pfam out by Evalue')
parser.add_argument('-i', '--input', required=True, help="pfam results")
parser.add_argument('-e', '--evalue', required=True, help="Evalue cutoff")
parser.add_argument('-o', '--output', required=True, help="pfam out filter results")
args = parser.parse_args()

Input = args.input
Evalue = float(args.evalue)
Output = args.output

F = open(Input,'r')
O = open(Output,'w')
for line in F:
    line = line.strip()
    e = float(line.split()[12])
    if e < Evalue:
        print(line,file=O)
O.close()
F.close()

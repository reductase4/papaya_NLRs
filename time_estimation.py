#!/usr/bin/python
#Author: Jiang Qian
#Data: 2023.01.19

import argparse,os

parser = argparse.ArgumentParser(description='extract ks values and calculate divergence times (Mya) for DupGens')
parser.add_argument('-i', '--input', required=True, help="DupGens file")
parser.add_argument('-k', '--ks', required=True, help="Ks file")
parser.add_argument('-r', '--rate', required=True, help="mutation rate")
parser.add_argument('-o', '--output', required=True, help="output file")
args = parser.parse_args()

inputFile = os.path.abspath(args.input)
ksFile = os.path.abspath(args.ks)
output = os.path.abspath(args.output)
rate = float(args.rate)

f = open(ksFile,'r')
ksDict={}
for line in f:
	line=line.strip()
	if line.startswith('#'):
		continue
	else:
		pairName=line.split()[0]
		ks=line.split()[2]
		ksDict[pairName]=ks
f.close()

o = open(output,'w')
f = open(inputFile,'r')
for line in f:
	line=line.strip()
	if line.startswith('#'):
		titles=line
		print(titles,'Ks','Time',sep='\t',file=o)
	else:
		dup1=line.split()[0]
		dup2=line.split()[2]
		pair=dup1+'-'+dup2
		ksValue=ksDict[pair]
		ksTime=float(ksValue)*1000/(2*rate)
		print(line,ksValue,str(ksTime),sep='\t',file=o)
o.close()
f.close()






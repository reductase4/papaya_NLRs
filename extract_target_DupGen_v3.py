#!/usr/bin/python
# -*- coding: UTF-8 -*-
'''
wensiyu v1.0 2022.1.5
python */dupgen_finder_more_outgroup/05.focus_result_tidy.py -d */result/Cinnamomum_micranthum -o */result -g foucs_gene.txt
python */dupgen_finder_more_outgroup/05.focus_result_tidy.py -d */result/Cinnamomum_micranthum -o */result -l gene1008,gene10081,gene10114
'''
#modified by Jiangqian on 2022.1.13: 1.add 'key' argument; 2.target_genes_N column
#v2: add one result file--cluster genes
#v3: add gene type in OUT2 (focus_cluster_dupgene_clean.txt) --2022.12.22
import argparse,sys,glob
import os
####传参
parser = argparse.ArgumentParser(description='dupgene_finder result tidy')
parser.add_argument("-d","--dir",required=True,help="dupgene_finder result:*.pairs.stats-unique dir path")
parser.add_argument("-o","--outdir",required=True,help="outdir")
#parser.add_argument("-l","--list",help="focus gene ID list,like -l gene1,gene2,gene 3")
parser.add_argument("-g","--gene",help="focus gene ID file")
parser.add_argument("-k","--key",help="species name")
args = parser.parse_args()

key = args.key
targetFile = os.path.abspath(args.gene)
dup_result_dir=os.path.abspath(args.dir)
outdir=os.path.abspath(args.outdir)

if not os.path.exists(outdir):
    os.system("mkdir -p "+outdir)

#tidy_result=glob.glob(dup_result_dir+"/*.pairs.stats-unique")
tidy_result=glob.glob(dup_result_dir+"/*.pairs.stats")
Tidy_result=open(tidy_result[0],'r')
tidy_d={}
for line in Tidy_result:
    if line.startswith('Types'):
        continue
    linelist=line.strip().split('\t')
    tidy_d[linelist[0]]=linelist[1]
Tidy_result.close()

'''
if args.list != None and args.gene ==None:
    focus_genes=(args.list).split(',')
## just gene file first cols
if args.list == None and args.gene !=None:
    focus_genes=[]
    GENE=open(os.path.abspath(args.gene),'r')
    for line in GENE:
        if line.startswith('#') or len(line.strip())==0:
            continue
        if line.strip().split()[0] not in focus_genes:
            focus_genes.append(line.strip().split()[0])
    GENE.close()
if args.list != None and args.gene !=None:
    print("no focus argument,check !")
    exit()
#print(focus_genes)
'''
focus_genes={}
f=open(targetFile,'r')
for line in f:
    line = line.strip()
    if line.startswith('#') or len(line)==0:
        continue
    geneID=line.split()[0]
    geneType=line.split()[1]
    focus_genes[geneID]=geneType
f.close()

#all_pair_file=glob.glob(dup_result_dir+"/*.pairs-unique")
all_pair_file=glob.glob(dup_result_dir+"/*.pairs")
Pair_OUT=open(outdir + os.sep + key + "_focus_gene.pairs",'w+')
print("Duplicate 1","Location","Duplicate 2","Location","E-value","dup_type","gene_type1","gene_type2",sep='\t',file=Pair_OUT)
IDlist=[]
pairs = []
for pair_file in all_pair_file:
    pair_type=os.path.basename(pair_file).rsplit('.',2)[-2]
    Pair=open(pair_file,'r')
    type_count=0
    for line in Pair:
        if line.startswith('Duplicate'):
            continue
        linelist=line.strip().split('\t')
        gene1=linelist[0]
        gene2=linelist[2]
        ## gene1,gene2 both in focus_gene
        if gene1 in focus_genes.keys() and gene2 in focus_genes.keys():
            pairs.append((gene1,gene2))
            type_count+=1
            IDlist.append(gene1)
            IDlist.append(gene2)
            print(line.strip(),pair_type,focus_genes[gene1],focus_genes[gene2],sep='\t',file=Pair_OUT)
    tidy_d["focus_"+pair_type]=type_count
    
    Pair.close()

Pair_OUT.close()

OUT=open(outdir + os.sep + key + "_target_genes.duplist",'w')
new_IDlist=[]
for id in IDlist:
    if id not in new_IDlist:
        new_IDlist.append(id)
        print(id,file=OUT)
OUT.close()

OUT=open(outdir + os.sep + key + "_focus_dupgene_stat.txt",'w+')
#species=os.path.basename(tidy_result[0]).split(".pairs.stats-unique")[0]
species=os.path.basename(tidy_result[0]).split(".pairs.stats")[0]
#print("species","WGD-pairs","TD-pairs","PD-pairs","TRD-pairs","DSD-pairs","focus-WGD-pairs","focus-TD-pairs","focus-PD-pairs","focus-TRD-pairs","focus-DSD-pairs",sep='\t',file=OUT)
linename=[]
stat=[]
for k in tidy_d.keys():
    linename.append(k)
    stat.append(str(tidy_d[k]))
    #linename += "\t"+k
    #stat += "\t"+str(tidy_d[k])

print("species","\t".join(linename),"target_genes_N",sep='\t',file=OUT)
print(species,"\t".join(stat),len(focus_genes.keys()),sep='\t',file=OUT)

OUT.close()

cluster_dict = {}
n=0
usedID=[]
for gene_pair in list(pairs):
    gene1 = gene_pair[0]
    gene2 = gene_pair[1]
    if gene1 in usedID:
        continue
    n+=1
    cluster_dict[n]=[] #cluster list
    cluster_dict[n].append(gene1)
    cluster_dict[n].append(gene2) 
    usedID.append(gene1)
    usedID.append(gene2)
    for gene in cluster_dict[n]:
        for p in list(pairs):
            if gene in p:
                if p[0] not in cluster_dict[n]:
                    cluster_dict[n].append(p[0])
                    usedID.append(p[0])
                elif p[1] not in cluster_dict[n]:
                    cluster_dict[n].append(p[1])
                    usedID.append(p[1])
                pairs.remove(p)

OUT1=open(outdir + os.sep + key + "_focus_cluster_dupgene.txt",'w+')
OUT2=open(outdir + os.sep + key + "_focus_cluster_dupgene_clean.txt",'w+')
for i in cluster_dict.keys():
    print("cluster"+str(i),"\t".join(cluster_dict[i]),sep = '\t',file =OUT1)
    for gene in cluster_dict[i]:
        print(gene,focus_genes[gene],"cluster"+str(i),sep = '\t',file =OUT2)
OUT1.close()
OUT2.close()


        

















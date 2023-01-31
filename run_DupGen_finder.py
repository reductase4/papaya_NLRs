#!/usr/bin/python
#Jiang Qian
#2022.12.21
import argparse
import os,gzip,re,time,sys

#usage: python /data01/jiangqian/NLR/papaya_NLR/scripts/run_DupGen_finder.py -d /data02/database/Genome/Genome_475 -o results -ts target.txt -os outgroup.txt
parser = argparse.ArgumentParser(description='run DupGen_finder')
parser.add_argument('-d', '--dataDir', required=True, help="input dataDir of gff and pep files ")
parser.add_argument('-o', '--outDir', required=True, help="output result project name")
parser.add_argument('-ts', '--targetSpecies', required=True, help="target species eg Ath")
parser.add_argument('-os', '--outgroup', required=True, help="outgroup species eg Ath")
args = parser.parse_args()
dataDir = os.path.abspath(args.dataDir)
outdir = os.path.abspath(args.outDir)
targetFile = os.path.abspath(args.targetSpecies)
outgroupFile = os.path.abspath(args.outgroup)

gffSuffix = '.gff.gz'
newGffSuffix = '.gff'
pepSuffix = '.longest_transcript.pep.fa'

def read_species(infile):
    f = open(infile,'r')
    species = []
    for line in f:
        line = line.strip()
        species.append(line)
    f.close()
    return species

def log_print(string):
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())) + '\t' + string)
    sys.stdout.flush()

def file_check(file):
    if not os.path.exists(file):
        print(file + " dose not exist!")
        exit()
    elif not os.path.getsize(file):
        print(file + " is empty!")
        exit()

def get_new_gff(species,gff,newGff):
    outfile = open(newGff,'w')
    gz = False
    if gff.endswith(".gz"):
        infile = gzip.open(gff,'rb')
        gz = True
    else:
        infile = open(gff,'r')

    for line in infile:
        if gz:
            line = line.decode()
        line = line.strip()
        content = line.split('\t')
        if line.startswith('#') or len(content) < 9:
            continue
        else:
            features = content[2]
            if features == "gene":
                desc = content[8]
                match = re.match(r'ID=(.*?);.*', desc)
                ID = match.group(1)
                Chr = content[0]
                start = content[3]
                end = content[4]
                print(species + '-' + Chr, ID, start, end, sep='\t', file=outfile)
    outfile.close()
    infile.close()

#mkdir outputDir
if not os.path.exists(outdir):
    os.system("mkdir -p " + outdir)
else:
    print(outdir + ' exists!')

#mkdir workDir
workDir = outdir + os.sep + "work_sh"
if not os.path.exists(workDir):
    os.system("mkdir " + workDir)
else:
    print(workDir + ' exists!')

#read species
targets = read_species(targetFile)
outgroups = read_species(outgroupFile)
speciesList = targets + outgroups

#store path of gff, new gff, and pep files
gffDir = outdir + os.sep + "new_gff"
dataPath = {}
for species in speciesList:
    gff = dataDir + os.sep + species + os.sep + species + gffSuffix
    newGff = gffDir + os.sep + species + newGffSuffix
    pep = dataDir + os.sep + species + os.sep + species + pepSuffix
    if species not in dataPath.keys():
        dataPath[species] = {}
    dataPath[species] = {'gff':gff,'newGff':newGff,'pep':pep}

####################### Step0: data check ######################
log_print("Step0: data checking...")
predataDir = outdir + os.sep + "pre_data"
if not os.path.exists(gffDir):
    os.system("mkdir " + gffDir)

if not os.path.exists(predataDir):
    os.system("mkdir " + predataDir)

#check data
sh = workDir + os.sep + "step0_datapath.txt"
SH = open(sh,'w')
for species in speciesList:
    gff = dataPath[species]['gff']
    newGff = dataPath[species]['newGff']
    pep = dataPath[species]['pep']

    file_check(gff)
    file_check(pep)

    if not os.path.exists(newGff):
        get_new_gff(species,gff,newGff)

    file_check(newGff)

    print(gff,pep,newGff,sep='\n',file=SH)
SH.close()

#prepare gff files for DupGen_finder
sh = workDir + os.sep + "step0_prepare_gff.txt"
SH = open(sh,'w')
print('cd ' + predataDir,file=SH)
for targetSpecies in targets:
    targetGff = dataPath[targetSpecies]['newGff'] 
    cmd = "cp " + targetGff + " ."
    print(cmd,file=SH)
    for outgroup in outgroups:
        outgroupGff = dataPath[outgroup]['newGff']
        mergeGff = targetSpecies + '_' + outgroup + '.gff'
        cmd = "cat " + targetGff + " " + outgroupGff + " > " + mergeGff
        print(cmd,file=SH)
SH.close()
os.system('sh ' + sh)

####################### Step1: makeblastdb ######################
log_print("Step1: makeblastdb running...")
databaseDir = outdir + os.sep + "database"
sh = workDir + os.sep + "step1_makeblastdb.sh"

if not os.path.exists(databaseDir):
    os.system("mkdir " + databaseDir)

SH = open(sh,'w')
print('cd ' + databaseDir,file=SH)
for species in speciesList:
    db1 = databaseDir + os.sep + species + '.phr'
    db2 = databaseDir + os.sep + species + '.pin'
    db3 = databaseDir + os.sep + species + '.psq'
    if os.path.exists(db1) and os.path.exists(db2) and os.path.exists(db3):
        continue
    cmd = "makeblastdb -in " + dataPath[species]['pep'] + " -dbtype prot -title " + species + " -out " + species
    print(cmd,file=SH)
SH.close()
os.system('sh ' + sh)

####################### Step2: run blastp ######################
log_print("Step2: blastp running...")
sh = workDir + os.sep + "step2_blastp.sh"

SH = open(sh,'w')
print('cd ' + predataDir,file=SH)
for targetSpecies in targets:
    query = dataPath[targetSpecies]['pep']
    db = databaseDir + os.sep + targetSpecies
    blastOut = targetSpecies + ".blast"
    if not os.path.exists(predataDir + os.sep + blastOut):
        cmd = "blastp -query " + query + " -db " + db + " -evalue 1e-10 -max_target_seqs 5 -num_threads 30 -outfmt 6 -out " + blastOut
        print(cmd,file=SH)
    for outgroup in outgroups:
        db = databaseDir + os.sep + outgroup
        blastOut = targetSpecies + "_" + outgroup + ".blast"
        if os.path.exists(predataDir + os.sep + blastOut):
            continue
        cmd = "blastp -query " + query + " -db " + db + " -evalue 1e-10 -max_target_seqs 5 -num_threads 30 -outfmt 6 -out " + blastOut
        print(cmd,file=SH)
SH.close()
os.system('sh ' + sh)

####################### Step3: run DupGen_finder ######################
log_print("Step3: DupGen_finder running...")
resultsDir = outdir + os.sep + "results"

if not os.path.exists(resultsDir):
    os.system("mkdir " + resultsDir)

sh = workDir + os.sep + "step3_run_DupGen_finder.sh"
SH = open(sh,'w')
print('cd ' + resultsDir,file=SH)
for targetSpecies in targets:
    outgroupSpecies = ','.join(outgroups)
    DupGen_out = resultsDir + os.sep + targetSpecies
    cmd = "DupGen_finder.pl -i " + predataDir + " -t " + targetSpecies + " -c " + outgroupSpecies + " -o " + DupGen_out
    print(cmd,file=SH)
SH.close()
os.system('sh ' + sh)

log_print("Done")




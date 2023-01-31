# papaya_NLRs

**<left>Description</left>**
---------------
## 1. Identify and classify NBS-type resistance genes in papaya

To identify and classify NBS-type resistance genes in papaya, we modified the <a href="https://bitbucket.org/yaanlpc/rgaugury/src/master/">RGAugury</a> project for identification of resistance gene analogs and their domain-related genes.  

Cite:
Li, P., Quan, X., Jia, G., Xiao, J., Cloutier, S., & You, F. M. (2016). RGAugury: a pipeline for genome-wide prediction of resistance gene analogs (RGAs) in plants. BMC genomics, 17(1), 852. https://doi.org/10.1186/s12864-016-3197-x

Usage: 
1. RGAugury_modified.pl -p longest_transcript.pep.fa
-d SUPERFAMILY,SMART,gene3d,Pfam

## 2. Target duplication gene pairs extraction

extract_target_DupGen_v3.py is used to extract NLRs involved in duplication modes identified by <a href="https://github.com/qiao-xin/DupGen_finder">DupGen_finder</a> and cluster duplicated NLRs into groups.

Cite: 
Qiao, X., Li, Q., Yin, H., Qi, K., Li, L., Wang, R., Zhang, S., & Paterson, A. H. (2019). Gene duplication and evolution in recurring polyploidization-diploidization cycles in plants. Genome biology, 20(1), 38. https://doi.org/10.1186/s13059-019-1650-2

Usage: 
1. python extract_target_DupGen_v3.py -d /DupGen/results/Carica_papaya -g Carica_papaya/Carica_papaya.RGA.candidates.lst -k Carica_papaya -o results

## 3. Estimation of divergent times

We run calc_kaks_4dtv_1sp_dup.sh with prepared cds, pep, and homology-pair files to calculate Ks values. Then estimate divergent times with time_estimation.py.

Cite:
Zhang, Z., Xiao, J., Wu, J., Zhang, H., Liu, G., Wang, X., & Dai, L. (2012). ParaAT: a parallel tool for constructing multiple protein-coding DNA alignments. Biochemical and biophysical research communications, 419(4), 779–781. https://doi.org/10.1016/j.bbrc.2012.02.101

Wang, D., Zhang, Y., Zhang, Z., Zhu, J., & Yu, J. (2010). KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. Genomics, proteomics & bioinformatics, 8(1), 77–80. https://doi.org/10.1016/S1672-0229(10)60008-3

Usage: 
1. sh calc_kaks_4dtv_1sp_dup.sh Carica_papaya

2. python time_estimation.py -i DupGene_input.txt -k all-results_Carica_papaya.txt -r 12 -o DupGene_ks_time_rate12.txt

DupGene_input.txt: output of extract_target_DupGen_v3.py.

all-results_Carica_papaya.txt: output of calc_kaks_4dtv_1sp_dup.sh.

r: mutation rate.
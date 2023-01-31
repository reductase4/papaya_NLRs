#!/bin/bash
set -e
set -u
##files used in this script: .pep, .cds
##usage: nohup sh calc_kaks_4dtv_1sp.sh Carica_papaya &

Species1=$1
# 对参考物种的蛋白序列构建索引
#echo "makeblastdb -in ${Species1}.pep -dbtype prot -out ${Species1}"
#makeblastdb -in ${Species1}.pep -dbtype prot -out ${Species1}
# 将目标物种的蛋白序列与参考物种进行比对，并保留最优匹配结果
#echo "blastp -query ${Species1}.pep -db ${Species1} -evalue 1e-5 -max_target_seqs 2 -num_threads 10 -out ${Species1}.blastp_out.m6 -outfmt 6"
#blastp -query ${Species1}.pep -db ${Species1} -evalue 1e-5 -max_target_seqs 2 -num_threads 10 -out ${Species1}.blastp_out.m6 -outfmt 6
# 提取最优匹配的同源序列基因对
#echo "cut -f1-2 ${Species1}.blastp_out.m6|sort|uniq|awk '{if ($1 != $2) print $0}' > ${Species1}.homolog"
#cut -f1-2 ${Species1}.blastp_out.m6|sort|uniq|awk '{if ($1 != $2) print $0}' > ${Species1}.homolog
# 合并目标物种和参考物种的蛋白序列和cds序列
#cat ${Species1}.cds ${Species2}.cds >${Species1}_${Species2}.cds
#cat ${Species1}.pep ${Species2}.pep >${Species1}_${Species2}.pep
# 使用ParaAT程序将蛋白序列比对结果转化为cds序列比对结果
echo "12" >proc
echo "ParaAT.pl -h ${Species1}.homolog -n ${Species1}.cds -a ${Species1}.pep -p proc -m muscle -f axt -o ${Species1}_out"
ParaAT.pl -h ${Species1}.homolog -n ${Species1}.cds -a ${Species1}.pep -p proc -m muscle -f axt -o ${Species1}_out
echo "cd ${Species1}_out"
cd ${Species1}_out
# 使用KaKs_Calculator计算kaks值
for i in `ls *.axt`;do KaKs_Calculator -i $i -o ${i}.kaks -m YN;done
# 将多行axt文件转换成单行
for i in `ls *.axt`;do axt2one-line.py $i ${i}.one-line;done
# 使用calculate_4DTV_correction.pl脚本计算4dtv值
ls *.one-line|while read id;do calculate_4DTV_correction.pl $id >${id%%one-line}4dtv;done
# 合并所有同源基因对的4dtv值
for i in `ls *.4dtv`;do awk 'NR>1{print $1"\t"$2}' $i >>all-4dtv.txt;done
# 合并所有同源基因对的kaks值
for i in `ls *.kaks`;do awk 'NR>1{print $1"\t"$3"\t"$4"\t"$5}' $i >>all-kaks.txt;done
# 排序并去冗余
sort all-4dtv.txt|uniq >all-4dtv.results
sort all-kaks.txt|uniq >all-kaks.results
# 清除中间文件
#rm *one-line
#rm all-4dtv.txt
#rm all-kaks.txt
# 将kaks结果文件和4dtv结果文件进行合并
join -1 1 -2 1 all-kaks.results all-4dtv.results > all-results_${Species1}.txt
# 给结果文件添加标题
sed -i '1i\Seq\tKa\tKs\tKa/Ks\t4dtv_corrected' all-results_${Species1}.txt

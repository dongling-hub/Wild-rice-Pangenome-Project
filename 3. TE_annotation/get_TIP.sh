#!/bin/bash

# Ensure sample and thread variables are set
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Please set the sample (first positional argument) and thread variables."
    exit 1
fi

sample=$1
thread=$2

# Identify INS-TIP
grep -v "#" ${sample}.INS.vcf | awk '{print ">"$3"\n"$5}' > ${sample}.INS.fa
grep -v "#" ${sample}.INS.vcf | awk '{print $3"\t"${sample}"\t"$2"\t"$2}' > ${sample}.INS.bed
blastn -db panEDTA.TElib.fa -query ${sample}.INS.fa -out ${sample}.INS.panTE.blastn -num_threads ${thread} -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore' -perc_identity 80 -word_size 50 -evalue 1e-10 -max_target_seqs 1
awk '{if ($4*${sample}2/${sample}3>=80) print ${sample}"\t"$2}' ${sample}.INS.panTE.blastn | uniq | awk -F '[\t|#]' '{print ${sample}"\t"$2"\t"$3}' > ${sample}.INS.panTE.txt
python combine_rows.py ${sample}.INS.panTE.txt ${sample}.INS.bed -o ins
awk '{if (NF==8) print $6"\t"$7"\t"$8"\t"$3"\t"$4}' ins > ${sample}.INS.TIP.txt

# Identify DEL-TIP
grep -v "#" ${sample}.DEL.vcf | awk '{print ">"$3"\n"$4}' > ${sample}.DEL.fa
grep -v "#" ${sample}.DEL.vcf | awk -F '[;|\t]' '{print $3"\t"${sample}"\t"$2"\t"$(NF-5)}' | sed 's/END=//g' > ${sample}.DEL.bed
blastn -db panEDTA.TElib.fa -query ${sample}.DEL.fa -out ${sample}.DEL.panTE.blastn -num_threads ${thread} -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore' -perc_identity 80 -word_size 50 -evalue 1e-10 -max_target_seqs 1
awk '{if ($4*${sample}2/${sample}3>=80) print ${sample}"\t"$2}' ${sample}.DEL.panTE.blastn | uniq | awk -F '[\t|#]' '{print ${sample}"\t"$2"\t"$3}' > ${sample}.DEL.panTE.txt
python combine_rows.py ${sample}.DEL.panTE.txt ${sample}.DEL.bed -o del
awk '{if (NF==8) print $6"\t"$7"\t"$8"\t"$3"\t"$4}' del > ${sample}.DEL.TIP.txt
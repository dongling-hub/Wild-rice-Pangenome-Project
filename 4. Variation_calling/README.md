# Variation calling

- [Variation calling](#variation-calling)
  - [SNPs calling](#snps-calling)
  - [SV calling](#sv-calling)
  - [Indels calling](#indels-calling)

## SNPs calling

```shell
# SNPs calling using longshot
# For HiFi reads
minimap2 -a -x map-hifi -t ${thread} ${reference} ${sample}.hifi.fastq --MD -o ${sample}.sam
samtools sort -@ ${thread} -o ${sample}.bam ${sample}.sam
samtools index ${sample}.bam
longshot --bam ${sample}.bam -c 3 -D 3:10:50 --ref ${reference} --out ${sample}.vcf

# For ONT reads
minimap2 -a -x map-ont -t ${thread} ${reference} ${sample}.hifi.fastq --MD -o ${sample}.sam
samtools sort -@ ${thread} -o ${sample}.bam ${sample}.sam
samtools index ${sample}.bam
longshot --bam ${sample}.bam -c 10 -D 3:10:50 --ref ${reference} --out ${sample}.vcf

# SNPs calling using show-snps
nucmer -p ${sample}_Nip -t ${thread} ${reference} ${sample}.genome.fa
delta-filter -1 ${sample}_Nip.delta >${sample}_Nip.filtered.delta
show-snps -C ${sample}_Nip.filtered.delta -T >${sample}_Nip.filtered.delta_C_tab.vcf
show-snps -C -I ${sample}_Nip.filtered.delta -T >${sample}_Nip.filtered.delta_C_SNP_tab.vcf
```

## SV calling

- pbsv

```shell
#for HiFi reads
pbmm2 align ${reference} ${sample}.reads.fasta ${sample}.reads_pbmm2.bam --sort --preset CCS --sample ${sample} -j 5 -J 5
pbsv discover ${sample}.reads_pbmm2.bam ${sample}.svsig.gz -s ${sample} --tandem-repeats ${reference}.trf.bed
pbsv call ${reference} ${sample}.svsig.gz ${sample}_pbsv.vcf --min-sv-length 30 --max-ins-length 100K --max-dup-length 100K -j 10 --ccs
awk -F '[;\t]'  '{if ($1 ~ /^#/) {print $0} else {if ($8=="SVTYPE=INS") {print $0}}}' ${sample}_pbsv.vcf >${sample}_pbsv_INS.vcf
awk -F '[;\t]'  '{if ($1 ~ /^#/) {print $0} else {if ($8=="SVTYPE=DEL") {print $0}}}' ${sample}_pbsv.vcf >${sample}_pbsv_DEL.vcf
```

- CuteSV

```shell
#for HiFi reads
minimap2 -a -x map-hifi -t ${thread} ${reference} ${sample}.reads.fasta --MD -o ${sample}.reads_minimap2.sam
samtools sort -@ ${thread} -o ${sample}.reads_minimap2.sorted.bam ${sample}.reads_minimap2.sam
samtools index ${sample}.reads_minimap2.sorted.bam
cuteSV ${sample}.reads_minimap2_sorted.bam ${reference} ${sample}.cuteSV.vcf ${path} --sample ${sample}.cuteSV --min_size 30 --max_size 100000 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_support 3 -t ${thread}

#for ONT reads
minimap2 -a -x map-ont -t ${thread} ${reference} ${sample}.reads.fasta --MD -o ${sample}.reads_minimap2.sam
samtools sort -@ ${thread} -o ${sample}.reads_minimap2.sorted.bam ${sample}.reads_minimap2.sam
samtools index ${sample}.reads_minimap2.sorted.bam
cuteSV ${sample}.reads_minimap2_sorted.bam ${reference} ${sample}.cuteSV.vcf ${path} --sample ${sample}.cuteSV  --min_size 30 --max_size 100000 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --min_support 10 -t ${thread}

awk -F '[;\t]'  '{if ($1 ~ /^#/) {print $0} else {if ($9=="SVTYPE=INS") {gsub(/END=[0-9]+/,"END="$2) ;print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}}' ${sample}.cuteSV.vcf >${sample}.cuteSV_INS.vcf #cutsSV从0开始计数，需要将坐标-1
awk -F '[;\t]'  '{if ($1 ~ /^#/) {print $0} else {if ($9=="SVTYPE=DEL") {gsub(/END=[0-9]+/,"END="$2) ;print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}}' ${sample}.cuteSV.vcf >${sample}.cuteSV_DEL.vcf #cutsSV从0开始计数，需要将坐标-1
```

- SVIM-asm

```shell
minimap2 -a -x asm5 --cs -r 2k -t ${thread} ${reference} ${sample}.contig.fasta -o ${sample}.contig_minimap2.sam
samtools sort -@ ${thread} -o ${sample}.contig_minimap2.sorted.bam ${sample}.contig_minimap2.sam
samtools index ${sample}.contig_minimap2.sorted.bam
svim-asm haploid ${path} ${sample}.contig_minimap2.sorted.bam ${reference} --min_sv_size 30 --max_sv_size 100000 --sample ${sample}.SVIM_asm
awk -F '[;\t]'  '{if ($1 ~ /^#/) {print $0} else {if ($8=="SVTYPE=INS") {print $0}}}' variants.vcf >${sample}.SVIM_asm_INS.vcf
awk -F '[;\t]'  '{if ($1 ~ /^#/) {print $0} else {if ($8=="SVTYPE=DEL") {print $0}}}' variants.vcf >${sample}.SVIM_asm_DEL.vcf
```

- SyRI

```shell
nucmer -p ${sample}_Nip -t ${thread} ${reference} ${sample}.transfer.merge.chr.fasta
delta-filter -m ${sample}_Nip.delta >${sample}_Nip.filtered.delta
show-coords -THrd ${sample}_Nip.filtered.delta >${sample}_Nip.filtered.coords
syri -c ${sample}_Nip.filtered.coords -d ${sample}_Nip.filtered.delta -r ${reference} -q ${sample}.transfer.merge.chr.fasta --prefix ${sample}.syri -s show-snps
python3.8 plotsr ${sample}.syri.out ${reference} ${sample}.transfer.merge.chr.fasta -H 8 -W 5

#inversion
grep "<INV>" ${sample}.syri.vcf |awk '{print $0"\tGT\t1/1"}' |bgzip >${sample}.syri.INV.vcf.gz
tabix -p vcf ${sample}.syri.INV.vcf.gz

#translocation
grep -P "<TRANS>|<INVTR>" ${sample}.syri.vcf |awk '{print $0"\tGT\t1/1"}' |bgzip >${sample}.syri.TRANS.vcf.gz
tabix -p vcf ${sample}.syri.TRANS.vcf.gz
```

- SURVIVOR

```shell
#insertion
#for HiFi assemblies
echo -e "${sample}_pbsv_INS.vcf\n${sample}.cuteSV_INS.vcf\n${sample}.SVIM_asm_INS.vcf" >${sample}_INS.vcf.files
SURVIVOR merge ${sample}_INS.vcf.files 50 2 1 1 0 30 ${sample}_INS.vcf #INS的merge参数为50bp
#for ONT assemblies
echo -e "${sample}.cuteSV_INS.vcf\n${sample}.SVIM_asm_INS.vcf" >${sample}_INS.vcf.files
SURVIVOR merge ${sample}_INS.vcf.files 50 1 1 1 0 30 ${sample}_INS.vcf #INS的merge参数为50bp

awk '/^##/ {print $0} /^#CHROM/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10} !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT\t1/1"}' ${sample}_INS.vcf | bcftools norm -c x -f ${reference} -d all -Ov --threads ${thread} | bcftools sort -Ov -o ${sample}_INS_norm.vcf #修改格式，只取基因型；矫正坐标，去除不正确的坐标

#deletion
#for HiFi assemblies
echo -e "${sample}_pbsv_DEL.vcf\n${sample}.cuteSV_DEL.vcf\n${sample}.SVIM_asm_DEL.vcf" >${sample}_DEL.vcf.files
SURVIVOR merge ${sample}_DEL.vcf.files 0.5 2 1 1 0 30 ${sample}_INS.vcf #DEL的merge参数为length的50%
#for ONT assemblies
echo -e "${sample}.cuteSV_DEL.vcf\n${sample}.SVIM_asm_DEL.vcf" >${sample}_DEL.vcf.files
SURVIVOR merge ${sample}_DEL.vcf.files 0.5 1 1 1 0 30 ${sample}_INS.vcf #DEL的merge参数为length的50%

awk '/^##/ {print $0} /^#CHROM/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10} !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT\t1/1"}' ${sample}_INS.vcf | bcftools norm -c x -f ${reference} -d all -Ov --threads ${thread} | bcftools sort -Ov -o ${sample}_INS_norm.vcf #修改格式，只取基因型；矫正坐标，去除不正确的坐标
```

## Indels calling

```shell
grep -v "#" ${sample}.syri.vcf |grep -P "INS|DEL" |grep -v -P "<INS|<DEL>" |awk '{if ((length($5)-length($4))*(length($5)-length($4)) < 900) print $1"\t"$2"\t"tolower(substr($3,1,3))"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT\t1/1"}' >${sample}.syri.InDel.vcf
```

# Indica-japonica Differentiated Analysis

- [Indica-japonica Differentiated Analysis](#indica-japonica-differentiated-analysis)
  - [1. The identify of indica-japonica differentiated variations](#1-the-identify-of-indica-japonica-differentiated-variations)
  - [2. Ascertain the origin of indica–japonica differentiated variations](#2-ascertain-the-origin-of-indicajaponica-differentiated-variations)

## 1. The identification of indica-japonica differentiated variations

```shell
vcftools --vcf ${merge}.vcf.gz --keep ${japonica}.list --freq --out ${japonica}
vcftools --vcf ${merge}.vcf.gz --keep ${indica}.list --freq --out ${indica}
vcftools --vcf ${merge}.vcf.gz --keep ${Or-IIIa}.list --freq --out ${Or-IIIa}
vcftools --vcf ${merge}.vcf.gz --keep ${Or-Ia}.list --freq --out ${Or-Ia}

# Filter SNPs with greater than 0.9 for japonica and indica
grep -v nan ${japonica}.frq | awk -v per=0.9 -F '[:|\t]' '{if ($6>=per) print $1"\t"$2"\t"$5 ;else if ($8>=per) print $1"\t"$2"\t"$7}' > ${japonica}.90.SNP
grep -v nan ${indica}.frq | awk -v per=0.9 -F '[:|\t]' '{if ($6>=per) print $1"\t"$2"\t"$5 ;else if ($8>=per) print $1"\t"$2"\t"$7}' > ${indica}.90.SNP

# Identify indica–japonica differentiated SNPs
python combine_rows.py ${japonica}.90.SNP ${indica}.90.SNP -o ${japonica-indica}
awk '{if (NF==8) print $0}' ${japonica-indica} | cut -f 1,2,5,8 | sed '1d' | awk '{if ($3!=$4) print $0}' | sed '1i CHROM\tPOS\tjaponica\tindica' > ${japonica-indica}-diff.SNP_90.txt
```

## 2. Ascertain the origin of indica–japonica differentiated variations

Table 1. Classification of origins of indica-japonica differentiated SNPs

<img src="./SupplementaryTable20.png">

```shell
python combine_rows.py ${japonica-indica}-diff.SNP_90.txt ${Or-IIIa}.frq ${Or-Ia}.frq

awk '{if (NF==18) print $0}' output.txt | \
cut -f 3-6,11,12,17,18 | \
grep -v "nan" | \
tr ":" "\t" | \
awk '{if ($6>=$8) print $0; else print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12}' | \
awk '{if ($10>=$12) print $0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$11"\t"$12"\t"$9"\t"$10}' | \
sed '1i Chr\tPos\tjaponica\tindica\tOr-IIIa_major\tfrequency\tOr-IIIa_alt\tfrequency\tOr-Ia_major\tfrequency\tOr-Ia_alt\tfrequency' > indica-japonica-diff.frq

# Divide the ancestral origin of Indica–Japonica Differentiated variations into 6 categories
# ancestor differentiation
awk -v per=0.6 '{if ($5!=$9 && $3==$5 && $4==$9 && $6>=per && $10>=per) print $1"\t"$2}' indica-japonica-diff.frq > ancestor-diff.list
# indica novel mutations
awk -v per=1 '{if ($5==$9 && $3==$9 && $6>=per && $10>=per) print $1"\t"$2}' indica-japonica-diff.frq > indica.new.list
# japonica novel mutations
awk -v per=1 '{if ($5==$9 && $4==$9 && $6>=per && $10>=per) print $1"\t"$2}' indica-japonica-diff.frq > japonica.new.list
# japonica preference
awk -v per=0.6 '{if ($5==$9 && $4==$9 && $6>=per && $10>=per) print $0}' indica-japonica-diff.frq | \
awk '{if ($6!=1 || $10!=1) print $1"\t"$2}' > japonica.preference.list
# indica preference
awk -v per=0.6 '{if ($5==$9 && $3==$9 && $6>=per && $10>=per) print $0}' indica-japonica-diff.frq | \
awk '{if ($6!=1 || $10!=1) print $1"\t"$2}' > indica.preference.list
```

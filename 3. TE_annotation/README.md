# TE Annotation

This document provides instructions and information on the TE (Transposable Element) annotation and analysis process.

- [TE Annotation](#te-annotation)
  - [1. TE annotation](#1-te-annotation)
  - [2. Construct a pan-genome TE library](#2-construct-a-pan-genome-te-library)
  - [3. The identify of Or-IIIa-lager Gypsy families](#3-the-identify-of-or-iiia-lager-gypsy-families)
  - [4. Analysis of the LTR family dynamics](#4-analysis-of-the-ltr-family-dynamics)
  - [5. The identification of TE insertion polymorphism (TIP)](#5-the-identification-of-te-insertion-polymorphism-tip)
  - [6. The identification of centromeres and telomeres](#6-the-identification-of-centromeres-and-telomeres)

## 1. TE annotation

```shell
perl EDTA.pl --genome ${sample}.chr.fa --species Rice --overwrite 0 --cds ${sample}.evm.out.cds --curatedlib rice6.9.5.liban --sensitive 1 --anno 1 --evaluate 0 --threads ${thread}
```

## 2. Construct a pan-genome TE library

```shell
sh panEDTA_construct.sh -g genome_list.txt -t ${thread}
sh panEDTA_reanno.sh -g genome_list.txt -t ${thread}
```

## 3. The identify of Or-IIIa-lager Gypsy families

```shell
grep Gypsy panTE.anno.TEfamily.txt | cut -f 2 >panTE.anno.Gypsyfamily.list
python extract_rows.py panTE.anno.TEfamily.bygroup.bp.sum panTE.anno.Gypsyfamily.list panTE.anno.Gypsyfamily.bygroup.bp.sum keep_header=true

#Or-IIIa/Os.japonica
# Filter and process Gypsy family data
cut -f 1,70-100,135-148 panTE.anno.Gypsyfamily.bygroup.bp.sum | awk '{sum=0; for(i=2;i<=NF;i++) sum+=$i; if(sum != 0) print $0}' | awk '{count=0; for(i=2;i<=NF;i++) if($i != 0) count++; if(count != 1) print $0}' > IIIa_japonica-Gypsyfamily.diff.bp.sum  # Remove family specific to varieties

# Calculate average length for each family in Or and Os groups
sed '1d' IIIa_japonica-Gypsyfamily.diff.bp.sum | awk '{sumOr=0; sumOs=0; for(i=2;i<=32;i++) sumOr+=$i; for(i=33;i<=46;i++) sumOs+=$i; avgOr=sumOr/31; avgOs=sumOs/14; print $1"\t"avgOr"\t"avgOs;}' | sed '1i TE_fam\tOr\tOs' > IIIa_japonica-Gypsyfamily.diff.bp.avg
$cat <(head -1 IIIa_japonica-Gypsyfamily.diff.bp.avg) <(sed '1d' IIIa_japonica-Gypsyfamily.diff.bp.avg  |sed 's/_INT//g;s/_LTR//g' |awk '{Or[$1]+=$2;Os[$1]+=$3}END{for(c in Or){print c"\t"Or[c]"\t"Os[c]}}') >IIIa_japonica-Gypsyfamily.merge.diff.avg

# List Gypsy superfamilies where Or is larger by at least 250000
awk '{if ($2-$3>=250000) print $1}' IIIa_japonica-Gypsyfamily.merge.diff.avg  >Or-IIIa.Gypsyfamily-lager.list
```

## 4. Analysis of the LTR family dynamics

```shell
# Get LTR information from panEDTA
perl ~/panEDTA_script/PopTEvo-main/TE_annotation/bin/find_LTR.pl -lib <(awk '{print ${sample}}' panEDTA.TElib.fa) >panEDTA.TElib.LTR.info

# Identify solo LTRs
perl ~/panEDTA_script/PopTEvo-main/TE_annotation/bin/solo_finder.pl -i ${sample}.transfer.merge.fasta.mod.out -info panEDTA.TElib.LTR.info > ${sample}.transfer.merge.fasta.mod.out.solo.LTR.list

# Identify intact LTRs
grep struc ${sample}.transfer.merge.fasta.mod.EDTA.TEanno.gff3 | grep LTR_retrotransposon | perl -nle 'my ($chr, $str, $end, $info) = (split)[0,3,4,8]; my ($id, $iden) = (${sample}, $2) if $info =~ /Name=(.*);Classification.*ltr_identity=([0-9.]+);/; print "$chr\t$str\t$end\t$id\t$iden"' > ${sample}.transfer.merge.fasta.intact.LTR.bed

# Summarize solo and intact LTRs
for i in `cat Or-IIIa.Gypsyfamily-lager.list`; do
    solo=$(sed 's/_INT//g;s/_LTR//g' ${sample}/*solo.LTR.list | grep $i | wc -l)
    intact=$(sed 's/_INT//g;s/_LTR//g' ${sample}/*intact.LTR.bed | grep $i | wc -l)
    ratio=`echo "scale=3;$solo/$intact" | bc`
    size=$(sed 's/_INT//g;s/_LTR//g' ${sample}/*mod.EDTA.TEanno.sum.fam | grep $i | awk '{sum+=$3} END {print sum}')
    echo -e $i"\t"$solo"\t"$intact"\t"$ratio"\t"$size
done | sed '1i Gypsyfamily\tsolo\tintact\tsolo:ratio\tSize' > ${sample}.solo_intact.txt

# Perform t-test
for Or in ${Or-IIIa}.list; do
    cat $Or.solo_intact.txt | awk '{print "'$Or'""\tOr-IIIa\t"$0}'
done | grep -v "Gypsyfamily" | awk '{if (NF==7) print $0}' | cut -f 1-3,6-7 | sed '1i Accesstion\tGroup\tGypsyfamily\tRatio\tSize' > Or-IIIa-japonica.ratio-size.summary

for Os in ${japonica}.list; do
    cat $Os.solo_intact.txt | awk '{print "'$Os'""\tjaponica\t"$0}'
done | grep -v "Gypsyfamily" | awk '{if (NF==7) print $0}' | cut -f 1-3,6-7 >> Or-IIIa-japonica.ratio-size.summary
```

```r
library(dplyr)

# Read data
data <- read.table("Or-IIIa-japonica.ratio-size.summary", header = TRUE)

# Split data into two groups
group_or_iiia <- filter(data, Group == "Or-IIIa")
group_os_japonica <- filter(data, Group == "japonica")

# Ensure each Gypsyfamily has data in both groups
common_families <- intersect(group_or_iiia$Gypsyfamily, group_os_japonica$Gypsyfamily)

# Initialize a dataframe to store Gypsyfamily and corresponding p-values
results_df <- data.frame(Gypsyfamily = character(), Ratio = numeric(), Size = numeric(), Diff = numeric(), stringsAsFactors = FALSE)

for (family in common_families) {
    # Extract data for each family
    or_iiia_ratios <- filter(group_or_iiia, Gypsyfamily == family)$Ratio
    os_japonica_ratios <- filter(group_os_japonica, Gypsyfamily == family)$Ratio
    or_iiia_size <- filter(group_or_iiia, Gypsyfamily == family)$Size
    os_japonica_size <- filter(group_os_japonica, Gypsyfamily == family)$Size

    # Initialize variables to store p-values
    ratio_p_value <- NA
    size_p_value <- NA

    # Ensure each family has at least one data point in both groups
    if (length(or_iiia_ratios) > 0 & length(os_japonica_ratios) > 0) {
        # Perform t-test
        ratio_result <- t.test(or_iiia_ratios, os_japonica_ratios)
        ratio_p_value <- ratio_result$p.value
    }

    # Ensure each family has at least one data point in both groups
    if (length(or_iiia_size) > 0 & length(os_japonica_size) > 0) {
        # Perform t-test
        size_result <- t.test(or_iiia_size, os_japonica_size)
        size_p_value <- size_result$p.value
    }

    # Add family name and p-values to the results dataframe
    results_df <- rbind(results_df, data.frame(Gypsyfamily = family, Ratio = ratio_p_value, Size = size_p_value, Diff = mean(or_iiia_size) - mean(os_japonica_size)))
}

# Replace NaN values in Ratio column with 1
results_df$Ratio <- ifelse(results_df$Ratio == "NaN", 1, results_df$Ratio)

# Plot the results
ggplot(results_df, aes(x = -log10(Size), y = -log10(Ratio), size = Diff / 1e+6, color = ifelse(-log10(Ratio) > -log10(0.01), "Removal", "Amplification"))) +
    geom_point() +
    scale_size(range = c(1, 10)) + # Adjust the size range of points
    scale_color_manual(values = c("Removal" = "#00B050", "Amplification" = "#00B0F0")) +
    theme_bw() + xlim(0, 20) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", size = 2) +
    geom_vline(xintercept = -log10(0.01), linetype = "dashed", color = "black", size = 2) +
    theme(legend.position = 'none')
```

## 5. The identification of TE insertion polymorphism (TIP)

```shell
#Identify TIP
sh get_TIP.sh ${sample} ${thread}

ls *.DEL.TIP.txt > DEL.TIP.list
ls *.INS.TIP.txt > INS.TIP.list
# Extract TE superfamily type
python merge.TIP.py DEL.TIP.list DEL.TIP.superfamily.txt
python merge.TIP.py INS.TIP.list INS.TIP.superfamily.txt

# Merge adjacent TIPs that meet the requirements
python combines.adjacent.TIP.py DEL.TIP.superfamily.txt 0.3 DEL.TIP.superfamily.merge.txt
python combines.adjacent.TIP.py INS.TIP.superfamily.txt 50 INS.TIP.superfamily.merge.txt

# Count the types of TIPs (the most common type in the population)
python sta.TIP_type.py DEL.TIP.superfamily.merge.txt >DEL.TIP.merge.superfamily.list
python sta.TIP_type.py INS.TIP.superfamily.merge.txt >INS.TIP.merge.superfamily.list
```

## 6. The identification of centromeres and telomeres

```shell
# Telomere
perl telomere_finder.pl --repeat-unit TTTAGGG ${sample}.genome.fasta >${sample}_telomere.info

# Centromere
clustalw2 -INFILE=rice_CentO_2002.fasta -TYPE=DNA -OUTFILE=rice_CentO_2002.aln
hmmbuild --dna rice_CentO_2002.hmm rice_CentO_2002.stockholm
nhmmer -o ${sample}_CentO.out --tblout ${sample}_CentO.tblout -E 1e-5 --cpu ${thread} rice_CentO_2002.hmm ${sample}.genome.fasta
mafft --maxiterate 1000 --globalpair --thread ${thread} all_162_CentO_155_cons.fa >all_162_CentO_155_cons.mafft.accuracy.aln.fa
iqtree -s all_162_CentO_155_cons.mafft.accuracy.aln.fa -nt AUTO -b 100 -pre all_162_CentO_155_cons.mafft.accuracy.aln
```

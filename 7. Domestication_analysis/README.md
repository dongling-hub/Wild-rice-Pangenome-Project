# Domestication analysis

- [Domestication analysis](#domestication-analysis)
  - [1. The identification of genomic selective sweep](#1-identification-of-genomic-selective-sweep)
  - [2. The identification domesticated-related PAVs (domPAVs)](#2-the-identification-domesticated-related-pavs-dompavs)
  - [3. Analysis of domestication genes](#3-analysis-of-domestication-genes)

## 1. The identification of genomic selective sweep

```shell
# piRatio
vcftools --gzvcf ${merge}.vcf.gz --keep ${Or-IIIa_Or-Ib}.list --window-pi 100000 --window-pi-step 10000 --out ${Or-IIIa_Or-Ib}
vcftools --gzvcf ${merge}.vcf.gz --keep ${Os}.list --window-pi 100000 --window-pi-step 10000 --out ${Os}
python combine_rows.py ${Or-IIIa_Or-Ib}.windowed.pi ${Os}.windowed.pi -o ${Or-IIIa_b_Os}.merge
sed '1d' ${Or-IIIa_b_Os}.merge |awk '{if (NF==12 && !($7<0.0005 && $12<0.0005) && $7/$12>2.45485) print $3"\t"$4"\t"$5}' >${Or-IIIa_b_Os}.merge.piRatio.bed

# Fst
vcftools --gzvcf ${merge}.vcf.gz --weir-fst-pop ${Or-IIIa_Or-Ib}.list --weir-fst-pop ${Os}.list --fst-window-size 100000 --fst-window-step 10000 --out ${Or-IIIa_Or-Ib}.merge
sed '1d' ${Or-IIIa_Or-Ib}.merge.windowed.weir.fst |awk '{if ($5>=0.305923) print $1"\t"$2"\t"$3}' >${Or-IIIa_Or-Ib}.merge.Fst.bed

# Intersect Fst and piRatio results
bedtools intersect -a ${Or-IIIa_Or-Ib}.merge.Fst.bed -b ${Or-IIIaIb_Os}.merge.piRatio.bed | bedtools sort | bedtools merge -d 30000 > ${Or-IIIa_Or-Ib}.merge.intersect.bed #allowing gaps of less than 30k
```

## 2. The identification domesticated-related PAVs (domPAVs)

We conducted a two-sided Fisher's test comparing wild and cultivated rice, considering PAVs with an FDR-adjusted p values < 0.05 as significant.

```shell
head -n4 ${Or-IIIaIb_Os}.INS.sta.txt
ID      Or_Ref  Or_Alt  Os_Ref  Os_Alt
INS_1   49      1       43      0
INS_2   50      0       43      0
INS_3   50      0       43      0

python3 fisher.py ${Or-IIIaIb_Os}.INS.sta.txt ${Or-IIIaIb_Os}.INS.FDR.txt
```

## 3. Analysis of domestication genes

We extracted gene sequences using `GMAP` and aligned these sequences using `MAFFT`.

```shell
# Build the GMAP database for the sample genome
gmap_build --dir ${gmap_db} --genomedb ${sample} ${sample}.genome.fa -t ${thread} -s

# Align the reference CDS to the sample genome and output in GFF3 format
gmap --dir=${path} -d ${gmap_db} -f gff3_gene ${gene.reference.cds} -t ${thread} > ${gene}.${sample}.gff3

# Extract CDS sequences from the GFF3 file
python gff3_to_cds.py ${gene}.${sample}.gff3 ${sample}.genome.fa ${gene}.${sample}.cds.fa
cat ${gene}*.cds.fa >${gene}.total.cds.fa

# Perform multiple sequence alignment on the concatenated CDS sequences
mafft --maxiterate 1000 --thread ${thread} --localpair ${gene}.total.cds.fa  >${gene}.total.cds.mafft.fa
```

- Haplotype analysis

```r
library(geneHapR)
AccINFO <- import_AccINFO("group.list")
seqs <- import_seqs(paste0(geneID,".total.cds.mafft.fa"))
hapResult <- seqs2hap(seqs,hetero_remove = TRUE, na_drop = FALSE,Ref="Nip_hifi",maxGapsPerSeq = 0.5)
hapSummary <- hap_summary(hapResult)
write.hap(hapSummary, file = paste0(geneID,".hapSummary"))
hap <- filter_hap(hapSummary,rm.mode = c("freq"),freq.min = 2)

# Display haplotypes in a table format
pdf(
  file = paste0(geneID,".hapTable.pdf"),
  width  = 20,
  height = 15,
)
plotHapTable(hap,
             angle = 0,
             displayIndelSize = 3,
             replaceMultiAllele = TRUE,
             title = paste0(geneID,".gene"))
dev.off()

# Display haplotypes in a network graph
hapNet <- get_hapNet(hap,AccINFO = AccINFO,groupName = "Group")
pdf(
  file = paste0(geneID,"hapNet.pdf"),
  width = 30,
  height = 50
)
plotHapNet(hapNet,
            size = "freq",
            scale = "log2",
            cex = 1,
            col.link = 1,
            link.width = 2,
            show.mutation = 2,
            pie.lim=c(3,1),
            labels.cex=2.5,
            labels.col="black",
            show_size_legend=TRUE,
            show_color_legend=FALSE,
)
dev.off()
```

- Gene phylogenetic tree construction 

```shell
python gff3_to_gene.py ${gene}.${sample}.gff3 ${sample}.genome.fa ${gene}.${sample}.gene.fa
cat ${gene}*.gene.fa >${gene}.total.gene.fa
mafft --maxiterate 1000 --thread ${thread} --localpair ${gene}.total.gene.fa  >${gene}.total.gene.mafft.fa
iqtree2 -s ${gene}.total.gene.fa --seqtype AA  --prefix ${gene}.total.gene --seed 12345 -T AUTO -B 1000 -m MFP
```

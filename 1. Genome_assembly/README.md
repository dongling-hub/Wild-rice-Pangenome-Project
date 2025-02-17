# Genome assembly and analysis

- [Genome assembly and analysis](#genome-assembly-and-analysis)
  - [1. Genome assembly](#1-genome-assembly)
  - [2. Extract nucleus and organelles genomes](#2-extract-nucleus-and-organelles-genomes)
  - [3. Evaluation of genome assemblies](#3-evaluation-of-genome-assemblies)
  - [4. The organelle genomes de novo assembly](#4-the-organelle-genomes-de-novo-assembly)

## 1. Genome assembly

- PacBio HiFi assembly

```shell
bedtools bamtofastq -i ${sample}.bam -fq ${sample}.hifi.fastq
hifiasm-0.16.0/hifiasm -t ${thread} --primary -o ${sample} ${sample}.hifi.fastq
awk '/^S/{print ">"$2;print $3}' ${sample}.p_ctg.gfa >${sample}.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${sample}.a_ctg.gfa >${sample}.a_ctg.fa
```

- Nanopore assembly

```shell
#Necat
necat.pl config ${sample}_config.txt
vim read_list.txt
necat.pl correct ${sample}_config.txt
necat.pl assemble ${sample}_config.txt
necat.pl bridge ${sample}_config.txt

#NextDenovo
nextDenovo run.cfg
minimap2 -t ${thread} -x map-ont nextgraph.assembly.contig.fasta ${sample}.fastq >${sample}.paf
racon -t ${thread} ${sample}.fastq ${sample}.paf nextgraph.assembly.contig.fasta >${sample}.consensus.fasta
#NextPolish
vi sgs.fofn
vi lgs.fofn
vi run.cfg
NextPolish-v1.0.2/nextPolish run.cfg
```

- Hi-C assembly

```shell
samtools faidx ${sample}.fasta
chromap -i -r ${sample}.fasta -o contigs.index

chromap \
    --preset hic \
    -r ${sample}.fasta \
    -x contigs.index \
    --remove-pcr-duplicates \
    -1 ${sample}_R1.fq.gz \
    -2 ${sample}_R2.fq.gz \
    --SAM \
    -o ${sample}.sam \
    -t ${thread}

samtools view -bh -u -F0xF0C -q ${thread} ${sample}.bam | bedtools bamtobed | \
  awk -v OFS='\t' '{$4=substr($4,1,length($4)-2); print}' > ${sample}.bed
yahs ${sample}.fasta ${sample}.bed

juicer pre -a -o out_JBAT \
    yahs.out.bin \
    yahs.out_scaffolds_final.agp \
    ${sample}.fasta.fai

JUICER=juicer_tools_1.19.02.jar
asm_size=$(awk '{s+=$2} END{print s}' ${sample}.fasta.fai)
java -Xmx36G -jar $JUICER \
    pre out_JBAT.txt out_JBAT.hic assembly ${asm_size}
```

- ALLMAPs assembly

```shell
grep "gene" ${sample}.gff3 |awk -F ";" '{print $1}'|awk '{print $1"\t"$4"\t"$5"\t"$9}'|sed 's/ID=//g' >${sample}.bed
sed 's/evm.model/evm.TU/g' ${sample}.cds >${sample}.cds
cp ref.bed ./
cp ref.cds ./
python -m jcvi.compara.catalog ortholog --nostdpf ref ${sample}
python -m jcvi.assembly.syntenypath bed ref.${sample}.anchors --switch -o ${sample}.transfer.bed
python -m jcvi.assembly.allmaps mergebed ${sample}.transfer.bed -o ${sample}.transfer.merge.bed
python -m jcvi.assembly.allmaps path ${sample}.transfer.merge.bed ${sample}.contig.fasta
liftOver -gff ${sample}.gff3 ${sample}.transfer.merge.chain ${sample}.genome.gff3 ${sample}.genome.unmapped;
```

## 2. Extract nucleus and organelles genomes

```shell
# Step 1: Align contigs to organelle reference genomes using nucmer
nucmer -t ${threads} ${IRGSP-1.0_organelles.fasta} ${sample}.contig.fasta --delta=${sample}_2organelles.delta
delta-filter -q ${sample}_2organelles.delta >${sample}_2organelles_q.delta
show-coords -clrT -I 80 -L 100 ${sample}_2organelles_q.delta >${sample}_2organelles_q.coords

# Step 2: Convert coordinates to BED format
grep -v "[\[|\/]" ${sample}_2organelles_q.coords |grep -v "NUCMER" |awk '{if ($3>$4) print $13"\t"$4"\t"$3; else print $13"\t"$3"\t"$4}' |awk NF |sort -k 1,2n  >${sample}_2organelles.bed
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' ${sample}.contig.fasta | awk '{print $1"\t"length($2)}' |sed 's/>//g' |awk '{print $1"\t1\t"$2}' >${sample}.bed

# Step 3: Calculate coverage of organelle alignments on contigs
# Criteria: Classify contigs as organelle if contig length <= 500kb and coverage >= 50%
bedtools coverage -a ${sample}_contig.bed -b ${sample}_2organelles.bed |sort -k1 >${sample}_2organelles_coverage.txt
awk '{if ($6<=500000 && $7>=0.5) print $0}' ${sample}_2organelles_coverage.txt |grep -v "[:|/]" | cut -f 1 > ${sample}_organelles_contig.list
awk '{if ($6>500000 || $7<0.5) print $0}' ${sample}_2organelles_coverage.txt |grep -v "[:|/]"  | cut -f 1 > ${sample}_nucleus_contig.list

# Step 4: Extract contig IDs for nuclear and organelle sequences
python get_fasta.py ${sample}_nucleus_contig.list ${sample}.contig.fasta ${sample}_nucleus.fasta
python get_fasta.py ${sample}_organelles_contig.list ${sample}.contig.fasta ${sample}_organelles.fasta
```

## 3. Evaluation of genome assemblies

- N50

```shell
perl get.chlength_N50.pl ${sample}_nucleus.fasta >${sample}.size.txt
```

- LAI

```shell
perl LTR_FINDER_parallel -seq ${sample}_nucleus.fasta -harvest_out -threads ${threads}
gt suffixerator -db ${sample}_nucleus.fasta -indexname ${sample} -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index ${sample} -maxlenltr 7000 >${sample}.harvest.scn
cat ${sample}.harvest.scn ${sample}.finder.combine.scn >${sample}_rawLTR_7k.scn
perl LTR_retriever -genome ${sample}_nucleus.fasta -inharvest ${sample}_rawLTR_7k.scn -threads ${threads}
```

- BUSCO

```shell
busco -m geno -i ${sample}_nucleus.fasta -l embryophyta_odb10 -o ${sample} --out_path ${path}/${sample} -c ${threads} --offline -f
```

- QUAST

```shell
quast.py ${sample}.fa -R ${reference} -g gene:${reference}.gff3 -l ${sample} -o result
```

- QV

```shell
inspector.py -c ${sample}.fasta -r ${sample}.fastq -o ${path} --datatype hifi -t ${thread}
```

- syntenic analysis at the genome level

```shell
nucmer -p ${sample}_Nip -t ${thread} ${reference} ${sample}.fasta
delta-filter -i 85 -l 1000 -o ${sample}_Nip.delta -1 >${sample}_Nip.filter_85_1000.delta
mummerplot -p ${sample}_Nip.filtered ${sample}_Nip.filter.delta --color --filter --fat --png
```

## 4. The organelle genomes de novo assembly

```shell
minimap2 -t ${threads} -x map-hifi ${IRGSP-1.0_chloroplast.fasta} ${sample}.hifi.fastq >${sample}.paf
perl paf_alignment_parameter_update.pl ${sample}.paf
awk '{if($5>0.7) print $1}' ${sample}.paf.aln >${sample}.paf.0.7.list
seqtk-1.3/seqtk subseq ${sample}.hifi.fastq ${sample}.paf.0.7.list >${sample}.chl.0.7.fastq
hifiasm-0.16.0/hifiasm -t {threads} -o ${sample}.asm ${sample}.chl.0.7.fastq
awk '/^S/{print ">"$2;print $3}' ${sample}.asm.bp.p_ctg.gfa >${sample}.p_ctg.fa

nucmer --prefix=${sample}.chl ${IRGSP-1.0_chloroplast.fasta} ${sample}.p_ctg.fa
show-coords -rcl ${sample}.chl.delta >${sample}.chl.coords
delta-filter -q ${sample}.chl.delta >${sample}.chl.filter
show-coords -rcl ${sample}.chl.filter >${sample}.chl.filter.coords
show-tiling -i 99 -v 50 -V 0 ${sample}.chl.filter >${sample}.chl.filter.tiling
show-tiling -a -i 99 -v 50 -V 0 ${sample}.chl.filter >${sample}.chl.filter.a.tiling
```
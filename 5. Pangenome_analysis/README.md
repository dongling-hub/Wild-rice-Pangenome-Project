# Pan-genome analysis

- [Pan-genome analysis](#pan-genome-analysis)
  - [1. Construct gene-based pan-genomes](#1-construct-gene-based-pan-genomes)
  - [2. Construct graph-based pan-genomes](#2-construct-graph-based-pan-genomes)
  - [3. Calling variations using graph-based pan-genomes](#3-calling-variations-using-graph-based-pan-genomes)

## 1. Construct gene-based pan-genomes

```shell
# Merge all CDS.fa files
cat *.CDS.fa > all.CDS.fa

# Perform BLAST search for each sample CDS.fa file against the merged CDS.fa file
formatdb -i all.CDS.fa -p F
blastall -i ${sample}.CDS.fa -d all.CDS.fa -p blastn -e 1e-10 -m 9 -b 1000 -a ${thread} -o ${sample}.out -F F
cat *.out > all.CDS.blastn.out # Combine all BLAST output files

# Parse the combined BLAST output to a table
perl get_length.pl all.CDS.fa # Get the length of each sequence
perl unique_merge.pl all.CDS.blastn.out ${list} all.CDS.fa.len 95

# Identify ortholog genes and remove redundant clusters
blast_del_v3 -s 0.5 all.CDS.blastn.95.merge all.CDS.blastn.95.merge.del
```

## 2. Construct graph-based pan-genomes

We used the [Minigraph-Cactus Pangenome Pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) with Cactus v2.2.1 to construct the assembled genomes of 16 *O. sativa*, 129 *O. rufipogon*, and a combined set of 145 accessions, respectively.

We need to create an input seqfile `${PREFIX}.seqfile`, firstly.

```shell
head -n4 ${PREFIX}.seqfile
IRGSP /data/reference/IRGSP-1.0_all.fasta
W0103 /data/assembly/W0103.contig.fasta
W0107 /data/assembly/W0107.contig.fasta
W0147 /data/assembly/W0147.contig.fasta
```

- 1.Creating the initial minigraph graph

```shell
cactus-minigraph ${jobstore} ${PREFIX}.seqfile.txt ${PREFIX}.gfa.gz \
    --mapCores ${thread} --realTimeLogging --reference IRGSP >nohup_cactus_minigraph.out
```

- 2.Mapping the Genomes Back to the Minigraph

```shell
cactus-graphmap ${jobstore} ${PREFIX}.seqfile.txt ${PREFIX}.gfa.gz ./${PREFIX}.paf \
    --outputFasta ${PREFIX}.gfa.fa.gz --reference IRGSP --delFilter 5000000 \
    --mapCores ${thread} --realTimeLogging >nohup_graphmap.out
```

- 3.Splitting by the chromosome

```shell
cactus-graphmap-split ${jobstore} ${PREFIX}.seqfile.txt ${PREFIX}.gfa.gz ${PREFIX}.paf \
    --reference IRGSP --outDir ./${PREFIX}.split >nohup_${PREFIX}.split.out
```

- 4.Creating the cactus alignment

```shell
cactus-align-batch ${jobstore} ./${PREFIX}.split/chromfile.txt ./${PREFIX}.hal \
    --alignCores ${thread} --alignOptions "--pangenome --maxLen 10000 --reference IRGSP --outVG" \
    >nohup_batch.out
```

- 5.Making the full graph

```shell
cactus-graphmap-join ${jobstore} --vg $(for j in `cut -f 1 ./${PREFIX}.split/chromfile.txt`; do echo ./${PREFIX}.hal/${j}.vg; done) \
    --hal $(for j in `cut -f 1 ./${PREFIX}.split/chromfile.txt`; do echo ./${PREFIX}.hal/${j}.hal; done) \
    --outDir ./${PREFIX}.full --outName ${PREFIX}.full --reference IRGSP --gfaffix --wlineSep "." \
    --indexCores ${thread} >nohup_join_full.out
```

- 6.Making the default graph

```shell
cactus-graphmap-join ${jobstore} --vg $(for j in `cut -f 1 ./${PREFIX}.split/chromfile.txt`; do echo ./${PREFIX}.full/clip-${PREFIX}.full/${j}.vg; done) \
    --outDir ./${PREFIX}.default --outName ${PREFIX}.default --reference IRGSP --wlineSep "." \
    --clipLength 10000 --clipNonMinigraph --vcf --giraffe --preserveIDs --indexCores ${thread} \
    --realTimeLogging >nohup_join_default.out
```

- 7.Filtering by Allele Frequency

```shell
cactus-graphmap-join ${jobstore} --vg $(for j in `cut -f 1 ./${PREFIX}.split/chromfile.txt`; do echo ./${PREFIX}.default/clip-${PREFIX}.default/${j}.vg; done) \
    --outDir ./${PREFIX}.filter --outName ${PREFIX}.filter --reference IRGSP --wlineSep "." \
    --vgClipOpts "-d 12 -m 1000" --preserveIDs --giraffe --indexCores ${thread} \
    --realTimeLogging >nohup_join_filter.out
```

- 8.Splitting the graph into components

```shell
vg chunk -x ./${PREFIX}.filter/${PREFIX}.filter.xg -M -O pg -t ${thread}
```

## 3. Calling variations using graph-based pan-genomes

- 1.mapping

We used Giraffe to map the Illumina short paired-end reads from each accession against the graph-based cultivated-wild pangenome.

```shell
vg giraffe -x ${PREFIX}.filter.xg -H ${PREFIX}.filter.gbwt \
    -m ${PREFIX}.filter.min -d ${PREFIX}.filter.dist \
    -N ${sample} -f ${sample}.R1.fastq.gz -f ${sample}.R2.fastq.gz \
    -P -o BAM -t ${thread} >${sample}.bam
samtools sort ${sample}.bam -@ ${thread} -m 2G -O BAM -o ${sample}.sort.bam
samtools index -@ ${thread} ${sample}.sort.bam
```

- 2.Calling variations

The variations were then called using DeepVariants with the NGS model, and all individual variants were merged using GLnexus.

```shell
# Call variations using DeepVariant
export SINGULARITYENV_CUDA_VISIBLE_DEVICES=${gpu_id}
INPUT_DIR="${path}/${sample}"
OUTPUT_DIR="${path}/${sample}"
singularity run --nv \
    -B "${INPUT_DIR}":"/input","${OUTPUT_DIR}":"/output" \
    ${path}/deepvariant_1.6.1-gpu.sif /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS --ref=${reference} --reads=/input/${sample}.sort.bam \
    --sample_name ${sample} --output_vcf=/output/${sample}.vcf.gz \
    --output_gvcf=/output/${sample}.g.vcf.gz \
    --make_examples_extra_args="min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true" \
    --intermediate_results_dir /output/${sample}.tmp \
    --num_shards=${thread}

# Merge variants using GLnexus
glnexus_cli --dir ${path} \
    --config DeepVariant --list ${merge}.list \
    --threads ${thread} > ${merge}.bcf
bcftools view ${merge}.bcf --threads ${thread} | bgzip -@ ${thread} -c > ${merge}.vcf.gz
tabix -p vcf ${merge}.vcf.gz
```
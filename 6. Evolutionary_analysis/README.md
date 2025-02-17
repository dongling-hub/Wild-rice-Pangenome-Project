# Evolutionary analysis

- [Evolutionary analysis](#evolutionary-analysis)
  - [1. Phylogenetic tree construction](#1-phylogenetic-tree-construction)
  - [2. Archetypal analysis](#2-archetypal-analysis)
  - [3. Calculating genetic index](#3-calculating-genetic-index)
  - [4. Demographic history inference](#4-demographic-history-inference)
  - [5. Introgression analysis](#5-introgression-analysis)

## 1. Phylogenetic tree construction

- single-copy gene phylogenetic tree

We used `DIAMOND` for protein comparison and `OrthoFinder` for orthogroups. A phylogenetic tree was built with `IQ-TREE` using single-copy orthologues.

```shell
# Run OrthoFinder with DIAMOND for protein comparison
orthofinder -S diamond -t ${thread} -og -f ${path} -b ${path}/OrthoFinder/Results_Nov18_1/WorkingDirectory/

# Build a phylogenetic tree with IQ-TREE using single-copy orthologues
iqtree2 -s ${SpeciesTreeAlignment.fasta} --seqtype AA --prefix ${prefix} --seed 12345 -T AUTO -B 1000 -m MFP
```

- SNP-based phylogenetic tree

```shell
# 1) Convert VCF file to TPED file
plink --noweb --vcf ${merge}.vcf.gz --recode 12 --output-missing-genotype 0 --transpose --biallelic-only --out ${merge} --threads ${thread}
# 2) Replace map file
awk 'BEGIN{OFS=""}{print $1,"\tchr",$1,"_",$4,"\t",$3,"\t",$4}' ${merge}.map >${merge}.map.tmp && mv ${merge}.map.tmp ${merge}.map
# 3) Calculate hIBS matrix (EMMAx-kin):
emmax-kin -v -h -s -d 10 ${merge}
# 4) Convert hIBS matrix to distance matrix for neighbor.exe
perl change.kinshipToNeighorInfile.pl ${merge}.hIBS.kinf ${merge}.infile ${merge}.tfam
# 5) Construct phylogenetic tree with neighbor.exe
```

## 2. Archetypal analysis

We identified core SNPs subsets and conducted archetypal analysis via the `archetypal-analysis`.

```shell
# Filter VCF file by minor allele frequency and missing data
vcftools --gzvcf ${merge}.vcf.gz --maf 0.05 --max-missing 0.8 --recode --out ${merge}.filter

# Replace 'chr0' with 'chr' in the filtered VCF file
sed -i 's/chr0/chr/g'  ${merge}.filter.recode.vcf

# Convert VCF to PLINK format
vcftools --vcf ${merge}.filter.recode.vcf --plink --out ${merge}.filter.recode

# Perform linkage disequilibrium pruning
plink --file ${merge}.filter.recode --indep-pairwise 100 10 0.5 --out ${merge}.filter.recode --noweb

# Extract pruned SNPs and convert to PED format
plink --file ${merge}.filter.recode --extract ${merge}.filter.recode.prune.in --recode --out ${merge}.filter.recode.prunein --noweb &

# Convert PED to BED format
plink --noweb --file ${merge}.filter.recode.prunein --make-bed --out ${merge}.filter.recode.prunein

# Convert BED to VCF format
plink2 --bfile ${merge}.filter.recode.prunein --export vcf --out ${merge}.filter.recode.prunein

# Perform archetypal analysis for k values from 2 to 10
for i in `seq 2 10`; do
  archetypal-analysis -i ${merge}.filter.recode.prunein.vcf -o ${prefix} -k $i --tolerance 0.0001 --max_iter 400
  
  # Generate bar plot for archetypal analysis results
  archetypal-plot -i ${prefix}.$i.Q -p bar_simple -s ${supplementary}.txt -title ${prefix}
  
  # Generate simplex plot for archetypal analysis results
  archetypal-plot -i ${prefix}.$i.Q -p plot_simplex -s ${supplementary}.txt -title ${prefix}
done
```

## 3. Calculating genetic index

- Fst

```shell
vcftools --gzvcf ${merge}.vcf.gz --weir-fst-pop ${group1}.list --weir-fst-pop ${group2}.list --fst-window-size 100000 --fst-window-step 10000 --out ${group1_group2}
```

- pi

```shell
vcftools --gzvcf ${merge}.vcf.gz --keep ${group}.list --window-pi 100000 --window-pi-step 10000 --out ${group}
```

- LD decay

```shell
PopLDdecay -InVCF ${merge}.vcf.gz -SubPop ${group}.list -OutStat ${group}.stat
```

- DST

```shell
plink --file ${merge}.vcf.gz --genome --genome-full
```

## 4. Demographic history inference

```shell
sh random_group.sh ${group1}.list ${group2}.list 50 >${group1_group2}.combine.list
for i in {1..50}
do
  sed -n "${i}p" ${group1_group2}.combine.list | cut -f 2- | tr "\t" "\n" | awk NF > group$i.list
  sh generate.pseudodiploid_input.run_msmc.sh group$i.list group$i
done
```

## 5. Introgression analysis

- TreeMix

```shell
sh vcf2treemix.sh ${merge}.vcf.gz ${population}.txt
for m in {0..10}; do
  for i in {1..10}; do
    # Generate random seed
    s=$RANDOM
    # Generate bootstrapped input file with ~80% of the SNP loci
    gunzip -c ${merge}.treemix.frq.gz | awk 'BEGIN {srand()} { if (NR==1) {print $0} else if (rand() <= .8) print $0}' | gzip > ${merge}.bs${i}.m${m}.treemix.frq.gz
    # Run treemix on bootstrapped input file
    treemix -i ${merge}.bs${i}.m${m}.treemix.frq.gz -o ${merge}.bs${i}.m${m} -root CC,O.meridionalis,O.longistaminata -global -m ${m} -k 500 -seed ${s} > treemix.bs${i}.m${m}_log 2> nohup_treemix.bs${i}.m${m}.out
  done
done
```

- f3-admixture

```shell
sh convertVCFtoEigenstrat.sh ${merge}.vcf.gz
cat f3.par
genotypename:   ${merge}.eigenstratgeno
snpname:        ${merge}.snp
indivname:      sample.group.list
popfilename:    f3.pop
inbreed: YES

cat f3.pop
${pop1} ${pop2} ${pop3}

head -n3 sample.group.list
02428 U japonica/East_Asia
534M  U indica/East_Asia
9311  U indica/East_Asia

qp3Pop -p f3.par >f3.out
```

- ABBA-BABA test

```shell
python ~/biosoft/genomics_general/VCF_processing/parseVCF.py -i ${merge}.vcf.gz | gzip > ${merge}.geno.gz

python ~/biosoft/genomics_general/freq.py -g ${merge}.geno.gz -p ${spc1} -p ${spc2} -p ${spc3} -p ${spc4} \
  --popsFile ${group}.list --target derived -t ${thread} | grep -v nan | gzip > ${spc1}_${spc2}_${spc3}_${spc4}.tsv.gz

Rscript ~/biosoft/genomics_general/calculate_abba_baba.r ${spc1}_${spc2}_${spc3}_${spc4}.tsv.gz \
  ${spc1}_${spc2}_${spc3}_${spc4}.abba_baba.txt ${spc1} ${spc2} ${spc3} ${spc4} chr_lengths.txt

python ~/biosoft/genomics_general/ABBABABAwindows.py -g ${merge}.geno.gz -o ${spc1}_${spc2}_${spc3}_${spc4}.abba_baba.windows.txt \
  -P1 ${spc1} -P2 ${spc2} -P3 ${spc3} -O ${spc4} --popsFile ${group}.list -f phased -w 100000 -s 10000 \
  --windType coordinate --minSites 250 --minData 0.4 --Threads ${thread}

Rscript ~/biosoft/genomics_general/plot_abbababa_windows.r ${spc1}_${spc2}_${spc3}_${spc4}.abba_baba.windows.txt \
  chr_lengths.txt ${spc1}_${spc2}_${spc3}_${spc4}.abba_baba.windows.pdf
```
#!/usr/bin/env bash

# This is the adapted version of panEDTA
#      Individual TE libraries will be combined by panEDTA to generate a pan-genome TE library.
#      The pan-genome TE library will be used to annotate TEs in all genomes.

## Test files explained
# Col.test.fa   Genome 1
# Ler.test.fa   Genome 2
# genome.cds.fa CDS sequence. Could be absent, or from all input genomes, or from one of the genomes, or from a closely related species.
#               If users provided only one CDS file, then use it for all genomes. Each genome could use a different CDS.
#               This file is optional, but having it will improve the annotation by purging gene sequences from TEs.
# athrep.ref.formatted  Existing TE library from this species. Optional. Providing a curated library will significantly improve
#               the annotation and also harvest existing knowledge.


# Help info
helpFunction()
{
   echo ""
   echo "Usage: $0 -genomes genome_list.txt -cds cds.fasta -threads 10"
   echo -e "\t-g        A list of genome files with paths accessible from the working directory.
                        Required: You can provide only a list of genomes in this file (one column, one genome each row).
                        Optional: You can also provide both genomes and CDS files in this file (two columns, one genome and one CDS each row).
                                  Missing of CDS files (eg, for some or all genomes) is allowed."
   echo -e "\t-c                Optional. Coding sequence files in fasta format.
                        The CDS file provided via this parameter will fill in the missing CDS files in the genome list.
                        If no CDS files are provided in the genome list, then this CDS file will be used on all genomes."
   echo -e "\t-l        Optional. A manually curated, non-redundant library following the RepeatMasker naming format."
   echo -e "\t-f        Minimum number of full-length TE copies in individual genomes to be kept as candidate TEs for the pangenome.
                        Lower is more inclusive, and will ↑ library size, ↑ sensitivity, and ↑ inconsistency.
                        Higher is more stringent, and will ↓ library size, ↓ sensitivity, and ↓ inconsistency.
                        Default: 3."
   echo -e "\t-t        Number of CPUs to run panEDTA. Default: 10."
   echo ""
   exit 1 # Exit script after printing help
}

# Preset parameters
genome_list=''
cds=''
curated_lib=''
fl_copy=3
threads=10

# Read user inputs
while getopts "g:c::l::f::t::" opt
do
   case "$opt" in
      g ) genome_list="$OPTARG" ;;
      c ) cds="$OPTARG" ;;
      l ) curated_lib="$OPTARG" ;;
      f ) fl_copy="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Check parameters
if [ ! -s "$genome_list" ]; then
   echo "ERROR: The genomes $genomes_list file is not found or is empty";
   helpFunction
fi

if [ "$cds" != '' ] & [ ! -s "$cds" ]; then
   echo "ERROR: The cds $cds file is not found or is empty"
   helpFunction
fi

if [ "$curated_lib" != '' ] & [ ! -s "$curated_lib" ]; then
   echo "ERROR: The curated library $curated_lib file is not found or is empty"
   helpFunction
fi

# Set paths
path=$(dirname "$0") #program path
dir=$(pwd) #current work directory


### Begin panEDTA and print all parameters
echo `date`
echo "Pan-genome Extensive de-novo TE Annotator (panEDTA)"
echo -e "\tGenome files: $genome_list"
echo -e "\tCoding sequences: $cds"
echo -e "\tCurated library: $curated_lib"
echo -e "\tCopy cutoff: $fl_copy"
echo -e "\tCPUs: $threads"

## Step 2, make pan-genome lib (quick step, use a single node is fine)
# get fl-TE with ≥ $fl_copy copies in each genome
for i in `cat $genome_list`; do   
        genome=`echo $i|awk '{print $1}'`
        accession=`echo $i|awk '{print $3}'`
        echo "Identify full-length TEs for genome $genome"
        perl $path/util/find_flTE.pl $accession/$genome.mod.EDTA.anno/$genome.mod.EDTA.RM.out | \
                awk '{print $10}'| \
                sort| \
                uniq -c |\
                perl -snale 'my ($count, $id) = (split); next if $count < $fl_copy; print $_' -- -fl_copy=$fl_copy | \
                awk '{print $2"#"}' > $accession/$genome.mod.EDTA.TElib.fa.keep.list
done

# extract pan-TE library candidate sequences
for i in `cat $genome_list`; do
        genome=`echo $i|awk '{print $1}'`
        accession=`echo $i|awk '{print $3}'`
        if [ -s "$curated_lib" ]; then

                # a) if --curated_lib is provided
                for j in `cat $accession/$genome.mod.EDTA.TElib.fa.keep.list`; do
                        grep $j $accession/$genome.mod.EDTA.TElib.novel.fa; 
                done | \
                        perl $path/util/output_by_list.pl 1 $accession/$genome.mod.EDTA.TElib.novel.fa 1 - -FA > $accession/$genome.mod.EDTA.TElib.fa.keep.ori
        else
                # b) if --curated_lib is not provided
                for j in `cat $accession/$genome.mod.EDTA.TElib.fa.keep.list`; do 
                        grep $j $accession/$genome.mod.EDTA.TElib.fa; 
                done | \
                        perl $path/util/output_by_list.pl 1 $accession/$genome.mod.EDTA.TElib.fa 1 - -FA > $accession/$genome.mod.EDTA.TElib.fa.keep.ori
        fi
done

# aggregate TE libs
i=0
for j in ./*/*keep.ori; do 
        i=$(($i+5000)); 
        perl $path/util/rename_TE.pl $j $i; 
done | perl $path/util/rename_TE.pl - > panEDTA.TElib.fa.raw

# remove redundant
echo "Generate the panEDTA library"
perl $path/util/cleanup_nested.pl -in panEDTA.TElib.fa.raw -cov 0.95 -minlen 80 -miniden 80 -t $threads
cp panEDTA.TElib.fa.raw.cln panEDTA.TElib.fa

# Extra step if --curated_lib is provided:
if [ -s "$curatedlib" ]; then
        cat $curated_lib >> panEDTA.TElib.fa
fi
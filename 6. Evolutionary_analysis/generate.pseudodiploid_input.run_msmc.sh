#!/bin/bash

### Variables:
if [ $# -lt 2 ]; then
    echo "Usage: $0 <species_list_file> <output_prefix>"
    exit 1
fi

IND=$(cat "$1") # Pass the species name list
POP=$2 # Output group name prefix
OUTDIR=$(pwd)
path="your_path_here" # Define the path variable

for s in {01..12}
do
    echo "working on chromosome $s"
    chr=chr$s
    MSMC_INPUT=${OUTDIR}/${POP}.$chr.multihetsep.txt
    MASK_GENOME=${path}/IRGSP-1.0_genome_${chr}.mask.bed.gz
    negative_mask=${path}/MSU_${chr}.bed.gz

    printf "\n \n \n \n"
    date
    echo "Script: msmc_2_generateInput_multiInd"
    echo "Individuals: ${IND}"
    echo "Population: $POP"
    echo "Chromosome: ${chr}"
    echo "MSMC input file: ${MSMC_INPUT}"

    mask_file=$(mktemp)
    vcf_file=$(mktemp)

    for ind in $IND
    do
        INDMASK=$(ls ~/pan_genome_of_rice/evolution/MSMC/$ind/$ind.mask.$chr.bed.gz)
        echo "--mask=$INDMASK " >> $mask_file
        INDVCF=$(ls ~/pan_genome_of_rice/evolution/MSMC/$ind/$ind.$chr.phased.vcf.gz)
        echo $INDVCF >> $vcf_file
    done

    ### Generate MSMC input files:
    generate_multihetsep.py $(tr "\n" " " < $mask_file) --negative_mask=${negative_mask} --mask=$MASK_GENOME $(tr "\n" " " < $vcf_file) > ${MSMC_INPUT}
    python ~/biosoft/msmc-tools-master/msmc2_scripts/generate_pseudodiploid.py ${MSMC_INPUT} ${OUTDIR}/${POP}.$chr.pseudodiploid.txt
    echo "Done with ${chr}; moving on to next chromosome"

    rm $mask_file $vcf_file
done

echo "Done with generate pseudodiploid input"

### Run MSMC
msmc2_Linux -t 10 -I 0,1,2,3 -o pop1 ${OUTDIR}/${POP}.*.pseudodiploid.txt
msmc2_Linux -t 10 -I 4,5,6,7 -o pop2 ${OUTDIR}/${POP}.*.pseudodiploid.txt
msmc2_Linux -t 10 -s -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 -o pop1_pop2 ${OUTDIR}/${POP}.*.pseudodiploid.txt
combineCrossCoal.py ${OUTDIR}/pop1_pop2.final.txt ${OUTDIR}/pop1.final.txt ${OUTDIR}/pop2.final.txt > ${OUTDIR}/${POP}.combined.final.txt
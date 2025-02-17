#!/bin/bash

# Set the sample name (Modify this variable as needed)
sample="your_sample_name"

# Step 1: Read heterozygous alleles from the allele list
echo "Reading heterozygous alleles..."
hetro_p_gene_list=($(awk '{print $2}' ${sample}.hetrozygous_allele.list))  # Get primary alleles
hetro_a_gene_list=($(awk '{print $1}' ${sample}.hetrozygous_allele.list))  # Get alternative alleles
hetro_num=${#hetro_p_gene_list[@]}  # Count the number of heterozygous alleles

# Check if the allele list is empty
if [[ $hetro_num -eq 0 ]]; then
    echo "Error: No heterozygous alleles found. Exiting..."
    exit 1
fi

# Step 2: Process each tissue type
for tissue in leaf seedling root fringe; do
    echo "Processing tissue: $tissue..."
    output_tmp="${sample}.${tissue}_hetro_transcript_expression_tmp"
    output="${sample}.${tissue}_hetro_transcript_expression.txt"

    # Check if the expression quantification file exists for this tissue
    quant_file="${sample}.all.${tissue}_quant/quant.sf"
    if [[ ! -f $quant_file ]]; then
        echo "Error: Missing quantification file for $tissue: $quant_file"
        continue  # Skip this tissue and move to the next one
    fi

    # Step 3: Extract expression values for heterozygous genes
    for ((i = 0; i < hetro_num; i++)); do
        gene_p="${hetro_p_gene_list[$i]}"
        gene_a="${hetro_a_gene_list[$i]}"

        # Extract TPM and Read Counts for each allele
        expression_p=$(grep -w "$gene_p" $quant_file)
        expression_a=$(grep -w "$gene_a" $quant_file)

        # Assign default values (0) if gene is not found
        p_TPM=$(echo "$expression_p" | awk '{print ($4 != "") ? $4 : 0}')
        p_ReadCounts=$(echo "$expression_p" | awk '{print ($5 != "") ? $5 : 0}')
        a_TPM=$(echo "$expression_a" | awk '{print ($4 != "") ? $4 : 0}')
        a_ReadCounts=$(echo "$expression_a" | awk '{print ($5 != "") ? $5 : 0}')

        # Write extracted values to temporary output file
        echo -e "$gene_p\t$p_TPM\t$p_ReadCounts\t$gene_a\t$a_TPM\t$a_ReadCounts" >> $output_tmp
    done

    # Step 4: Calculate log2FC and Allele-Specific Expression Ratio (ADER)
    echo "Calculating log2FC and ADER for $tissue..."
    awk '{
        if (NF == 6) print $0;  # Ensure only complete records are processed
    }' $output_tmp | awk '{
        if ($3 != 0 && $6 != 0) 
            print $0 "\t" (log($2 / $5) / log(2)) "\t" $3 / ($3 + $6) "\t" $6 / ($3 + $6);
        else if ($3 == 0 && $6 != 0) 
            print $0 "\t" (log($2 + 0.0000001 / $5) / log(2)) "\t0\t1";
        else if ($3 != 0 && $6 == 0) 
            print $0 "\t" (log($2 / ($5 + 0.0000001)) / log(2)) "\t1\t0";
        else 
            print $0 "\t-\t-\t-";
    }' > ${output}

    # Step 5: Filter ASE (Allele-Specific Expression) transcripts based on thresholds
    # Filtering Criteria:
    # - Sum of reads from both alleles (p_ReadCounts + a_ReadCounts) >= 10
    # - Absolute log2FC >= 1 (Fold Change >= 2)
    # - At least one allele has ADER (Allele-Derived Expression Ratio) >= 0.75
    awk '{
        if (($3 + $6 >= 10) && ($7 >= 1 || $7 <= -1) && ($8 >= 0.75 || $9 >= 0.75)) print $0;
    }' ${output} > ${sample}.${tissue}_hetro_ASE.txt
done

echo "All tissues processed successfully!"
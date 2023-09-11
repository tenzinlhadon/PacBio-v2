#!/bin/bash 
fasta_dir="/home/ubuntu/analysis_docker/analysis/unmasked_references"
fasta_files=("$fasta_dir"/*)
sample="m64409e_230203_155352.hifi_reads_0.0001_sampled"
echo "reference,aligned_reads" > "$sample"_primary_aligned_counts.temp.csv
for fasta_file in "${fasta_files[@]}"; do 
    reference="$(basename "$fasta_file" | cut -d '.' -f 1)"
    sequence_ids=($(grep '^>' "$fasta_file" | sed 's/^>//g' | tr '\t' ' ' |  cut -d ' ' -f 1 ))
    input_bam="m64409e_230203_155352.hifi_reads_0.0001_sampled.sorted.bam"
    output_bam="$reference"".bam"
    # echo $output_bam
    # echo "${sequence_ids[@]}"
    if [[ "$fasta_file" != "combined_ref.fa" ]]; then 
        samtools view -h -b -o $output_bam $input_bam "${sequence_ids[@]}"
    else
        output_bam="$input_bam"
    fi
           
    samtools view -h -b  -F 0x104 -o  "${sample}_${reference}_cov.bam"  "$output_bam"
    bedtools genomecov -ibam "${sample}_${reference}_cov.bam" -dz -split > "${sample}_${reference}.cov.bed"
    awk -F'\t' 'BEGIN { OFS = FS} { $2 = $2 + 1} 1' "${sample}_${reference}.cov.bed" > temp_cov  &&  mv temp_cov "${sample}_${reference}.cov.bed"
    bedtools genomecov -ibam "${sample}_${reference}_cov.bam" -bg > "${sample}_${reference}_coverage.txt"
    
    output_stats="${sample}${reference}.stats"
    echo "Metrics,${reference}" > "$output_stats"
    fasta_length=$(python -c "from Bio import SeqIO; total_length = sum(len(record.seq) for record in SeqIO.parse('${fasta_file}', 'fasta')); print(total_length)")
    echo $fasta_length
    echo "Read length:,${fasta_length}" >> "$output_stats"
    samtools stats $output_bam  -@ 12  | grep ^SN | cut -f 2- | awk -F'\t' -v OFS="," '{print $1,$2}' >> "$output_stats"
    awk 'BEGIN { FS=OFS="," } { gsub(":", "", $1) } 1' "$output_stats" > temp_stats && mv temp_stats  "$output_stats" 
    grep -f /home/ubuntu/docker/import/metrics.txt "$output_stats" > temp_stats &&  mv temp_stats "$output_stats"
    grep -v "readsmappedandpaired" "$output_stats" > tempfile && mv tempfile "$output_stats"
    ref_read_count=$(samtools view -c -F 0x104 "$output_bam")
    echo "${reference},${ref_read_count}" >> "${sample}_primary_aligned_counts.temp.csv"
done


# python3 percentage_sort.py "${sample}_primary_aligned_counts.temp.csv" "$read_count"
python3 percentage_sort.py "${sample}_primary_aligned_counts.temp.csv" 345678900
rm *cov.bed
python3 alignment_metrics.py "$sample"
rm *stats
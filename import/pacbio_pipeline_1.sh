#!/bin/bash 
## PacBio pipeline by Elucidata, 2022

# read arguments from JSON file


five_prime_ITR=`jq '.five_prime_ITR_region' /analysis/job_file.json`
three_prime_ITR=`jq '.three_prime_ITR_region' /analysis/job_file.json`
clustering_identity_threshold=`jq '.clustering_identity_threshold' /analysis/job_file.json`
subsampling_rate=`jq '.subsampling_rate' /analysis/job_file.json`
intact_cassette_lower_threshold=`jq '.intact_cassette_lower_threshold' /analysis/job_file.json` 
intact_cassette_upper_threshold=`jq '.intact_cassette_upper_threshold' /analysis/job_file.json`
minimum_chimeric_read_support=`jq '.minimum_chimeric_read_support' /analysis/job_file.json`
s3_output_folder=`jq '.output_folder' /analysis/job_file.json`
max_insertion_length_sniffles=`jq '.max_insertion_length_sniffles' /analysis/job_file.json`
max_splits_per_read_sniffles=`jq '.max_splits_per_read_sniffles' /analysis/job_file.json`
min_alignment_length_sniffles=`jq '.min_alignment_length_sniffles' /analysis/job_file.json`
usearch_cluster_sequence_type=`jq '.usearch_cluster_sequence_type' /analysis/job_file.json`


shopt -s nullglob


rm /analysis/unmasked_references/combined_ref.fa
rm /analysis/features_BED_files/combined_annotation.temp.bed 
rm /analysis/features_BED_files/combined_annotation.bed



#reference files to assign the chr to each ref and to caluclate the total length of combined references

array3=(/analysis/unmasked_references/*)
awk 1 /analysis/unmasked_references/* > /analysis/unmasked_references/combined_ref.fa
awk 1 /analysis/features_BED_files/*bed > /analysis/features_BED_files/combined_annotation.temp.bed
python /import/annotation.py 

list_samples=$(ls -d /analysis/sample_files/*/ | grep -v '/reference_masking_output/')


#####feature list#########

for i in $list_samples;do

    cd $i

    sample=""
    sample=$(echo *.bam | awk -F'.bam' '{print $1}')

    echo " ######################################################"
    echo " ##### processing sample $sample "
    echo " ######################################################"

    echo " ### converting BAM to FASTQ  ### "
    bam2fastq -u -o $sample "$sample".bam 

    echo "checking mem and diskspace before gunzip"
    free -h
    df -h 
    echo " ####### extracting fastq.zip file ##### "
    
    # gunzip "$sample".fastq.gz
    sample=$(echo *.fastq | awk -F'.fastq' '{print $1}')
    rm "$sample".bam
    
    /tools/bbmap/reformat.sh in1="$sample".fastq out1="$sample"_"$subsampling_rate"_sampled.fastq  samplerate=$subsampling_rate
    rm "$sample".fastq
    [ ! -f "$sample"_1_sampled.fastq ] || mv "$sample"_1_sampled.fastq "$sample".fastq

    sample=$(echo *.fastq | awk -F'.fastq' '{print $1}')


    python3 /import/read_length.py "$intact_cassette_lower_threshold" "$intact_cassette_upper_threshold" "$sample.fastq"



    mkdir FASTQC_output
    echo "###### Performing FASTQC on the reads  #####"
    fastqc  -o FASTQC_output "$sample".fastq


    echo " ###### Checking DNA Contamination ##################################################################################"

    echo "######## aligning reads to  masked references #######" 

    /tools/bwa.kit/bwa mem -t 24  -x pacbio /analysis/reference_masking_output/bwa_index/masked_ref.idx "$sample".fastq > "$sample".sam
    
    #function to convert SAM file into a sorted and indexed BAM file  
    sam_2_sortedBAI(){
     base_name=$(basename $1 .sam)
     /tools/sambamba-0.8.2 view -S -f  bam -o "$base_name.bam" "$1"  
     /tools/sambamba-0.8.2 sort --tmpdir=$(pwd) -m 5G -t 9 -o "$base_name.sorted.bam" "$base_name.bam"
     samtools index  -@ 20 "$base_name.sorted.bam" "$base_name.sorted.bai" 
     return 0
    }
    
    #function to convert BAM file into a sorted and indexed BAM file 
    bam_2_sortedBAI(){
     base_name=$(basename $1 .bam)
     /tools/sambamba-0.8.2 sort --tmpdir=$(pwd) -m 5G -t 9 -o "$base_name.sorted.bam" "$base_name.bam" 
     samtools index -@ 20 "$base_name.sorted.bam" "$base_name.sorted.bai"
     return 0
    }

    sam_2_sortedBAI "$sample".sam
    samtools view -H -o "$sample".hdr "$sample".sam
    rm  "$sample".sam

     
    #extract primary aligned reads from the raw BAM file 
    samtools view  -@ 24 -h  -b -F 0x904 -o  "$sample"_primary.bam "$sample".sorted.bam 
    bam_2_sortedBAI "$sample"_primary.bam
    read_count=$(samtools view -c  "$sample"_primary.sorted.bam)
    samtools view  -@ 24 -h -b  -o "$sample"_primary_cassette.bam "$sample"_primary.sorted.bam cassette
    bedtools bamtobed -i "$sample"_primary_cassette.bam > "$sample"_cassette.bed
    python /import/cassette_breakpoint.py "$sample"_cassette.bed  "$read_count"

    
    mkdir Cassette_breakpoint 
    mv *cassette_combined_read_counts.csv ./Cassette_breakpoint/
    mv *cassettebreakpoint_read_count.png ./Cassette_breakpoint/
    cur_dir=$(pwd) 
    sample_name=$(basename $cur_dir)

    
    python /import/upload_output_files.py $sample_name "Cassette_breakpoint" 


    rm -rf Cassette_breakpoint
    
    #extract primary supplementary aligned reads from the raw BAM file 
    samtools view -@ 24  -h -F 0x104 -o  "$sample"_primary_sup.sam "$sample".sorted.bam 
    sam_2_sortedBAI "$sample"_primary_sup.sam 
    rm "$sample"_primary_sup.sam 
    prim_sup_read_count=$(samtools view -c "$sample"_primary_sup.sorted.bam)
    
    #alignment metrics and aligned reads to the references 

    
    fasta_dir="/analysis/unmasked_references"
    fasta_files=("$fasta_dir"/*)
    echo "reference,aligned_reads" > "$sample"_primary_aligned_counts.temp.csv
    echo "reference,aligned_reads" > "$sample"_primary_sup_aligned_counts.temp.csv

    for fasta_file in "${fasta_files[@]}"; do 
        reference="$(basename "$fasta_file" | cut -d '.' -f 1)"
        sequence_ids=($(grep '^>' "$fasta_file"  | sed -E 's/^>[[:space:]]*([^[:space:]]+).*$/\1/' | cut -d ' ' -f 1))
        input_bam="$sample.sorted.bam"
        output_bam="${reference}.bam"
        if [[ "$fasta_file" != "combined_ref.fa" ]]; then
            samtools view  -@ 24  -h -b -o $output_bam $input_bam "${sequence_ids[@]}"
        else
            output_bam="$input_bam"
        fi
           
        samtools view  -@ 24 -h -b  -F 0x104 -o  "${sample}_${reference}_cov.bam"  "$output_bam"
        bedtools genomecov -ibam "${sample}_${reference}_cov.bam" -dz -split > "${sample}_${reference}.cov.bed"
        awk -F'\t' 'BEGIN { OFS = FS} { $2 = $2 + 1} 1' "${sample}_${reference}.cov.bed" > temp_cov  &&  mv temp_cov "${sample}_${reference}.cov.bed"
        bedtools genomecov -ibam "${sample}_${reference}_cov.bam" -bg > "${sample}_${reference}_coverage.txt"
        
        output_stats="${sample}${reference}.stats"
        echo "Metrics,${reference}" > "$output_stats"
        fasta_length=$(python -c "from Bio import SeqIO; total_length = sum(len(record.seq) for record in SeqIO.parse('${fasta_file}', 'fasta')); print(total_length)")
        echo $fasta_length
        echo "Read length:,${fasta_length}" >> "$output_stats"
        samtools stats $output_bam  -@ 20  | grep ^SN | cut -f 2- | awk -F'\t' -v OFS="," '{print $1,$2}' >> "$output_stats"
        awk 'BEGIN { FS=OFS="," } { gsub(":", "", $1) } 1' "$output_stats" > temp_stats && mv temp_stats  "$output_stats" 
        grep -f /import/metrics.txt "$output_stats" > temp_stats &&  mv temp_stats "$output_stats"
        grep -v "readsmappedandpaired" "$output_stats" > tempfile && mv tempfile "$output_stats"
        ref_read_count=$(samtools view -c -F 0x904 "$output_bam")
        ref_prim_sup_read_count=$(samtools view -c -F 0x104 "$output_bam")
        echo "${reference},${ref_read_count}" >> "${sample}_primary_aligned_counts.temp.csv"
        echo "${reference},${ref_prim_sup_read_count}" >> "${sample}_primary_sup_aligned_counts.temp.csv"
       
    done

    find . -type f -empty  -delete

    python3 /import/percentage_sort.py "${sample}" "_primary_aligned"  "$read_count"
    python3 /import/percentage_sort.py  "$sample" "_primary_sup_aligned" "$prim_sup_read_count"

    # Define the input parameters
    sample="$sample"
    python3 /import/alignment_metrics.py "$sample"
    # rm "$sample"_primary_sup.sorted.bam
    rm *stats

    #getting coverage plots for each reference
    Rscript /import/coverage_plot.r "$sample"_

    rm *.cov.bed 

    echo " ###### Checking chimeric species ##################################################################################"

    #identifying chimeric alignments
    samtools view  -@ 24  -h "$sample".sorted.bam | grep -e '^@' -e 'SA:Z'  > "$sample"_temp.sam
    samtools view -b "$sample"_temp.sam > "$sample"_SA.bam
    rm "$sample"_temp.sam 
    bam_2_sortedBAI "$sample"_SA.bam

    #Convert BAM to BED to get the coordinates of the chimeric alignments

    bedtools bamtobed -i "$sample"_SA.sorted.bam -cigar > "$sample"_supplementary.bed
    rm "$sample"_SA.bam


    # Look for chimeric species only if supplementary alignments are present 

    if [[ -s "$sample"_supplementary.bed ]]; then
        echo "Chimeric species found, identifying constructs"
        python3 /import/chimeric_species.py "$sample"_supplementary.bed "$read_count" "$minimum_chimeric_read_support"
    else
        echo " ##############   No chimeric species found     ####################"
    fi 

    # Sniffles structural variant calling
    mkdir Sniffles_output 

    sniffles -i "$sample"_primary_sup.sorted.bam  -v Sniffles_output/"$sample"_without_qc_sniffles.vcf  --no-qc 
    sniffles -i "$sample"_primary_sup.sorted.bam  -v Sniffles_output/"$sample"_with_qc_sniffles.vcf --reference /analysis/unmasked_references/combined_ref.fa --minsupport $minimum_chimeric_read_support --long-ins-length $max_insertion_length_sniffles --max-splits-base $max_splits_per_read_sniffles -t 12  --min-alignment-length $min_alignment_length_sniffles
 
    rm "$sample"_primary_sup.sorted.bam


    python /import/upload_output_files.py $sample_name "Sniffles_output" 



    mkdir -p ./Contamination_detection/{alignment_metrics,BAM_files,coverage_plots,coverage_txt}
    mv "$sample"_primary.bam  "$sample".sorted.bam "$sample".sorted.bai   *cov.bam  ./Contamination_detection/BAM_files/
    mv *_alignment_per_ref_metrics.csv *_reference_percentage.csv ./Contamination_detection/alignment_metrics/
    mv *_coverage.png ./Contamination_detection/coverage_plots/
    mv *coverage.txt ./Contamination_detection/coverage_txt/   



    ###Chimeric species######
    mkdir Chimeric_species
    mv "$sample"_SA.sorted.bam  ./Chimeric_species/
    mv "$sample"_SA.sorted.bai ./Chimeric_species/
    [ ! -f *annotated_chimeric_species.csv ] || mv *annotated_chimeric_species.csv *annotated_chimeric_species_coordinates.csv ./Chimeric_species/



    python /import/upload_output_files.py $sample_name "Contamination_detection" 

    python /import/upload_output_files.py $sample_name  "Chimeric_species" 

    find . -mindepth 1 ! -name "*.fastq" -exec rm -rf {} +

    mkdir -p Cassette_heterogeneity/{Consensus_configs,read_length_histogram_plots,read_length_csv_files}
    mkdir -p ITR_heterogeneity/{ITR_sequence,ITR_consensus_sequence,Consensus_config_alignment_csv,Consensus_config_alignment_txt}



    echo " ###### Checking ITR heterogeneity ##################################################################################"

    sed -i '/^>/s/[^a-zA-Z>0-9  ]/_/g' /analysis/cassette_configurations.fa
    /tools/bwa.kit/bwa index /analysis/cassette_configurations.fa

    #align reads to the wild cassette configurations

    echo " ########### Aligning reads to cassette cofigurations ##############"

    /tools/bwa.kit/bwa mem -t 24  -x pacbio /analysis/cassette_configurations.fa "$sample".fastq > "$sample"_heterogeneity.sam
    samtools view -h -b -o "$sample"_heterogeneity.bam "$sample"_heterogeneity.sam
    rm "$sample"_heterogeneity.sam
    samtools view -H "$sample"_heterogeneity.bam | grep "@SQ" | cut -f 2 -d ":" | cut -f1 > ref_heterogeneity.txt

    samtools view -@ 24 -h -b -F 0x904 -o "$sample"_heterogeneity_primary.bam "$sample"_heterogeneity.bam  
    cassette_heterogeneity_read_count=$(samtools view -c "$sample"_heterogeneity_primary.bam)
    bam_2_sortedBAI "$sample"_heterogeneity_primary.bam
    rm "$sample"_heterogeneity.bam
    readarray -t list_ref_heterogeneity < ref_heterogeneity.txt

    
    # Cassette configuration wise separation 

    while read ref_seq; do

        echo "split sam reference wise"
        samtools view -@ 24 "$sample"_heterogeneity_primary.sorted.bam "$ref_seq" | awk '{print length($10)}' > "$sample"_"$ref_seq"_read_length.txt


    done < ref_heterogeneity.txt


    awk '{if(NR%4==2) print length($1)}' "$sample".fastq > "$sample"_all_reads_read_length.txt 
    rm "$sample".fastq 
    intact_cassette_lower_threshold=$intact_cassette_lower_threshold
    intact_cassette_upper_threshold=$intact_cassette_upper_threshold
    partial_genome="length(seq)<$intact_cassette_lower_threshold"
    complete_genome="length(seq)>=$intact_cassette_lower_threshold&&length(seq)<=$intact_cassette_upper_threshold"

    #sorting reads into partial and complete based on the threshold provided
    samtools view -e $complete_genome -O BAM -o "$sample"_heterogeneity.intact.bam "$sample"_heterogeneity_primary.sorted.bam  
    samtools view -e $partial_genome -O BAM -o "$sample"_heterogeneity.partial.bam "$sample"_heterogeneity_primary.sorted.bam

    # rm "$sample"_heterogeneity_primary.sam    
    bam_2_sortedBAI "$sample"_heterogeneity.intact.bam
    rm "$sample"_heterogeneity.intact.bam
    bam_2_sortedBAI "$sample"_heterogeneity.partial.bam
    rm "$sample"_heterogeneity.partial.bam


    #delete any empty file generated
    find . -type f -empty  -delete

    ####read length histogram #########
    echo " ## generating histogram plots ##"
    python /import/read_length_histogram.py $sample

    #extracting ITRs from the reads
    Rscript /import/extract_ITR.r  $sample $three_prime_ITR $five_prime_ITR
    sed -i '/^[^>]/s/[^ATGCNatgc]//g' *.fa
    
    #usearch clustering of the ITRs 

    echo " ####### Clustering ITRs ########"
    for FILE in ./*.fa
    do
    prefix=$(basename $FILE .fa)
    /tools/usearch -sortbylength "$prefix".fa -fastaout "$prefix".sorted.fa
    if [ "$usearch_cluster_sequence_type" = "consensus" ]; then
        echo "### generating consensus for ITRs  #####" 
        /tools/usearch  -cluster_fast "$prefix".sorted.fa -id $clustering_identity_threshold -consout "$prefix"_consensus.fasta -sizeout -threads 24
    else 
        /tools/usearch  -cluster_fast "$prefix".sorted.fa -id $clustering_identity_threshold -centroids "$prefix"_consensus.fasta -sizeout -threads 24
        echo "### generating centroids for ITRs  #####" 

    fi

    rm "$prefix".sorted.fa
    # /tools/usearch -sort length -cluster_fast "$prefix".fa -id "$clustering_identity_threshold"  -consout"$prefix"_consensus.fasta -sizeout -threads 24 
    done

    #aligning ITR consensus configurations to wild ITRs

    for i in  ./*_consensus.fasta
    do
    prefix=$(basename $i _consensus.fasta)
    echo "$prefix"
    /tools/blat -t=dna -q=dna -out=blast8 /analysis/ITR.fa  "$prefix"_consensus.fasta "$prefix"_consensus_aln.temp.txt
    /tools/blat -t=dna -q=dna -out=blast /analysis/ITR.fa  "$prefix"_consensus.fasta "$prefix"_consensus_aln.txt

    done

    #calculating percentage of reads present in each ITR cluster and assigning header

    for file in  ./*_consensus_aln.temp.txt
    do 
    
    prefix=$(basename $file _consensus_aln.temp.txt)
    sed -i "s/;size=/`echo "\t"`/g;s/;//g" "$prefix"_consensus_aln.temp.txt
    done

    echo "checking mem and diskspace before genome clustering"
    free -h
    df -h 

    echo " ###### Checking cassette heterogeneity ##################################################################################"
    
    /tools/bbmap/reformat.sh  in="$sample"_heterogeneity_primary.bam out="$sample"_heterogeneity.fa

    /tools/usearch -sortbylength "$sample"_heterogeneity.fa -fastaout "$sample"_heterogeneity.sorted.fa

    if [ "$usearch_cluster_sequence_type" = "consensus" ]; then
        echo " #### generate consensus sequences for cassettes #### "
        /tools/usearch -cluster_fast "$sample"_heterogeneity.sorted.fa -id $clustering_identity_threshold  -consout "$sample"_cassette_consensus.fasta -sizeout -threads 24
    else
        echo " #### generate centroid sequences for cassettes #### "
        /tools/usearch -cluster_fast "$sample"_heterogeneity.sorted.fa -id $clustering_identity_threshold  -centroids "$sample"_cassette_consensus.fasta -sizeout -threads 24
    fi

    rm "$sample"_heterogeneity.fa "$sample"_heterogeneity.sorted.fa
    awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' "$sample"_cassette_consensus.fasta > temp && mv temp "$sample"_cassette_consensus.fasta
    /tools/blat -t=dna -q=dna -out=blast8 /analysis/cassette_configurations.fa "$sample"_cassette_consensus.fasta  "$sample"_heterogeneity_consensus_aln.temp.txt
    /tools/blat -t=dna -q=dna -out=blast /analysis/cassette_configurations.fa "$sample"_cassette_consensus.fasta  "$sample"_heterogeneity_consensus_aln.txt
    sed -i "s/;size=/`echo "\t"`/g;s/;//g" "$sample"_heterogeneity_consensus_aln.temp.txt
 
    python /import/cluster_sort.py  "$cassette_heterogeneity_read_count" 

    echo "######### Analysis complete, Now moving files #################"  


   ###genome heterogeneity######
    
    mv "$sample"_heterogeneity_primary.bam   ./Cassette_heterogeneity/
    mv "$sample"_cassette_consensus.fasta ./Cassette_heterogeneity/Consensus_configs/
    mv "$sample"_heterogeneity_consensus_aln.csv ./Cassette_heterogeneity/Consensus_configs/
    mv "$sample"_heterogeneity_consensus_aln.txt ./Cassette_heterogeneity/Consensus_configs/
    mv *histogram.png ./Cassette_heterogeneity/read_length_histogram_plots/
    mv *read_length_range.csv ./Cassette_heterogeneity/read_length_csv_files/
    
    ###ITR####
    mv *aln.txt ./ITR_heterogeneity/Consensus_config_alignment_txt/
    mv *_aln.csv ./ITR_heterogeneity/Consensus_config_alignment_csv/
    mv *prime.fa ./ITR_heterogeneity/ITR_sequence/  
    mv *prime_consensus.fasta ./ITR_heterogeneity/ITR_consensus_sequence/
    

    
    find . -maxdepth 1 -type f -delete 
    python /import/upload_output_files.py $sample_name "ITR_heterogeneity" 
    python /import/upload_output_files.py $sample_name "Cassette_heterogeneity" 

    rm -rf * 

done
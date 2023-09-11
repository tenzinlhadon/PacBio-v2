#!/bin/bash 
# Pacbio pipeline by Elucidata, 2022

#read arguments from JSON file
#lscpu
five_prime_ITR=`jq '.five_prime_ITR_region' /analysis/job_file.json`
three_prime_ITR=`jq '.three_prime_ITR_region' /analysis/job_file.json`
intact_genome_threshold=`jq '.intact_genome_threshold' /analysis/job_file.json`
clustering_identity_threshold=`jq '.clustering_identity_threshold' /analysis/job_file.json`
echo $clustering_identity_threshold
intact_lower_threshold=`jq '.intact_lower_threshold' /analysis/job_file.json` 
intact_upper_threshold=`jq '.intact_upper_threshold' /analysis/job_file.json`


subsampling_rate=`jq '.subsampling_rate' /analysis/job_file.json`
echo $subsampling_rate
# shopt -s nullglob
rm /analysis/unmasked_references/combined_ref.fa
rm /analysis/features_BED_files/combined_annotation.temp.bed
rm /analysis/features_BED_files/combined_annotation.bed

#reference files to assign the chr to each ref and to caluclate the total length of combined references
array3=(/analysis/unmasked_references/*)
awk 1 /analysis/unmasked_references/* > /analysis/unmasked_references/combined_ref.fa
awk 1 /analysis/features_BED_files/* > /analysis/features_BED_files/combined_annotation.temp.bed
python /import/annotation.py 

# list_samples=$(ls -d /analysis/sample_files/*/)
list_samples=$(ls -d /analysis/sample_files/*/ | grep -v '/reference_masking_output/')




mkdir /analysis/streamlined_output 



####feature list#########

for i in $list_samples;do

    cd $i
    echo $i
    # rm *sample*
    rm *sample*
    rm *txt*

    # sample=""
    # sample=$(echo *.bam | awk -F'.bam' '{print $1}')
    # echo " ######################################################"
    # echo " ##### processing sample $sample ######################"
    # echo " ######################################################"
    # echo " ### converting BAM to FASTQ  ### "
    # bam2fastq -o $sample "$sample".bam 
    # echo "checking mem and diskspace before gunzip"
    # free -h
    # df -h 
    # echo " ####### extracting fastq.zip file ##### "
    # gunzip "$sample".fastq.gz
    # echo "fastq list"
    # echo *.fastq
    sample=$(echo *.fastq | awk -F'.fastq' '{print $1}')
    ls
    pwd
    # sample=$(ls *.fastq 2>/dev/null | head -n 1 | sed 's/\.fastq$//')
    echo $sample

    
    /tools/bbmap/reformat.sh in1="$sample".fastq out1="$sample"_"$subsampling_rate"_sampled.fastq  samplerate=$subsampling_rate
    #rm "$sample".fastq
    [ ! -f "$sample"_1_sampled.fastq ] || mv "$sample"_1_sampled.fastq "$sample".fastq

    python3 /import/read_length.py "$intact_lower_threshold" "$intact_upper_threshold"

    sample=$(echo *.fastq | awk -F'.fastq' '{print $1}')
    ref_list=()
    length_all_ref=$(grep -v ">" /analysis/unmasked_references/combined_ref.fa| tr -d '\n' | wc -c)
    ref_list+=('Ref_length_(bp)')
    ref_list+=($length_all_ref)
    mkdir FASTQC_output
    echo "###### Performic FASTQC on the reads  #####"
    fastqc  -o FASTQC_output "$sample".fastq
    echo "$sample" >  "$sample"_alignment_metrics_per_ref.stats

    echo " ###### Checking DNA Contamination ##################################################################################"

    echo "######## aligning reads to  masked references #######"

    /tools/bwa.kit/bwa mem -t 15 -x pacbio /analysis/reference_masking_output/bwa_index/masked_ref.idx "$sample".fastq > "$sample".sam
    #samtools view -h input.bam | sed -e 's/pHELP-KanV4.SGR.(S.04MAY19.N.4Dec19)/pHELP-KanV4.SGR._S.04MAY19.N.4Dec19_/g' -e 's/pNLRep2-Caprh74kan.SGR.(S.30MAY19.N.4Dec19)/pNLRep2-Caprh74kan.SGR._S.30MAY19.N.4Dec19_/g' | samtools view -b -o output.bam 
    
    #function to convert SAM file into a sorted and indexed BAM file  
    sam_2_sortedBAI(){
     base_name=$(basename $1 .sam)
     samtools view -bh "$1" | \
     samtools sort -l 9 -@ 10 -o "$base_name.sorted.bam"
     samtools index -@ 10 "$base_name.sorted.bam" "$base_name.sorted.bai"
     return 0
    }

    #function to convert BAM file into a sorted and indexed BAM file 
    bam_2_sortedBAI(){
     base_name=$(basename $1 .bam)
     samtools sort -l 9 -@ 12 -o "$base_name.sorted.bam" "$base_name.bam"
     samtools index -@ 3 "$base_name.sorted.bam" "$base_name.sorted.bai"
     return 0
    }

    
    sam_2_sortedBAI "$sample".sam
    samtools view -H -o "$sample".hdr "$sample".sam

    #extract primary aligned reads from the raw BAM file  
    
    samtools view -h -F 0x904 -o  "$sample"_primary.sam "$sample".sorted.bam 
    sam_2_sortedBAI "$sample"_primary.sam 
    
    #getting alignment metrics from the BAM file 
    samtools stats -@ 12 "$sample".sorted.bam | grep ^SN | cut -f 2- | awk -F'\t' '{print $1,$2}' > "$sample".stats 
    awk -F'\t' '{print $1,$2}' "$sample".stats > "$sample"_merged.stats
    echo -e "Metrics,$sample" | cat - "$sample".stats > "$sample"_merged.stats 
    str_file="$sample"_merged.stats"   "
    echo "reference,primary_aligned_reads" > "$sample"_primary_aligned_counts.temp.csv
    bedtools genomecov -ibam "$sample".sorted.bam -dz -split > "$sample"_all_references.cov.bed
    bedtools genomecov -ibam "$sample".sorted.bam -bg > "$sample"_coverage.txt
    echo "######## bed file generated#######"
    
    echo " reference wise coverage and alignment metrics generation" 
    #finding reads aligning to each reference 
    for i in "${array3[@]}";do

        ref=$(basename $i | sed -r 's|^(.*?)\.\w+$|\1|')
        #sed -i '/^>/s/[^a-zA-Z>0-9  ]/_/g' $i
        grep '^>' $i | sed 's/>//g' | awk '{print $1; }' > "$sample"_"$ref".txt
        readarray -t list_ref < "$sample"_"$ref".txt
        for j in "${list_ref[@]}";do
            awk '{if($3=="'"$j"'")print}' "$sample".sam > "$sample"_"$j".sam  
            cat "$sample"_"$j".sam >> "$sample"_"$ref".temp.sam
        done

        cat "$sample".hdr "$sample"_"$ref".temp.sam > "$sample"_"$ref"_ref.sam 
        sam_2_sortedBAI "$sample"_"$ref"_ref.sam
        bedtools genomecov -ibam "$sample"_"$ref"_ref.sorted.bam -dz -split > "$sample"_"$ref".cov.bed
        bedtools genomecov -ibam "$sample"_"$ref"_ref.sorted.bam -bg > "$sample"_"$ref"_coverage.txt
        primary_reads=$(samtools view -c -F 0x904 "$sample"_"$ref"_ref.sorted.bam)
        echo "$ref,$primary_reads" >> "$sample"_primary_aligned_counts.temp.csv
        samtools stats -@ 12 "$sample"_"$ref"_ref.sorted.bam | grep ^SN | cut -f 2- | awk -F'\t' '{print $2}' > "$sample"_"$ref".temp.stats
        echo -e "$ref" | cat - "$sample"_"$ref".temp.stats > "$sample"_"$ref".stats
        str_file+="$sample"_"$ref".stats" "
        paste -d , $str_file > "$sample"_alignment_metrics_per_ref.temp.stats 
        count_ref="$(grep -v ">" $i| tr -d '\n' | wc -c)"
        ref_list+=($count_ref)

        echo  "$count_ref" >> ref_lengths.txt 

    done
        
    primary_count="$sample"_primary_aligned_counts.temp.csv
    echo $str_file 
    grep -f /import/metrics.txt "$sample"_alignment_metrics_per_ref.temp.stats > "$sample"_alignment_metrics_per_ref.temp.2.stats
    cp "$sample"_alignment_metrics_per_ref.temp.2.stats "$sample"_alignment_metrics_per_ref.temp.check.stats
    grep -v "readsmappedandpaired" "$sample"_alignment_metrics_per_ref.temp.2.stats > tempfile && cp tempfile "$sample"_alignment_metrics_per_ref.temp.2.stats
    sed '2,$s/:/,/;s/ //g;s/total length/totalbasessequenced/g' "$sample"_alignment_metrics_per_ref.temp.2.stats > "$sample"_alignment_metrics_per_ref.stats
    find . -type f -empty  -delete
    echo ${ref_list[@]} | tr ' ' , >> "$sample"_alignment_metrics_per_ref.stats

    #calulating percentage for each reference and shifting read length row to the top
    
    python3 /import/percentage_sort.py "$primary_count"

    #getting coverage plots for each reference
    Rscript /import/coverage_plot.r "$sample"_




    echo " ###### Checking chimeric species ##################################################################################"

    #identifying chimeric alignments
    samtools view -h "$sample".sorted.bam | grep -e '^@' -e 'SA:Z'  > "$sample"_temp.sam
    samtools view -b "$sample"_temp.sam > "$sample"_SA.bam
    bam_2_sortedBAI "$sample"_SA.bam

    #Convert BAM to BED to get the coordinates of the chimeric alignments

    bedtools bamtobed -i "$sample"_SA.bam  > sup_"$sample"_temp.bed
    sort sup_"$sample"_temp.bed | uniq > sup_"$sample".bed
    cut -f 1-4 sup_"$sample".bed > sup_"$sample"_temp_2.bed

    #Annotating chimeric species 
    bedtools intersect -a sup_"$sample"_temp_2.bed -b /analysis/features_BED_files/combined_annotation.bed -wa -wb > "$sample"_annotation.txt
    read_count=$(samtools view -c -F 0x904 "$sample".sorted.bam )
    find . -type f -empty  -delete

    #creating a csv with chimeric species and their read support
    python /import/chimera_coordinate.py $read_count "$sample"
    [ ! -f "$sample"_chimeric_species.csv ] || sed -i '/^$/d;s/[[:blank:]]//g' "$sample"_chimeric_species.csv 

    echo " ###### Checking ITR heterogeneity ##################################################################################"

    sed -i '/^>/s/[^a-zA-Z>0-9  ]/_/g' /analysis/cassette_configurations.fa
    /tools/bwa.kit/bwa index /analysis/cassette_configurations.fa

    #align reads to the wild cassette configurations
    /tools/bwa.kit/bwa mem -t 32 -x pacbio /analysis/cassette_configurations.fa "$sample".fastq > "$sample"_heterogeneity.sam

    sam_2_sortedBAI "$sample"_heterogeneity.sam
    
    samtools view -h -b -F 0x904 -o "$sample"_heterogeneity_primary.bam "$sample"_heterogeneity.sorted.bam  
    samtools sort -@ 10 -o "$sample"_heterogeneity_primary.sorted.bam "$sample"_heterogeneity_primary.bam
    samtools index -@ 10 "$sample"_heterogeneity_primary.sorted.bam
    samtools idxstats "$sample"_heterogeneity_primary.sorted.bam | cut -f 1 | sed 's/*//g' > ref_heterogeneity.txt
    sed -i '/^[[:space:]]*$/d' ref_heterogeneity.txt
    readarray -t list_ref_heterogeneity < ref_heterogeneity.txt

    #function to extract primary aligned reads to each cassette configuration

    primary_ref_wise_read_length(){
         samtools view -h -b -F 0x904 "$1"_"$2"_heterogeneity.sorted.bam > "$1"_"$2"_primary_ref_wise_heterogeneity.temp.bam
         samtools view "$1"_"$2"_primary_ref_wise_heterogeneity.temp.bam | awk '{print length($10)}' > "$1"_"$2"_read_length.txt        
    }
    
    # Cassette configuration wise separation 

    for i in "${list_ref_heterogeneity[@]}"
    do  
        echo $i
        samtools sort -@ 16 -o "$sample"_heterogeneity_primary.sorted.bam "$sample"_heterogeneity_primary.bam
        samtools index -@ 16 "$sample"_heterogeneity_primary.sorted.bam
        samtools view -h  "$sample"_heterogeneity_primary.sorted.bam | awk '{if($3=="'"$i"'")print}' > "$sample"_"$i"_heterogeneity.temp.sam
        samtools view -H -o "$sample"_heterogeneity.hdr "$sample"_heterogeneity_primary.sorted.bam
        cat "$sample"_heterogeneity.hdr "$sample"_"$i"_heterogeneity.temp.sam > "$sample"_"$i"_heterogeneity.sam 
        sam_2_sortedBAI "$sample"_"$i"_heterogeneity.sam 
        rm *sam
        primary_ref_wise_read_length "$sample" "$i"
 
    done

    #find read length distribution of raw reads 

    awk '{if(NR%4==2) print length($1)}' "$sample".fastq > "$sample"_All_reads_read_length.txt 
    # partial_genome="length(seq)>$intact_genome_threshold"
    # complete_genome="length(seq)<=$intact_genome_threshold"

    # Extracting read length information from the FASTQ file
    awk '{if(NR%4==2) print length($1)}' "$sample".fastq > "$sample"_All_reads_read_length.txt
    intact_lower_threshold=$intact_lower_threshold
    intact_upper_threshold=$intact_upper_threshold
    partial_genome="length(seq) < $intact_lower_threshold"
    complete_genome="length(seq) >= $intact_lower_threshold && length(seq) <= $intact_upper_threshold"
    
    #sorting reads into partial and complete based on the threshold provided
    samtools view -e $partial_genome -O BAM -o "$sample"_heterogeneity.intact.bam "$sample"_heterogeneity_primary.sorted.bam    
    samtools view -e $complete_genome   -O BAM -o "$sample"_heterogeneity.partial.bam "$sample"_heterogeneity_primary.sorted.bam    
    bam_2_sortedBAI "$sample"_heterogeneity.intact.bam
    bam_2_sortedBAI "$sample"_heterogeneity.partial.bam

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

    /tools/usearch -cluster_fast "$prefix".fa -id $clustering_identity_threshold -centroids "$prefix"_centroids_size.fasta -consout "$prefix"_consensus_size.fasta -sizeout -threads 12
    done

    #aligning ITR consensus configurations to wild ITRs

    for i in  ./*_consensus_size.fasta 
    do
    prefix=$(basename $i _consensus_size.fasta)

    /tools/blat -t=dna -q=dna -out=blast8 /analysis/ITR.fa  "$prefix"_consensus_size.fasta  "$prefix"_consensus_aln.temp.txt
    /tools/blat -t=dna -q=dna -out=blast /analysis/ITR.fa  "$prefix"_consensus_size.fasta  "$prefix"_consensus_aln.txt

    done

    #calculating percentage of reads present in each ITR cluster and assigning header

    for file in  ./*_consensus_aln.temp.txt
    do 
    
    prefix=$(basename $file _consensus_aln.temp.txt)
    sed -i "s/;size=/`echo "\t"`/g;s/;//g" "$prefix"_consensus_aln.temp.txt
    awk 'NR==FNR{sum+= $2; next}{printf("%0.4f\n", $2/sum  * 100)}' "$prefix"_consensus_aln.temp.txt "$prefix"_consensus_aln.temp.txt >> temp 
    paste "$prefix"_consensus_aln.temp.txt temp  > "$prefix"_consensus_size_aln.temp2.txt && mv "$prefix"_consensus_size_aln.temp2.txt "$prefix"_consensus_aln.temp.txt
    echo -e "Cluster\tsize\tITR\t%_identity\talignment_length\tmismatches\tgap_openings\tquery_start\tquery_end\tsubject_start\tsubject_end\tE_value\tbit_score\tsize_%" | cat - "$prefix"_consensus_aln.temp.txt > temp && mv temp  "$prefix"_consensus_aln.temp.txt    
    done



    echo "checking mem and diskspace before genome clustering"
    free -h
    df -h 

    echo " ###### Checking cassette heterogeneity ##################################################################################"

    /tools/bbmap/reformat.sh in="$sample"_heterogeneity_primary.sorted.bam out="$sample"_heterogeneity.fa
    /tools/usearch -cluster_fast "$sample"_heterogeneity.fa -id 0.95 -consout "$sample"_cassette_consensus.fasta -sizeout -threads 13
    awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' "$sample"_cassette_consensus.fasta > temp && mv temp "$sample"_cassette_consensus.fasta
    /tools/blat -t=dna -q=dna -out=blast8 /analysis/cassette_configurations.fa "$sample"_cassette_consensus.fasta  "$sample"_heterogeneity_consensus_aln.temp.txt
    /tools/blat -t=dna -q=dna -out=blast /analysis/cassette_configurations.fa "$sample"_cassette_consensus.fasta  "$sample"_heterogeneity_consensus_aln.txt

    sed -i "s/;size=/`echo "\t"`/g;s/;//g" "$sample"_heterogeneity_consensus_aln.temp.txt
    awk 'NR==FNR{sum+= $2; next}{printf("%0.4f\n", $2/sum  * 100)}' "$sample"_heterogeneity_consensus_aln.temp.txt "$sample"_heterogeneity_consensus_aln.temp.txt >> temp 
    paste "$sample"_heterogeneity_consensus_aln.temp.txt temp  > "$sample"_heterogeneity_consensus_size_aln.temp2.txt && mv "$sample"_heterogeneity_consensus_size_aln.temp2.txt   "$sample"_heterogeneity_consensus_aln.temp.txt
    echo -e "Cluster\tsize\tcassette\t%_identity\talignment_length\tmismatches\tgap_openings\tquery_start\tquery_end\tsubject_start\tsubject_end\tE_value\tbit_score\tsize_%" | cat -   "$sample"_heterogeneity_consensus_aln.temp.txt > temp && mv temp "$sample"_heterogeneity_consensus_aln.temp.txt  
    python /import/cluster_sort.py   

    echo "######### Analysis complete, Now moving files #################"  


    mkdir -p ./Contamination_detection/{alignment_metrics,BAM_files,coverage_plots,coverage_txt}
    mkdir ./Chimeric_species
    mkdir -p Cassette_heterogeneity/{Consensus_configs,read_length_histogram_plots,read_length_csv_files}
    mkdir -p ITR_heterogeneity/{ITR_sequence,Consensus_config_alignment_csv,Consensus_config_alignment_txt}
    ###DNA contamination#########
    
    mv "$sample"_primary.sorted.bam "$sample"_primary.sorted.bai  "$sample".sorted.bam "$sample".sorted.bai   *ref.sorted.bam *ref.sorted.bai ./Contamination_detection/BAM_files/
    mv "$sample"_alignment_metrics_per_ref.csv "$sample"_primary_aligned_counts.csv ./Contamination_detection/alignment_metrics/
    mv *_coverage.png ./Contamination_detection/coverage_plots/
    mv *coverage.txt ./Contamination_detection/coverage_txt/   


   ###genome heterogeneity######

    mv "$sample"_cassette_consensus.fasta ./Cassette_heterogeneity/Consensus_configs/"$sample"_cassette_consensus.fasta
    mv "$sample"_heterogeneity_consensus_aln.csv ./Cassette_heterogeneity/Consensus_configs/"$sample"_heterogeneity_consensus_aln.csv
    mv "$sample"_heterogeneity_consensus_aln.txt ./Cassette_heterogeneity/Consensus_configs/"$sample"_heterogeneity_consensus_aln.txt
    mv *histogram.png ./Cassette_heterogeneity/read_length_histogram_plots/
    mv *read_length_range.csv ./Cassette_heterogeneity/read_length_csv_files/
    ##ITR####
    mv *aln.txt ./ITR_heterogeneity/Consensus_config_alignment_txt/
    mv *_aln.csv ./ITR_heterogeneity/Consensus_config_alignment_csv/
    mv *prime.fa ./ITR_heterogeneity/ITR_sequence/

    ###Chimeric species######
    mv "$sample"_SA.bam  ./Chimeric_species/"$sample"_primary_supplementary_alignment.bam
    mv "$sample"_SA.sorted.bai ./Chimeric_species/"$sample"_primary_supplementary_alignment.bai
    [ ! -f "$sample"_chimeric_species.csv ] || mv "$sample"_chimeric_species.csv ./Chimeric_species/"$sample"_chimeric_species.csv

    find . -maxdepth 1 -type f ! -name '*.fastq' -delete
done
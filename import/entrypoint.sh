#!/bin/bash
# mkdir analysis
ls /analysis/
# rm -rf /analysis/*


aws s3 cp $S3_PATH /analysis/job_file.json
python3 /import/download_input_files.py


sample_folder="/analysis/sample_files"
unmasked_folder="/analysis/unmasked_references"
masked_folder="/analysis/reference_masking_output/bwa_index"
feature_BED_files="/analysis/features_BED_files"

# folders_to_check=("$sample_folder" "$unmasked_folder" "$features_BED_files")
# for folder in "${folders_to_check[@]}"; do
#     if [ ! -d "$folder" ] || [ -z "$(ls -A "$folder")" ]; then
#         echo "Folder '$folder' is either not found or empty. Please check the path in the JSON config file and ensure the files are present in the S3 bucket."
#         exit 1
#     fi
# done
check_folder() {
    local folder="$1"
    if [ ! -d "$folder" ]; then
        echo " $folder does not exist."
        exit 1
    fi
    if find "$folder" -maxdepth 0 -empty | grep -q "."; then
        echo " $folder exists but is empty."
        exit 1
    else
        echo " $folder exists and is not empty."
    fi
}


check_folder "$sample_folder"
check_folder "$unmasked_folder"
check_folder "$feature_BED_files" 





check_file(){
file_name=$1 
if [ -f "${file_name}" ];then
    if [ -s "${file_name}" ];then
        echo " ${file_name} found"
    else
        echo " ${file_name} is empty, please check the path in JSON config file and make sure the file is present in s3 bucket"
        exit 1
    fi
else
    echo "${file_name} not found, please check the path in JSON config file and make sure the file is present in s3 bucket"
    exit 1
fi
}

check_file "/analysis/cassette_configurations.fa"
check_file "/analysis/ITR.fa"

perform_ref_masking=`jq -r '.ref_masking' /analysis/job_file.json`
path_to_indexed_genome=`jq -r '.masked_references_folder_path' /analysis/job_file.json`
s3_output_folder=`jq '.output_folder' /analysis/job_file.json`

echo $perform_ref_masking
echo $path_to_indexed_genome

mkdir /analysis/ContaVect-1.0.0-Python3/
mkdir /analysis/ContaVect_hm_coordinates/
cp -r /import/ContaVect-1.0.0-Python3/*  /analysis/ContaVect-1.0.0-Python3/
cp -r /import/ContaVect_hm_coordinates/* /analysis/ContaVect_hm_coordinates/

# rm /analysis/unmasked_references/combined_ref.fa

if [[ "$perform_ref_masking" = "True" ]]
then

filenames=`ls /analysis/unmasked_references/`
for i in $filenames;
do
grep ">" /analysis/unmasked_references/$i  | sed s/">"// | sed 's/\r//'   |  awk -F  " "  '$0 = $1 ' | sed "s/$/\t$i/"  >> /analysis/seq_headers.tsv
done

# export PATH="/analysis/ContaVect-1.0.0-Python3:$PATH"

cd /import/ContaVect-1.0.0-Python3/

ls -l
mkdir /analysis/reference_masking_output/
mkdir /analysis/reference_masking_output/Homologous_sequences

python3 /import/ContaVect-1.0.0-Python3/ContaVect.py

cd /analysis/




# python3 /import/ContaVect_hm_coordinates/ContaVect.py
python3 /import/hm_coordinates_counts.py
rm -rf /analysis/reference_masking_output/results 
rm -rf /analysis/reference_masking_output/blast_db 
reference_masking_output_s3="${s3_output_folder}reference_masking_output/"

python3 /import/upload_reference_masking_output.py

else

aws s3 cp $path_to_indexed_genome /analysis/reference_masking_output/bwa_index/ --recursive  --quiet

fi

masked_folder="/analysis/reference_masking_output/bwa_index/"

if find "$masked_folder" -maxdepth 0 -empty | grep -q "."; then
    echo " Masked references folder is empty"
    exit 1
fi


time /import/pacbio_pipeline_1.sh


python3 /import/upload_multisample_files.py "Contamination_percentage_primary_aligned.csv"
python3 /import/upload_multisample_files.py "Contamination_percentage_primary_sup_aligned.csv"
python3 /import/upload_multisample_files.py "Intact_cassette_percentage.csv"


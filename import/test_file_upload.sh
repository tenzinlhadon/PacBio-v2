#!/bin/bash 
sample="sample_1"
s3_output_folder=`jq '.output_folder' /home/ubuntu/analysis_docker/analysis/job_file.json`
# upload_files_s3(){
#     echo $s3_output_folder
#     destination_path="${s3_output_folder}${sample}/$1"
#     command="aws s3 cp "$1" $destination_path --recursive"
#     echo $command 
#     return 0   
# }

# upload_files_s3 Cassette_breakpoint 
    if  [[ $s3_output_folder != */ ]]; then 
        s3_output_folder="${s3_output_folder%/}/"
        echo $s3_output_folder
    fi 

upload_files_s3() {
    echo "$s3_output_folder"
    destination_path="${s3_output_folder}${sample}/$1"
    command="aws s3 cp $1 \"$destination_path\" --recursive"
    echo $command
    return 0
}

upload_files_s3 Cassette_breakpoint
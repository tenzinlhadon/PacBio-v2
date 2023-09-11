#!/bin/bash 
s3_output_folder=`jq '.output_folder' /home/ubuntu/analysis_docker/analysis/job_file.json`
path=\"${s3_output_folder}m64409e_230127_174950/\"
echo $path
aws s3 ls $path

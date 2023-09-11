#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from os import listdir
from os.path import isfile, join
import subprocess as sb
import sys,os, shutil,time
import glob
import json



with open('/analysis/job_file.json') as f:
    data_parsed = json.load(f)

ref_masking = data_parsed["ref_masking"]

# print("Downloading Input files")

paths = {
    'unmasked_ref': {
        'src': data_parsed["unmasked_references_folder_path"],
        'dst': "/analysis/unmasked_references/"
    },
    'cassette_file': {
        'src': data_parsed["cassette_configuration_FASTA_path"],
        'dst': "/analysis/cassette_configurations.fa"
    },
    'feature_bed': {
        'src': data_parsed["features_BED_files_folder_path"],
        'dst': "/analysis/features_BED_files/"
    },
    'ITR':{
        'src': data_parsed["ITR_FASTA_path"],
        'dst':"/analysis/ITR.fa"
    }
}
# Download masked references only if ref_masking is set as True

if not ref_masking:
    masked_ref_src = data_parsed["masked_references_folder_path"]
    masked_ref_dst = "/analysis/reference_masking_output/bwa_index/"
    download_files(masked_ref_src,asked_ref_dst )

# Function to download input files

def download_files(copy_path, paste_path):
    if paste_path.endswith("/"):
        download_cmd = f'aws s3 cp "{copy_path}" {paste_path} --recursive --quiet'
    else: 
        download_cmd = f'aws s3 cp "{copy_path}" {paste_path} --quiet'
    print(download_cmd)
    sb.run(download_cmd, shell = True)
    print(f"Downloaded {copy_path} to {paste_path}")


# Iterate over the paths and print them
folder_lists = ['unmasked_ref', 'feature_bed']

for path_name, path_info in paths.items():
    if path_name in folder_lists:
        if not path_info['src'].endswith("/"):
            path_info['src'] = path_info['src'] + "/"
    print(f"downloading {path_name}")
    download_files(path_info['src'], path_info['dst'])

s3_input_path = data_parsed["sample_folder_path"]
if s3_input_path.endswith("/"):
    print("Sample folder path found")
else: 
    s3_input_path = s3_input_path + "/"
subfolder_lines = sb.check_output(f"aws s3 ls {s3_input_path}", shell=True, text=True).split()
subfolders = [line.split("/")[-2] for line in subfolder_lines if line.endswith('/')]
for subfolder in subfolders:
    subfolder_path = f"\"{s3_input_path}{subfolder}/\""
    subfolder_paste_path= f"/analysis/sample_files/{subfolder}/"
    download_cmd = f"aws s3 cp {subfolder_path} {subfolder_paste_path} --recursive --quiet --exclude '*' --include '*.bam' --include '*.pbi'"
    test_cmd = f"aws s3 ls {subfolder_path}"
    # sb.run(test_cmd, shell=True)
    print(download_cmd)
    sb.run(download_cmd, shell=True)

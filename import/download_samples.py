import os
from os import listdir
from os.path import isfile, join
import subprocess as sb
import sys,os, shutil,time
import glob
import json
f = open('/analysis/job_file.json')
# returns JSON object as
# a dictionary
data_parsed = json.load(f)
s3_input_path = data_parsed["sample_folder_path"]
subfolder_lines = sb.check_output(f"aws s3 ls {s3_input_path}", shell=True, text=True).split()
subfolders = [line.split("/")[-2] for line in subfolder_lines if line.endswith('/')]
for subfolder in subfolders:
    subfolder_path = f"{s3_input_path}{subfolder}"
    subfolder_paste_path= f"/home/sarepta/sample_files/{subfolder}"
    download_cmd = f"aws s3 cp {subfolder_path} {subfolder_paste_path} --recursive --exclude '*' --include '*.bam' --include '*.pbi'"
    print(download_cmd)
    sb.run(download_cmd, shell=True)
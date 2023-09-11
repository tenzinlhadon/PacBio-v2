#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from os import listdir
from os.path import isfile, join
import subprocess as sb
import sys,os, shutil,time
import glob
import json
# f = open('/analysis/job_file.json')
# # returns JSON object as
# # a dictionary
# data_parsed = json.load(f)

# ref_paste_path = data_parsed['output_folder']
# if not ref_paste_path.endswith("/"):
#     ref_paste_path = ref_paste_path + "/"
    
# print("exporting multisample result files to s3 bucket")
# def upload_multisample_result_files(file_name, ref_paste_path):
#     ref_copy_path = "/analysis/sample_files/" + file_name
#     ref_paste_path_dst =  ref_paste_path + file_name
#     download_cmd = "aws s3 cp %(ref_copy_path)s \"%(ref_paste_path_dst)s\""%globals()
#     print(download_cmd)
#     sb.check_output(download_cmd, shell = True)

# upload_multisample_result_files("Contamination_percentage_primary_aligned.csv", ref_paste_path)
# upload_multisample_result_files("Contamination_percentage_primary_sup_aligned.csv", ref_paste_path)
# upload_multisample_result_files("Intact_cassette_percentage.csv", ref_paste_path)
# aws s3 cp /analysis/sample_files/Contamination_percentage_primary_aligned.csv 

# f.close()

f = open('/analysis/job_file.json')
file_name = sys.argv[1]
# returns JSON object as
# a dictionary
data_parsed = json.load(f)
ref_paste_path = data_parsed['output_folder']
if not ref_paste_path.endswith("/"):
    ref_paste_path = ref_paste_path + "/"
ref_copy_path= "/analysis/sample_files/"  + file_name
ref_paste_path = ref_paste_path + file_name
print("exporting output folder to s3 bucket")
download_cmd = "aws s3 cp %(ref_copy_path)s %(ref_paste_path)s"%globals()
print(download_cmd)
sb.check_output(download_cmd, shell = True)
f.close()

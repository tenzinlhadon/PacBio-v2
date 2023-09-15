#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from os import listdir
from os.path import isfile, join
import subprocess as sb
import sys,os, shutil,time
import glob
import json


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

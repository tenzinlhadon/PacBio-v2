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
# returns JSON object as
# a dictionary
data_parsed = json.load(f)

ref_paste_path= data_parsed['output_folder']
if not ref_paste_path.endswith("/"):
    ref_paste_path = ref_paste_path + "/"

print("exporting output folder to s3 bucket")
ref_paste_path_dst = ref_paste_path + "reference_masking_output/"
download_cmd = "aws s3 cp /analysis/reference_masking_output \"%(ref_paste_path_dst)s\" --recursive"%globals()
print(download_cmd)
sb.check_output(download_cmd, shell = True)
f.close()
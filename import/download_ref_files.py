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

for ref in data_parsed["Running_params"]:
    ref_copy_path = ref['fasta']
    ref_paste_path = "/import/ref_path/" + ref['fasta'].split("/")[-1]
    ref['fasta'] = ref_paste_path
    print(ref_copy_path)
    print(ref_paste_path)
    download_cmd = "aws s3 cp %(ref_copy_path)s %(ref_paste_path)s"%globals()
    print(download_cmd)
    sb.check_output(download_cmd, shell = True)
# Closing file
f.close()
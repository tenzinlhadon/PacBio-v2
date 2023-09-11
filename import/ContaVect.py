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
ref_masking= eval(data_parsed["ref_masking"])
f.close()

# Function to Get the current
# working directory
def current_path():
    print("Current working directory before")
    print(os.getcwd())
    print()


# Driver's code
# Printing CWD before
current_path()

# Changing the CWD
os.chdir('/import/ContaVect-1.0.0-Python3/')

# Printing CWD after
current_path()

try:
# Standard library packages import
    os.chdir('/import/ContaVect-1.0.0-Python3/')
    from os import path, remove, environ # Mandatory package
    from time import time # Mandatory package
    import configparser # Mandatory package
    from sys import argv # Mandatory package
    import csv # Mandatory package
    import optparse # Mandatory package

    # Third party packages
    import pysam # Mandatory package
    import Bio # Mandatory package

# Local Package import
    os.chdir('/import/ContaVect-1.0.0-Python3/')
    from pyDNA.Utilities import mkdir, file_basename, file_name, expand_file, rm_blank, is_gziped
    # Mandatory package
    from pyDNA.Blast import Blastn # if not imported = not ref masking
    from pyDNA.RefMasker import mask # if not imported = not ref masking
    from pyDNA.FastqFT.FastqFilter import FastqFilter # if not imported = not fasta filter
    from pyDNA.FastqFT.QualityFilter import QualityFilter # if not imported = not fasta filter
    from pyDNA.FastqFT.AdapterTrimmer import AdapterTrimmer # if not imported = not fasta filter
    from pyDNA.Ssw import ssw_wrap # if not imported = not fasta filter
    from pyDNA.Bwa import Mem # Mandatory package
    from pyDNA.pySamTools import Bam, Coverage, Variant # if not imported = not requested output
    from ContaVect.Reference import Reference, Sequence # Mandatory package
    from ContaVect.Conf_file import write_conf_file #if not  imported = not creation of configuration file

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    exit()


class Main(object):
    """
    Main class of the program. In a first time the class is initialize with values parsed from the
    configuration file. Then, the pipeline of analysis is launch through the 'run' method
    """
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "ContaVect 0.3.0"
    USAGE = "Usage: %prog -c Conf.txt [-i -h]"

#~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        pass
        return Main()

#~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self):
        try:
            # Mandatory paramaters
            import json

            # Opening JSON file
            f = open('/analysis/job_file.json')

            # returns JSON object as
            # a dictionary
            data_parsed = json.load(f)



            self.outdir = "/analysis/reference_masking_output"
            if not self.outdir:
                self.outdir = "./"
            if self.outdir[-1] != "/":
                self.outdir += "/"
            self.outprefix = "masked_ref"
            if not self.outdir:
                self.outdir = "out"

            self.ref_masking = True
            self.R1 = "/import/ContaVect-1.0.0-Python3/test/dataset/S1_R1.fastq"
            self.R2 = "/import/ContaVect-1.0.0-Python3/test/dataset/S1_R2.fastq"
            self.input_qual = "fastq-sanger"
            self.quality_filtering = True
            self.adapter_trimming = False
            self.bwa_index = ""
            self.bwa_mem_opt = "-M"
            self.bwa_threads = 20
            self.bwa_index_opt = ""
            self.bwa_aligner = "/tools/bwa.kit/bwa mem"
            self.bwa_indexer = "/tools/bwa.kit/bwa index"
            self.min_mapq = 25
            self.min_size = 20
            self.unmapped_bam = True
            self.unmapped_sam = False
            self.cov_min_depth = 4
            self.var_min_depth = 500
            self.var_min_freq = 2

            # Conditional paramaters
            if self.ref_masking:
                self.blastn_opt = str(data_parsed["blastn_opt"])
                self.blastn_threads = str(data_parsed["blastn_threads"])
                self.mkblastdb_opt = ""
                self.blastn = "blastn"
                self.mkblastdb = "mkblastdb"

            # More complicated import in a list of dictionnary for references informations
            self.raw_ref_list =[]
            ##################################


            for ref in data_parsed["Running_params"]:
                ref_copy_path = ref['fasta']
                ref_paste_path = "/analysis/unmasked_references/" + ref['fasta'].split("/")[-1]
                ref['fasta'] = ref_paste_path
                ref['output'] = "bam bed"
                print(ref_copy_path)
                print(ref_paste_path)
                self.raw_ref_list.append(ref)
            # Closing file
            f.close()
            ##################################
        except configparser.NoOptionError as E:
            print (E)
            print ("An option is missing in the configuration file")
            print ("Please report to the descriptions in the configuration file\n")
            exit()
        except configparser.NoSectionError as E:
            print (E)
            print ("An section is missing in the configuration file")
            print ("Please report to the descriptions in the configuration file\n")
            exit()
        except ValueError as E:
            print (E)
            print ("One of the value in the configuration file is not in the correct format")
            print ("Please report to the descriptions in the configuration file\n")
            exit()

    def __repr__(self):
        msg = "MAIN CLASS\n"
        msg+= "\tParameters list\n"
        for i, j in list(self.__dict__.items()):
            msg+="\t{}\t{}\n".format(i, j)
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value

    def __call__(self):
        """
        Launch the complete pipeline of analyse:

        * Reference importation/parsing
        * Facultative step of reference masking to remove homologies between reference sequences
        * Facultative step of Fastq quality Filtering/ adapter trimming
        * Facultative step of reference indexing for bwa from merged references
        * Short read alignment with /tools/bwa.kit//tools/bwa.kit/bwa mem
        * Spliting of sam to attribute reads to each original references (or unmmapped)
        * Output per reference bam, sam, bedgraph, bed, covgraph, variant call
        * Output distribution table and graph
        """
        stime = time()
        self.outdir = mkdir(path.abspath(self.outdir))

        print ("\n##### PARSE REFERENCES #####\n")
        # Create CV_Reference.Reference object for each reference easily accessible through
        # Reference class methods

        if self.ref_masking or not self.bwa_index:
            self.ref_dir = mkdir(path.join(self.outdir, "references/"))
            self.index_dir = mkdir(path.join(self.outdir, "bwa_index/"))
            self._extract_ref(expand=True)
        else:
            self.ref_dir = ""
            self.index_dir = ""
            self._extract_ref(expand=False)

        # Reference Masking
        if self.ref_masking:
            print ("\n##### REFERENCE HOMOLOGIES MASKING #####\n")
            self.db_dir = mkdir(path.join(self.outdir, "blast_db/"))
            ref_list = self._iterative_masker()
            # Erase existing index value if ref masking was performed
            bwa_index = None


            # BWA alignment
            print ("\n##### READ REFERENCES AND ALIGN WITH BWA #####\n")
            # An index will be generated if no index was provided
            self.result_dir = mkdir(path.join(self.outdir, "results/"))

            self.sam = Mem.align (
                self.R1, self.R2,
                index = self.bwa_index,
                ref = Reference.allFasta(),
                align_opt = self.bwa_mem_opt,
                index_opt = self.bwa_index_opt,
                aligner = self.bwa_aligner,
                align_threads = self.bwa_threads,
                indexer = self.bwa_indexer,
                align_outdir = self.result_dir,
                index_outdir = self.index_dir,
                align_outname = self.outprefix+".sam",
                index_outname = self.outprefix+".idx")


        print ("\n##### DONE #####\n")
        print(("Total execution time = {}s".format(round(time()-stime, 2))))

    ##~~~~~~~PRIVATE METHODS~~~~~~~#

    def _extract_ref(self, expand=True):
        """
        Import and expand fasta references and associated flags in a Reference object
        expand the file if Gziped to avoid multiple compression/decompression during execution
        if require for next operations
        """
        for ref in self.raw_ref_list:
            # Expand fasta if needed
            if expand:
                ref_fasta = expand_file(infile=ref['fasta'], outdir=self.ref_dir)
            else:
                ref_fasta = ref['fasta']

            # Create a Reference object
            Ref = Reference(
                name = ref['name'],
                ref_fasta = ref_fasta,
                bam_maker = Bam.BamMaker(
                    make_bam = 'bam' in ref['output'],
                    make_sam = 'sam' in ref['output']),
                cov_maker = Coverage.CoverageMaker(
                    min_depth=self.cov_min_depth,
                    make_bedgraph = 'bedgraph' in ref['output'],
                    make_bed = 'bed' in ref['output'],
                    make_covgraph = 'covgraph' in ref['output']),
                var_maker = Variant.VariantMaker(
                    min_depth=self.var_min_depth,
                    min_freq=self.var_min_freq,
                    make_freqvar = 'variant' in ref['output']))

            ## Test if all seq in ref are longer than 3000 for compatibility with bwa
            #for seq in Ref.seq_dict.values():
                #if seq.length < 3000:
                    #import_and_pad (

            print((repr(Ref)))

    def _iterative_masker (self): #### TODO The fuction directly manipulate reference field= change that
        """
        Mask references homologies iteratively, starting by the last reference which is masked by
        all the others then to the penultimate masked by all others except the last and and so
        forth until there is only 1 reference remaining
        """
        # Iterate over index in Reference.instances staring by the last one until the 2nd one
        for i in range(Reference.countInstances()-1, 0, -1):

            # Extract subject and query_list from ref_list
            subject = Reference.Instances[i]
            query_list = Reference.Instances[0:i]
            print(("\n# PROCESSING REFERENCE {} #\n".format(subject.name)))

            # Perform a blast of query list against subject
            hit_list = Blastn.align (
                query_list = [ref.ref_fasta for ref in query_list],
                subject_fasta = subject.ref_fasta,
                align_opt = self.blastn_opt,
                num_threads = self.blastn_threads,
                db_opt = self.mkblastdb_opt,
                db_outdir = self.db_dir,
                db_outname = subject.name)

            # Masking hits in suject fasta if hits in hit_list
            subject.ref_fasta = mask (
                subject_fasta= subject.ref_fasta,
                hit_list = hit_list,
                ref_outdir = self.ref_dir,
                ref_outname = "masked_{}.fa".format(subject.name),
                compress_ouput = False)




#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#

main = Main()
main.class_init()
main.__call__()



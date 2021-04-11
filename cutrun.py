#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 19:16:03 2021

@author: niek
"""

import os
import argparse
import argcomplete
import subprocess
import multiprocessing
import yaml
import timeit
import time

#start timer
start = timeit.default_timer()


ap = argparse.ArgumentParser()
ap.add_argument("-t", "--threads", required=False, default=1,
   help="<INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
ap.add_argument("-r", "--rename", required=False, action='store_true',
   help="Rename fq files")
ap.add_argument("-f", "--fastqc", required=False, action='store_true',
   help="Perform FASTQC")
ap.add_argument("-a", "--align", choices=['hisat-se','hisat-se','bwa-pe','bwa-se'], required=False,
   help="Trim and align raw data to index.")
ap.add_argument("-d", "--deduplication", required=False, action='store_true',
   help="Perform deduplication of BAM files")
ap.add_argument("-s", "--downsample", required=False, action='store_true',
   help="Perform downsampling of BAM files")
ap.add_argument("-b", "--bigwig", required=False, action='store_true',
   help="Create BigWig files")
ap.add_argument("-q", "--qc", required=False, action='store_true',
   help="Perform QC analysis of BAM files")
ap.add_argument("-p", "--peaks", required=False, action='store_true',
   help="Call and annotate peaks")
ap.add_argument("-n", "--ngsplot", required=False, action='store_true',
   help="Generate metageneplots and heatmaps with ngs.plot")
args = ap.parse_args()


#Setting variablea for parsing to bash scripts
script_dir=os.path.abspath(os.path.dirname(__file__))
current_dir=os.getcwd()
picard=[line[0:] for line in subprocess.check_output("find $HOME -name picard.jar", shell=True).splitlines()]
picard=picard[0].decode("utf-8")

with open(script_dir+"/settings.yaml") as file: #reads yaml setting file into a dictionary
        settings=yaml.full_load(file)

genome=settings.get("Genome") #sets genome from settings.yaml

#set thread count for processing
max_threads=str(multiprocessing.cpu_count())
threads=args["threads"]
if threads == "max":
    threads=max_threads

rename=args["rename"]
rename_script=script_dir+"/rename.sh"
if rename == True:
    subprocess.run([rename_script])

fastqc=args["fastqc"]
fastqc_script=script_dir+"/fastqc.sh"
if fastqc == True:
    subprocess.run([fastqc_script,str(threads)])

align=args["align"]
align_options=["hisat2-se","hisat2-pe","bwa-se","bwa-pe"]
align_script=script_dir+"/align.sh"
if align in align_options:
    subprocess.run([align_script,str(threads),align,genome])
    
dedup=args["deduplication"]
dedup_script=script_dir+"/dedup.sh"
if dedup == True:
    subprocess.call([dedup_script,picard])

downsample=args["downsample"]
downsample_script=script_dir+"/downsample.sh"
if downsample == True:
    subprocess.run([downsample_script,str(threads)])
    
bigwig=args["bigwig"]
bigwig_script=script_dir+"/bigwig.sh"
if bigwig == True:
    subprocess.run([bigwig_script,str(threads)])
    
qc=args["qc"]
qc_script=script_dir+"/qc.sh"
if qc == True:
    subprocess.run([qc_script,str(threads)])
    
peaks=args["peaks"]
peaks_script=script_dir+"/peaks.sh"
if peaks == True:
    subprocess.run([peaks_script])
    
ngsplot=args["ngsplot"]
ngsplot_script=script_dir+"/ngsplot.sh"
if ngsplot == True:
    subprocess.run([ngsplot_script])

#print total run time
stop = timeit.default_timer()
total_time = stop - start
ty_res = time.gmtime(total_time)
res = time.strftime("%H:%M:%S",ty_res)
print('Total run time: ', res)
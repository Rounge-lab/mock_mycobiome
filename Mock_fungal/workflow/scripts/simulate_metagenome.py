#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 12:49:24 2024

Simulate metagenome fastq files based on the profile and reads

@author: ec-ekateria
"""

import pandas as pd
import random
import os
import subprocess as sp
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

# Add arguments
parser.add_argument('--profile', help='Give full path to the profile .csv file')
parser.add_argument('--simfol', help='Give full path to the folder with simulated data')
parser.add_argument('--outfol', help='Specify the name for the metagenome output folder')
parser.add_argument('--seed', help='Specify the seed for random sampling')

# Parse the arguments
args = parser.parse_args()


prof='/'.join(args.profile.split('/')[:-1])
prof_name=args.profile.split('/')[-1].replace('.csv','')
simfol=args.simfol
outfol=args.outfol
seed=args.seed

profile=pd.read_csv('/'.join([prof,prof_name+'.csv']))
random.seed(seed)
outpath='/'.join(['/'.join(simfol.split('/')[:-1]),outfol,prof_name])
for ix,row in profile.iterrows():
    file1=list(SeqIO.parse('/'.join([simfol,row.Accession+'.fna_1.fq']),'fastq'))
    file2=list(SeqIO.parse('/'.join([simfol,row.Accession+'.fna_2.fq']),'fastq'))
    if len(file1)>int(row.NumReads):
        idc=random.sample(range(0,len(file1)),int(row.NumReads))
    else:
        idc=list(range(0,len(file1)))
    seqs1=[file1[i] for i in idc]
    seqs2=[file2[i] for i in idc]
    if not os.path.exists(outpath) or not os.path.isdir(outpath):
        os.mkdir(outpath)
    with open('/'.join([outpath, row.Accession+'.fna_1.fq']), "w") as output_handle:
        SeqIO.write(seqs1, output_handle, "fastq")
    with open('/'.join([outpath, row.Accession+'.fna_2.fq']), "w") as output_handle:
        SeqIO.write(seqs2, output_handle, "fastq")

fn1=outpath+'_1.fq'
fn2=outpath+'_2.fq'
sp.run(f"cat {outpath}/*_1.fq > {fn1}", shell=True)
sp.run(f"cat {outpath}/*_2.fq > {fn2}", shell=True)



    

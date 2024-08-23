#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 11:11:21 2024

Download and format fungal genomes, create profiles for metagenome datasets
NCBI_DATASETS need to be installed as a conda env 

@author: ekateria
"""

import pandas as pd
import subprocess as sp
import os
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random

wdir='/FULL_PATH_TO/Mock_fungal/mock_profiles' # NB! CHANGE TO LOCAL PATH
file=pd.read_csv('/'.join([wdir,'unicellular_classes.tsv']),sep='\t')
file=file.rename(columns={'Tax name':'Tax_name'})

os.chdir(wdir)

global condaenvs, conda_env
condaenvs="PATH_TO_CONDA_ENVS" # NB! CHANGE TO LOCAL PATH
conda_env="ncbi_datasets" #use name of the environment where ncbi_datasets is installed

taxons=pd.DataFrame()
for ix,f in file.iterrows():
    output=f.Tax_name+'.zip'
    command = f"{condaenvs}/{conda_env}/bin/datasets download taxonomy taxon {f.Taxid} --children --filename {output}"
    sp.run(command, shell=True)
    sp.run(f"unzip {output} -d {f.Tax_name}", shell=True)
    
    tax=pd.read_csv('/'.join([f.Tax_name,'ncbi_dataset','data','taxonomy_summary.tsv']),sep='\t')
    taxons=pd.concat([taxons,tax],ignore_index=True)

#Keep only species
species=taxons.loc[taxons['Rank']=='SPECIES']

#Keep only those that were reported in publications
publ=species.dropna(subset='Authority')
publ=publ.dropna(subset='Genus name')
num_species=publ['Genus name'].value_counts()

#Add parasitic fungal species 
manual=pd.read_csv('/'.join([wdir,'manually_added_taxa.csv']),sep=';').dropna(subset='Species')
found=manual['Species'].tolist()
found=['"'+f+'"' for f in found]
found=','.join(found)

command = f"{condaenvs}/{conda_env}/bin/datasets download taxonomy taxon {found} --filename additional.zip"
sp.run(command, shell=True)
sp.run("unzip additional.zip -d additional", shell=True)
tax=pd.read_csv('/'.join(['additional','ncbi_dataset','data','taxonomy_summary.tsv']),sep='\t')

#Remove cordyceps since it is a mushroom
tax=tax.loc[~tax['Tax name'].str.contains('Cordyceps')]

tax=pd.concat([tax,publ]).drop_duplicates(subset='Tax name')
tax.to_csv('/'.join([wdir,'genomes_to_download.csv']),index=False)

#Download the genomes
tax=pd.read_csv('/'.join([wdir,'genomes_to_download.csv']))

def download_genomes(tax):
    for ix,f in tax.iterrows():
        command = f"{condaenvs}/{conda_env}/bin/datasets download genome taxon {int(f.Taxid)} --include genome --assembly-source 'RefSeq' --reference --filename genomes/{int(f.Taxid)}.zip"
        sp.run(command, shell=True)
        if os.path.isfile("/".join([wdir,'genomes',str(int(f.Taxid))+'.zip'])):
            sp.run(f"unzip genomes/{int(f.Taxid)}.zip -d genomes/{int(f.Taxid)}", shell=True)

download_genomes(tax)

##Format genome files to concatenate all contigs into one and add genome size summary to the taxonomy file

dg=os.listdir("/".join([wdir,'genomes']))
dg=[d for d in dg if '.zip' not in d]

#Make a summary of all assemblies:

final=pd.DataFrame()
final['Taxid']=dg

for ix,g in final.iterrows():
    summary=pd.read_json("/".join([wdir,'genomes',g.Taxid,'ncbi_dataset','data','assembly_data_report.jsonl']),lines=True)
    org=summary['organism'][0]
    stats=summary['assemblyStats'][0]
    final.at[ix,'Species']=org['organismName']
    final.at[ix,'Accession']=summary['accession'][0]
    final.at[ix,'Species']=org['organismName']
    final.at[ix,'GC']=stats['gcPercent']
    final.at[ix,'NumContigs']=stats['numberOfContigs']
    final.at[ix,'GenomeLength']=stats['totalSequenceLength']
    
    fasta=os.listdir("/".join([wdir,'genomes',g.Taxid,'ncbi_dataset','data',summary['accession'][0]]))[0]
    shutil.copy2("/".join([wdir,'genomes',g.Taxid, 'ncbi_dataset','data',summary['accession'][0],fasta]), "/".join([wdir,'fungal_genomes','unconcat',summary['accession'][0]+'.fna']))

final['Taxid']=final['Taxid'].astype(int)
final=final.merge(tax[['Taxid','Phylum name','Class name','Order name','Family name','Genus name','Species name']],on='Taxid',how='left')

final.to_csv('/'.join([wdir,'final_genomes_summary.csv']),index=False)

## Concatenate the files

files=os.listdir("/".join([wdir,'fungal_genomes','unconcat']))
outdir="/".join([wdir,'fungal_genomes','concatenated'])
os.mkdir(outdir)

cstr='N'*10
for fi in files:
        qfasta=SeqIO.parse('/'.join([wdir,'fungal_genomes',fi]),'fasta')
        concat=str('')
        for fasta in qfasta:
            concat=cstr.join([concat,str(fasta.seq)])
        record = SeqRecord(Seq(concat), name=fi.replace('.fna',''),
                           id=fi.replace('.fna',''), description='')
        SeqIO.write(record,'/'.join([outdir,fi]),'fasta')

#Remove cordyceps since it is a mushroom
final=final.loc[~final['Species'].str.contains('Cordyceps')]

##Create profiles - even number of reads

def generate_profiles(final,num_fungal, num_profiles, numreads,toincl,profdir,prefix):
    tochoose=[d for d in final['Taxid'].tolist() if d not in toincl]
    
    #If there is not enough of species for several different profiles
    if num_fungal*num_profiles>len(final):
        species=[]
        for a in range(num_profiles):
            species.extend(random.sample(tochoose,k=(num_fungal-len(toincl))))
    elif num_fungal==len(final):
        species=final['Taxid'].tolist()
    else: #otherwise
        species=random.sample(tochoose,k=(num_fungal-len(toincl))*num_profiles)
        
    for i in range(num_profiles):
        ind=i*(num_fungal-len(toincl))
        if i!=0:
            ind+=1
        splist=species[ind:ind+num_fungal-len(toincl)]
        profile=final.loc[final['Taxid'].isin(toincl)]
        profile=pd.concat([profile,final.loc[final['Taxid'].isin(splist)]])
        profile=profile[['Taxid','Species','Accession','GenomeLength']]
        profile['NumReads']=numreads
        profile['AvgCoverage']=profile['NumReads']*150/profile['GenomeLength']
        profile=profile.drop_duplicates(subset='Accession')
        profile.to_csv("/".join([profdir,prefix+'_'+str(i+1)+'.csv']),index=False)
   

profdir="/".join([wdir,'fungal_genomes','profiles'])
os.mkdir(profdir)

toincl=[5476,4959,28985,4932,76775]

num_fungal=len(final)
num_profiles=1
numreads=100000

generate_profiles(final,num_fungal,num_profiles,numreads,toincl,profdir,'complete')    

##Create profiles - even coverage of genomes

def create_coverage_profile(profdir,coverage):
    proflist=os.listdir(profdir)
    outdir=profdir+'_equal_cov'
    os.mkdir(outdir)
    for p in proflist:
        prof=pd.read_csv('/'.join([profdir,p]))
        prof['NumReads']=round(prof['GenomeLength']*coverage/150)
        prof['AvgCoverage']=coverage
        prof.to_csv('/'.join([outdir,p.replace('.csv','_EC.csv')]),index=False)

create_coverage_profile(profdir,2)
        
    


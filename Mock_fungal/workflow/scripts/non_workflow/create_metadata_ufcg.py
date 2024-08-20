#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 10:01:45 2024

#Make a metadata file for UFCG core gene tree

@author: ec-ekateria
"""
import pandas as pd
import os

profdir='FULL_PATH_TO/Mock_fungal/mock_profiles'

full=pd.read_csv('/'.join([profdir,'profiles/complete_1.csv']))
ncbi=pd.read_csv('/'.join([profdir.replace('/fungal_genomes',''),'final_genomes_summary.csv']))

meta=full[['Accession']]
meta['accession']=meta['Accession']
meta=meta.rename(columns={'Accession':'filename'})
meta['filename']=meta['filename'].apply(lambda row: row+'.fna')

ncbi=ncbi.fillna('')
ncbi['taxonomy']=ncbi.apply(lambda row: ';'.join(['Fungi',row['Phylum name'],row['Class name'], row['Order name'],row['Family name'],
                                                  row['Genus name'],row['Species name']]), axis=1)
ncbi=ncbi.rename(columns={'Accession':'accession'})

meta=meta.merge(ncbi[['accession','Species','taxonomy']], on='accession',how='left')
meta=meta.rename(columns={'Species':'label'})

meta.to_csv('/'.join([profdir,'metadata_ufcg.tsv']),index=False, sep='\t')

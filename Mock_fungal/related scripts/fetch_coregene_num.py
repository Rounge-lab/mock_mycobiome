#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 14:03:50 2024

Fetch number of core genes detected by the UFCG

@author: ec-ekateria
"""

import pandas as pd
import json
import os
import seaborn as sb
from scipy.stats import pearsonr

profdir='FULL_PATH_TO/Mock_fungal/mock_profiles'

genomes=os.listdir('/'.join([profdir,'ufcg_coregenes']))

taxonomy=pd.read_csv('/'.join([profdir,'metadata_ufcg_full.csv']),sep='\t')

genomes=pd.DataFrame(genomes)
genomes.columns=['UCG']
targetcg=61
for ix,g in genomes.iterrows():
    with open('/'.join([profdir,'ufcg_coregenes',g.UCG]),'r') as file:
        data=json.load(file)
    
    genomes.at[ix,'label']=data['genome_info']['label']
    genomes.at[ix,'DetectedCoreGenes']=len(data['data'])
    
genomes['PercCoreGenes']=genomes['DetectedCoreGenes']/targetcg*100

taxonomy=taxonomy.merge(genomes[['label','DetectedCoreGenes','PercCoreGenes']], on='label',how='left')

summary=pd.read_csv('/'.join([profdir,'final_genomes_summary.csv']))
summary=summary.rename(columns={'Accession':'accession'})

taxonomy=taxonomy.merge(summary[['accession','GC','NumContigs','GenomeLength']], on='accession', how='left')

taxonomy.to_csv('/'.join([profdir,'metadata_ufcg_full.csv']),sep='\t')

fig=sb.scatterplot(data=taxonomy,x='GenomeLength',y='PercCoreGenes')
fig.set(xlabel='Genome length, bp',ylabel='Core genes, %', xscale='log',xlim=[2*10**6, 10**8])

r,p=pearsonr(taxonomy['GenomeLength'],taxonomy['PercCoreGenes'])

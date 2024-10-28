#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 11:11:49 2024

Create dataset descriptions for iTOL tree visualization

@author: ec-ekateria
"""

import pandas as pd
import os

profdir='FULL_PATH_TO/Mock_fungal/mock_profiles' # NB! CHANGE TO LOCAL PATH
prof=['small_1','small_2','small_3','median_1','median_2','median_3','large_1','large_2','complete_1']

os.mkdir('/'.join([profdir,'datasets_itol']))
ds=['small','median','large','complete']
tools=['Kraken','Metaphlan','Euk','HMS','FunOMIC','MiCoP']

for d in ds:
    for p in prof:
        if d in p:
            file=pd.read_csv('/'.join([profdir,'profiles',p+'.csv']))
            file=file[['Species','KrakenDetected_Species','EukDetected_Species','MetaphlanDetected_Species', 'HMSDetected_Species','FunOMICDetected_Species','MiCoPDetected_Species']]
            file[p]=file.apply(lambda row: (row == 'Yes').sum(), axis=1)
            file[p]=file[p]*10+10
            file[p]=file[p].astype(str)
            file=file[['Species', p]]
            if 'dataset' not in globals():
                dataset=file.copy()
            else:
                dataset=dataset.merge(file,on='Species',how='outer')
    dataset=dataset.fillna(0)
    dataset['Species']=dataset['Species'].apply(lambda row: row.replace(' ','_'))
    dataset['Species']=dataset['Species'].apply(lambda row: '_subhashii' if '_subhashii' in row else row)
    dataset.to_csv('/'.join([profdir,'datasets_itol_new','dataset_'+d+'.txt']),sep=',',index=False)
    del dataset

def count_yes_no(row):
    counts = row.value_counts()
    return pd.Series({'Yes': counts.get('Yes', 0), 'No': counts.get('No', 0)})

for t in tools:
    for p in prof:
        file=pd.read_csv('/'.join([profdir,'profiles',p+'.csv']))
        file=file[['Species',t+'Detected_Species']]
        if 'dataset' not in globals():
            dataset=file.copy()
        else:
            dataset=dataset.merge(file,on='Species',how='outer')
        dataset=dataset.rename(columns={t+'Detected_Species':p})
        
    #Count number of Yes and No
    yes_no = dataset.apply(count_yes_no, axis=1)
    dataset = dataset[['Species']].join(yes_no)
    dataset[t]=dataset.apply(lambda row: 100 if row.No==0 and row.Yes>0 else (50 if row.No>0 and row.Yes>0 else 25), axis=1)
        
    if 'final' not in globals():
        final=dataset[['Species',t]].copy()
    else:
        final=final.merge(dataset[['Species',t]],on='Species',how='outer')

final['Species']=final['Species'].apply(lambda row: row.replace(' ','_'))
final['Species']=final['Species'].apply(lambda row: '_subhashii' if '_subhashii' in row else row)
final.to_csv('/'.join([profdir,'datasets_itol','Tools_dataset2.txt']),sep=',',index=False)


#create taxonomy color dataset
                
taxonomy=pd.read_csv('/'.join([profdir,'metadata_ufcg.csv']),sep='\t')
taxonomy[['domain','phylum','class','order','family','genus','species']] = taxonomy['taxonomy'].str.split(';', expand=True)
taxonomy=taxonomy[['label','class']]
taxonomy.loc[taxonomy['label'].str.contains('Encephalitozoon'),'class']='Microsporea'
classes=taxonomy['class'].unique().tolist()
colors=['#c7c6c8']*len(classes)
colcodes=dict(zip(classes,colors))

colcodes['Saccharomycetes']='#f4c9ae'
colcodes['Sordariomycetes']='#aeccf4'
colcodes['Wallemiomycetes']='#e7d4b5'
colcodes['Eurotiomycetes']='#cb6de6'
colcodes['Tremellomycetes']='#ffd866'
colcodes['Dothideomycetes']='#0097b2'
colcodes['Malasseziomycetes']='#ff5859'
colcodes['Pneumocystidiomycetes']='#ff66c4'
colcodes['Mucoromycetes']='#03bf62'


taxonomy['range']='range'
taxonomy['color']=taxonomy['class'].apply(lambda row: colcodes[row])
taxonomy['label']=taxonomy['label'].apply(lambda row: row.replace(' ','_'))
taxonomy['label']=taxonomy['label'].apply(lambda row: '_subhashii' if '_subhashii' in row else row)
taxonomy=taxonomy.drop(columns='class')
taxonomy.to_csv('/'.join([profdir,'datasets_itol','class_colors.txt']),sep=' ',index=False)


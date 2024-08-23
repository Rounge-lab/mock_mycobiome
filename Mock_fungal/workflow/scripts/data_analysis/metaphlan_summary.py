#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 13:10:32 2024

Analyze metaphlan reports

@author: ec-ekateria
"""

import pandas as pd
import os
import seaborn as sb
import matplotlib.pyplot as plt

global suffix 

profdir='/FULL_PATH_TO/Mock_fungal/mock_profiles/profiles_equal_cov' # if equal coverage; else profiles AND remember to change the last line of the script; NB! CHANGE TO LOCAL PATH
mdir='/FULL_PATH_TO/Mock_fungal/data/metaphlan/equal_coverage_metagenomes' # NB! CHANGE TO LOCAL PATH
taxonomy=pd.read_csv('/FULL_PATH_TO/Mock_fungal/mock_profiles/final_genomes_summary.csv') # NB! CHANGE TO LOCAL PATH

suffix='_metarep.txt'
replist=os.listdir(mdir)
replist=[l for l in replist if suffix in l]


def summary_metaphlan(replist,taxonomy):
    summary=pd.DataFrame()
    allprofiles=pd.DataFrame()

    for r in replist:
        rep=pd.read_csv('/'.join([mdir,r]),sep='\t',skiprows=4)
        rep=rep.rename(columns={'#clade_name':'Taxonomy','relative_abundance': 'RelAb'})
        rep=rep.drop(columns={'NCBI_tax_id','additional_species'})

        prof=pd.read_csv('/'.join([profdir,r.replace(suffix,'.csv')]))
        
        if 'Species name' not in prof.columns.tolist():
            prof=prof.merge(taxonomy[['Taxid','GC','NumContigs','Family name','Genus name','Species name']], on='Taxid',how='left')
        
        rep['Species name']=rep['Taxonomy'].apply(lambda row: row.split('s__')[-1].split('|t__')[0] if 's__' in row else None)
        rep['Species name']=rep['Species name'].apply(lambda row: ' '.join(row.split('_')[0:2]) if row is not None else row)

        rep['Genus name']=rep['Taxonomy'].apply(lambda row: row.split('g__')[-1].split('|s__')[0] if 'g__' in row else None)
        rep['Family name']=rep['Taxonomy'].apply(lambda row: row.split('f__')[-1].split('|g__')[0] if 'f__' in row else None)
        
        colname=r.replace(suffix,'')
        summary.at['GeneratedReads',colname]=prof['NumReads'].sum()
        summary.at['NumSpecies',colname]=len(prof['Species name'].unique().tolist())
        summary.at['NumGenera',colname]=len(prof['Genus name'].unique().tolist())
        summary.at['NumFamily',colname]=len(prof['Family name'].unique().tolist())

        if 'Homo sapiens' in rep['Species name'].tolist():
            summary.at['ClassifiedHuman',colname]=rep.loc[rep['Species name']=='Homo sapiens','NumReads'].tolist()[0]
        
        rep['TDSpecies']=rep['Species name'].apply(lambda row: 'Yes' if row in prof['Species name'].tolist() else 'No')
        rep['TDGenus']=rep['Genus name'].apply(lambda row: 'Yes' if row in prof['Genus name'].tolist() else 'No')
        rep['TDFamily']=rep['Family name'].apply(lambda row: 'Yes' if row in prof['Family name'].tolist() else 'No')
        
        prof['MetaphlanDetected_Species']=prof['Species name'].apply(lambda row: 'Yes' if row in rep['Species name'].tolist() else 'No')
        prof['MetaphlanDetected_Genus']=prof['Genus name'].apply(lambda row: 'Yes' if row in rep['Genus name'].tolist() else 'No')
        prof['MetaphlanDetected_Family']=prof['Family name'].apply(lambda row: 'Yes' if row in rep['Family name'].tolist() else 'No')
        
        summary.at['TDSpecies',colname]=len(prof.loc[prof['MetaphlanDetected_Species']=='Yes','Species name'].unique().tolist())
        summary.at['TDGenera',colname]=len(prof.loc[prof['MetaphlanDetected_Genus']=='Yes','Genus name'].unique().tolist())
        summary.at['TDFamily',colname]=len(prof.loc[prof['MetaphlanDetected_Family']=='Yes','Family name'].unique().tolist())
        
        summary.at['FDSpecies',colname]=len(rep.dropna(subset='Species name').loc[rep['TDSpecies']=='No','Species name'].unique().tolist())
        summary.at['FDGenera',colname]=len(rep.dropna(subset='Genus name').loc[rep['TDGenus']=='No','Genus name'].unique().tolist())
        summary.at['FDFamily',colname]=len(rep.dropna(subset='Family name').loc[rep['TDFamily']=='No','Family name'].unique().tolist())
    
        prof.to_csv('/'.join([profdir,r.replace(suffix,'.csv')]),index=False)
        prof['MockCommunity']=colname
        allprofiles=pd.concat([allprofiles,prof],ignore_index=True)
    
    return summary,allprofiles

summary,allprofiles=summary_metaphlan(replist,taxonomy)

summary=summary.T.reset_index().rename(columns={'index':'MockCommunity'})

level=['Species','Genera', 'Family']
for l in level:
    summary['Recall'+l]=summary['TD'+l]/summary['Num'+l]*100 #how many of those that were there, were detected
    summary['Precision'+l]=summary['TD'+l]/(summary['TD'+l]+summary['FD'+l])*100 #how many of those that were detected, truly were added there
    
summary.to_csv('/'.join([mdir.replace('equal_coverage_metagenomes','summaries'),'EqualCoverage_Summary.csv'])) #Change for Equal reads




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 11:09:33 2024

Analyze EukDetect reports

@author: ec-ekateria
"""

import pandas as pd
import os
import seaborn as sb
import matplotlib.pyplot as plt

profdir='/FULL_PATH_TO/Mock_fungal/mock_profiles/profiles_equal_cov' # if equal coverage; else profiles AND remember to change the last line of the script; NB! CHANGE TO LOCAL PATH
eukdir='/FULL_PATH_TO/Mock_fungal/data/eukdetect/equal_coverage_metagenomes' # NB! CHANGE TO LOCAL PATH
taxonomy=pd.read_csv('FULL_PATH_TO/Mock_fungal/mock_profiles/final_genomes_summary.csv') # NB! CHANGE TO LOCAL PATH

replist=os.listdir(eukdir)
replist=[l for l in replist if '_filtered_hits_table' in l]

def summary_eukdetect(replist,taxonomy):
    summary=pd.DataFrame()
    allprofiles=pd.DataFrame()

    for r in replist:
        rep=pd.read_csv('/'.join([eukdir,r]),sep='\t')

        prof=pd.read_csv('/'.join([profdir,r.replace('_filtered_hits_table.txt','.csv')]))
        
        if 'Species name' not in prof.columns.tolist():
            prof=prof.merge(taxonomy[['Taxid','GC','NumContigs','Family name','Genus name','Species name']], on='Taxid',how='left')
        
        rep['Species name']=rep['Lineage'].apply(lambda row: row.split('|species-')[-1] if '|species-' in row else None)
        rep['Species name']=rep['Species name'].str.replace('_',' ')
        rep['Genus name']=rep['Lineage'].apply(lambda row: row.split('|genus-')[-1].split('|species-')[0] if '|genus-' in row else None)
        rep['Family name']=rep['Lineage'].apply(lambda row: row.split('|family-')[-1].split('|genus-')[0] if '|family-' in row else None)
        
        colname=r.replace('_filtered_hits_table.txt','')
        summary.at['GeneratedReads',colname]=prof['NumReads'].sum()
        summary.at['NumSpecies',colname]=len(prof['Species name'].unique().tolist())
        summary.at['NumGenera',colname]=len(prof['Genus name'].unique().tolist())
        summary.at['NumFamily',colname]=len(prof['Family name'].unique().tolist())

        summary.at['ClassifiedReads',colname]=rep['Read_counts'].sum()
        summary.at['UnclassifiedReads',colname]=summary.loc['GeneratedReads',colname]-summary.loc['ClassifiedReads',colname]

        rep['TDSpecies']=rep['Species name'].apply(lambda row: 'Yes' if row in prof['Species name'].tolist() else 'No')
        rep['TDGenus']=rep['Genus name'].apply(lambda row: 'Yes' if row in prof['Genus name'].tolist() else 'No')
        rep['TDFamily']=rep['Family name'].apply(lambda row: 'Yes' if row in prof['Family name'].tolist() else 'No')
        
        #some of lineages do not have information about the family, keep 'yes' if species/genus was detected
        rep['TDFamily']=rep.apply(lambda row: 'Yes' if row.TDSpecies=='Yes' or row.TDGenus=='Yes' else row.TDFamily, axis=1)
        
        prof['EukDetected_Species']=prof['Species name'].apply(lambda row: 'Yes' if row in rep['Species name'].tolist() else 'No')
        prof['EukDetected_Genus']=prof['Genus name'].apply(lambda row: 'Yes' if row in rep['Genus name'].tolist() else 'No')
        prof['EukDetected_Family']=prof['Family name'].apply(lambda row: 'Yes' if row in rep['Family name'].tolist() else 'No')
        
        summary.at['TDSpecies',colname]=len(prof.loc[prof['EukDetected_Species']=='Yes','Species name'].unique().tolist())
        summary.at['TDGenera',colname]=len(prof.loc[prof['EukDetected_Genus']=='Yes','Genus name'].unique().tolist())
        summary.at['TDFamily',colname]=len(prof.loc[prof['EukDetected_Family']=='Yes','Family name'].unique().tolist())
        
        summary.at['FDSpecies',colname]=len(rep.dropna(subset='Species name').loc[rep['TDSpecies']=='No','Species name'].unique().tolist())
        summary.at['FDGenera',colname]=len(rep.dropna(subset='Genus name').loc[rep['TDGenus']=='No','Genus name'].unique().tolist())
        summary.at['FDFamily',colname]=len(rep.dropna(subset='Family name').loc[rep['TDFamily']=='No','Family name'].unique().tolist())
    
        prof.to_csv('/'.join([profdir,r.replace('_filtered_hits_table.txt','.csv')]),index=False)
        prof['MockCommunity']=colname
        allprofiles=pd.concat([allprofiles,prof],ignore_index=True)
    
    return summary,allprofiles

summary,allprofiles=summary_eukdetect(replist,taxonomy)

summary=summary.T.reset_index().rename(columns={'index':'MockCommunity'})

summary['PropEukaryota']=summary['ClassifiedReads']/summary['GeneratedReads']*100
summary['PropUnclassified']=summary['UnclassifiedReads']/summary['GeneratedReads']*100


level=['Species','Genera', 'Family']
for l in level:
    summary['Recall'+l]=summary['TD'+l]/summary['Num'+l]*100 #how many of those that were there, were detected
    summary['Precision'+l]=summary['TD'+l]/(summary['TD'+l]+summary['FD'+l])*100 #how many of those that were detected, truly were added there
    
summary.to_csv('/'.join([eukdir.replace('equal_coverage_metagenomes','summaries'),'EqualCoverage_Summary.csv'])) #Change for Equal reads


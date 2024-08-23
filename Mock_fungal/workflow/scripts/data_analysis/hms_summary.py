#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 13:20:55 2024

Analyze HMS reports

@author: ec-ekateria
"""

import pandas as pd
import os
import seaborn as sb
import matplotlib.pyplot as plt

profdir='/FULL_PATH_TO/Mock_fungal/mock_profiles/profiles' #if equal reads; else profiles_equal_cov AND remember to change the last line of the script; NB! CHANGE TO LOCAL PATH
hmsdir='/FULL_PATH_TO/Mock_fungal/data/HMS/equal_reads_metagenomes' # NB! CHANGE TO LOCAL PATH
taxonomy=pd.read_csv('/FULL_PATH_TO/Mock_fungal/mock_profiles/final_genomes_summary.csv') # NB! CHANGE TO LOCAL PATH

replist=os.listdir(hmsdir)

def summary_hms(replist,taxonomy):
    summary=pd.DataFrame()
    allprofiles=pd.DataFrame()

    for r in replist:
        
        rep_sp=pd.read_csv('/'.join([hmsdir,r,r+'_Species_level_results-%.txt']),sep='\t').reset_index().rename(columns={'index':'Species','Unnamed: 0':'count'})
        rep_g=pd.read_csv('/'.join([hmsdir,r,r+'_Genera_level_results-%.txt']),sep='\t').reset_index().rename(columns={'index':'Genus','V1':'count'})
        rep_f=pd.read_csv('/'.join([hmsdir,r,r+'_Family_level_results-Counts.txt']),sep='\t',header=None).rename(columns={0:'Family',1:'count'})

        #Filter family detection, remove family with reads <'0.1% of sum of classified'
        
        if 'Unclassified' in rep_f['Family'].tolist():
            unclas=rep_f.loc[rep_f['Family']=='Unclassified','count'].values[0]
        else:
            unclas=0
        clas=rep_f.loc[rep_f['Family']!='Unclassified','count'].sum()
        filt=clas*0.001
        rep_f=rep_f.loc[rep_f['count']>filt]
        prof=pd.read_csv('/'.join([profdir,r+'.csv']))
        
        if 'Species name' not in prof.columns.tolist():
            prof=prof.merge(taxonomy[['Taxid','GC','NumContigs','Family name','Genus name','Species name']], on='Taxid',how='left')
        
        rep_sp['Species']=rep_sp['Species'].apply(lambda row: row.replace('_',' '))
        
        colname=r
        summary.at['GeneratedReads',colname]=prof['NumReads'].sum()
        summary.at['NumSpecies',colname]=len(prof['Species name'].unique().tolist())
        summary.at['NumGenera',colname]=len(prof['Genus name'].unique().tolist())
        summary.at['NumFamily',colname]=len(prof['Family name'].unique().tolist())

        summary.at['ClassifiedReads',colname]=clas
        summary.at['UnclassifiedReads',colname]=unclas

        rep_sp['TDSpecies']=rep_sp['Species'].apply(lambda row: 'Yes' if row in prof['Species name'].tolist() else 'No')
        rep_g['TDGenus']=rep_g['Genus'].apply(lambda row: 'Yes' if row in prof['Genus name'].tolist() else 'No')
        rep_f['TDFamily']=rep_f['Family'].apply(lambda row: 'Yes' if row in prof['Family name'].tolist() else 'No')
        
        prof['HMSDetected_Species']=prof['Species name'].apply(lambda row: 'Yes' if row in rep_sp['Species'].tolist() else 'No')
        prof['HMSDetected_Genus']=prof['Genus name'].apply(lambda row: 'Yes' if row in rep_g['Genus'].tolist() else 'No')
        prof['HMSDetected_Family']=prof['Family name'].apply(lambda row: 'Yes' if row in rep_f['Family'].tolist() else 'No')
        
        #some of lineages do not have information about the family, keep 'yes' if species/genus was detected
        prof['HMSDetectedFamily']=prof.apply(lambda row: 'Yes' if row.HMSDetected_Species=='Yes' or row.HMSDetected_Genus=='Yes' else row.HMSDetected_Family, axis=1)
        
        summary.at['TDSpecies',colname]=len(prof.loc[prof['HMSDetected_Species']=='Yes','Species name'].unique().tolist())
        summary.at['TDGenera',colname]=len(prof.loc[prof['HMSDetected_Genus']=='Yes','Genus name'].unique().tolist())
        summary.at['TDFamily',colname]=len(prof.loc[prof['HMSDetected_Family']=='Yes','Family name'].unique().tolist())
        
        summary.at['FDSpecies',colname]=len(rep_sp.loc[rep_sp['TDSpecies']=='No'])
        summary.at['FDGenera',colname]=len(rep_g.loc[rep_g['TDGenus']=='No'])
        summary.at['FDFamily',colname]=len(rep_f.loc[rep_f['TDFamily']=='No'])
    
        prof.to_csv('/'.join([profdir,r+'.csv']),index=False)
        prof['MockCommunity']=colname
        allprofiles=pd.concat([allprofiles,prof],ignore_index=True)
    
    return summary,allprofiles

summary,allprofiles=summary_hms(replist,taxonomy)

summary=summary.T.reset_index().rename(columns={'index':'MockCommunity'})

level=['Species','Genera', 'Family']
for l in level:
    summary['Recall'+l]=summary['TD'+l]/summary['Num'+l]*100 #how many of those that were there, were detected
    summary['Precision'+l]=summary['TD'+l]/(summary['TD'+l]+summary['FD'+l])*100 #how many of those that were detected, truly were added there
    
summary.to_csv('/'.join([hmsdir.replace('equal_reads_metagenomes','summaries'),'EqualReads_Summary.csv'])) ## Remember to switch for equal coverage


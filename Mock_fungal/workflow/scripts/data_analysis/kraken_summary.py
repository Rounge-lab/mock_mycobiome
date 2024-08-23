#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 10:05:56 2024

Analyze kraken reports

@author: ec-ekateria
"""

import pandas as pd
import os
import seaborn as sb
import matplotlib.pyplot as plt

profdir='/FULL_PATH_TO/Mock_fungal/mock_profiles/profiles' #if equal reads; else profiles_equal_cov AND remember to change the last line of the script; NB! CHANGE TO LOCAL PATH
krdir='/FULL_PATH_TO/Mock_fungal/data/kraken/equal_reads_metagenomes' # NB! CHANGE TO LOCAL PATH
taxonomy=pd.read_csv('/FULL_PATH_TO/Mock_fungal/mock_profiles/final_genomes_summary.csv')# NB! CHANGE TO LOCAL PATH

replist=os.listdir(krdir)
replist=[l for l in replist if '_report' in l]


def summary_kraken(replist,taxonomy):
    summary=pd.DataFrame()
    allprofiles=pd.DataFrame()

    domains=['d__Eukaryota','d__Bacteria','d__Viruses']
    for r in replist:
        rep=pd.read_csv('/'.join([krdir,r]),sep='\t',header=None)
        rep=rep.rename(columns={0:'Taxonomy',1: 'NumReads'})

        prof=pd.read_csv('/'.join([profdir,r.replace('_report.txt','.csv')]))
        
        if 'Species name' not in prof.columns.tolist():
            prof=prof.merge(taxonomy[['Taxid','GC','NumContigs','Family name','Genus name','Species name']], on='Taxid',how='left')
        
        rep['Species name']=rep['Taxonomy'].apply(lambda row: row.split('s__')[-1] if 's__' in row else None)
        rep['Genus name']=rep['Taxonomy'].apply(lambda row: row.split('g__')[-1].split('|s__')[0] if 'g__' in row else None)
        rep['Family name']=rep['Taxonomy'].apply(lambda row: row.split('f__')[-1].split('|g__')[0] if 'f__' in row else None)
        
        summary.at['GeneratedReads',r.replace('_report.txt','')]=prof['NumReads'].sum()
        summary.at['NumSpecies',r.replace('_report.txt','')]=len(prof['Species name'].unique().tolist())
        summary.at['NumGenera',r.replace('_report.txt','')]=len(prof['Genus name'].unique().tolist())
        summary.at['NumFamily',r.replace('_report.txt','')]=len(prof['Family name'].unique().tolist())


        if 'd__' in rep.loc[0,'Taxonomy']:
            totreads=rep.loc[rep['Taxonomy'].isin(domains)]
            d=[dm for dm in domains if dm in totreads['Taxonomy'].tolist()]
            
            for dm in d:
                summary.at['Classified'+dm.replace('d__',''),r.replace('_report.txt','')]=totreads.loc[totreads['Taxonomy']==dm,'NumReads'].tolist()[0]
        
        if 'Homo sapiens' in rep['Species name'].tolist():
            summary.at['ClassifiedHuman',r.replace('_report.txt','')]=rep.loc[rep['Species name']=='Homo sapiens','NumReads'].tolist()[0]
    
        #Keep only fungi
        fungi=rep.loc[rep['Taxonomy'].str.contains('k__Fungi')]
        
        fungi['TDSpecies']=fungi['Species name'].apply(lambda row: 'Yes' if row in prof['Species name'].tolist() else 'No')
        fungi['TDGenus']=fungi['Genus name'].apply(lambda row: 'Yes' if row in prof['Genus name'].tolist() else 'No')
        fungi['TDFamily']=fungi['Family name'].apply(lambda row: 'Yes' if row in prof['Family name'].tolist() else 'No')
        
        prof['KrakenDetected_Species']=prof['Species name'].apply(lambda row: 'Yes' if row in fungi['Species name'].tolist() else 'No')
        prof['KrakenDetected_Genus']=prof['Genus name'].apply(lambda row: 'Yes' if row in fungi['Genus name'].tolist() else 'No')
        prof['KrakenDetected_Family']=prof['Family name'].apply(lambda row: 'Yes' if row in fungi['Family name'].tolist() else 'No')
        
        summary.at['TDSpecies',r.replace('_report.txt','')]=len(prof.loc[prof['KrakenDetected_Species']=='Yes','Species name'].unique().tolist())
        summary.at['TDGenera',r.replace('_report.txt','')]=len(prof.loc[prof['KrakenDetected_Genus']=='Yes','Genus name'].unique().tolist())
        summary.at['TDFamily',r.replace('_report.txt','')]=len(prof.loc[prof['KrakenDetected_Family']=='Yes','Family name'].unique().tolist())
        
        summary.at['FDSpecies',r.replace('_report.txt','')]=len(fungi.dropna(subset='Species name').loc[fungi['TDSpecies']=='No','Species name'].unique().tolist())
        summary.at['FDGenera',r.replace('_report.txt','')]=len(fungi.dropna(subset='Genus name').loc[fungi['TDGenus']=='No','Genus name'].unique().tolist())
        summary.at['FDFamily',r.replace('_report.txt','')]=len(fungi.dropna(subset='Family name').loc[fungi['TDFamily']=='No','Family name'].unique().tolist())
    
        #prof.to_csv('/'.join([profdir,r.replace('_report.txt','.csv')]),index=False)
        prof['MockCommunity']=r.replace('_report.txt','')
        allprofiles=pd.concat([allprofiles,prof],ignore_index=True)
    
    return summary,allprofiles

summary,allprofiles=summary_kraken(replist,taxonomy)

summary=summary.T.reset_index().rename(columns={'index':'MockCommunity'})

summary['PropEukaryota']=summary['ClassifiedEukaryota']/summary['GeneratedReads']*100
summary['PropHuman']=summary['ClassifiedHuman']/summary['GeneratedReads']*100
summary['PropBacteria']=summary['ClassifiedBacteria']/summary['GeneratedReads']*100
summary['PropViruses']=summary['ClassifiedViruses']/summary['GeneratedReads']*100
summary['PropUnclassified']=100-summary['PropEukaryota']-summary['PropBacteria']-summary['PropViruses']


level=['Species','Genera', 'Family']
for l in level:
    summary['Recall'+l]=summary['TD'+l]/summary['Num'+l]*100 #how many of those that were there, were detected
    summary['Precision'+l]=summary['TD'+l]/(summary['TD'+l]+summary['FD'+l])*100 #how many of those that were detected, truly were added there
    
summary.to_csv('/'.join([krdir.replace('equal_reads_metagenomes','summaries'),'EqualReads_Summary.csv'])) ## Remember to switch for equal coverage


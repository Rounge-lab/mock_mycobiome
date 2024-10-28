#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 10:55:36 2024

Analyze MiCoP reports 

@author: ec-ekateria
"""

import pandas as pd
import os
import seaborn as sb
import matplotlib.pyplot as plt

bac='99'
mode=f'plusbacteria_metagenomes{bac}'
sumname=f'PlusBacteria{bac}_Summary.csv'

p='profiles' if '_reads' or 'plusbacteria' in mode else 'profiles_equal_cov' 
profdir=f'FULL_PATH_TO/Mock_fungal/mock_profiles/{p}'
micdir=f'FULL_PATH_TO/Mock_fungal/data/micop/{mode}'
taxonomy=pd.read_csv('FULL_PATH_TO/Mock_fungal/mock_profiles/final_genomes_summary.csv')

replist=os.listdir(micdir)
replist=[l for l in replist if '_abundance.txt' in l]


def summary_micop(replist,taxonomy):
    summary=pd.DataFrame()
    allprofiles=pd.DataFrame()

    for r in replist:
        r=r.replace('_abundance.txt','')
        rep=pd.read_csv('/'.join([micdir,r+'_abundance.txt']),sep='\t')
        
        rep=rep[['RANK','TAXPATHSN','PERCENTAGE']]
        rep['TAXPATHSN']=rep['TAXPATHSN'].apply(lambda row: row.split('|')[-1])
        
        ranks=['family','genus','species']
        rep=rep.loc[rep['RANK'].isin(ranks)]

        prof=pd.read_csv('/'.join([profdir,r.replace('b','')+'.csv']))
        
        if 'Species name' not in prof.columns.tolist():
            prof=prof.merge(taxonomy[['Taxid','GC','NumContigs','Family name','Genus name','Species name']], on='Taxid',how='left')
        
        summary.at['GeneratedReads',r]=prof['NumReads'].sum()
        summary.at['NumSpecies',r]=len(prof['Species name'].unique().tolist())
        summary.at['NumGenera',r]=len(prof['Genus name'].unique().tolist())
        summary.at['NumFamily',r]=len(prof['Family name'].unique().tolist())

        rep['TDSpecies']=rep['TAXPATHSN'].apply(lambda row: 'Yes' if row in prof['Species name'].tolist() else 'No')
        rep['TDGenus']=rep['TAXPATHSN'].apply(lambda row: 'Yes' if row in prof['Genus name'].tolist() else 'No')
        rep['TDFamily']=rep['TAXPATHSN'].apply(lambda row: 'Yes' if row in prof['Family name'].tolist() else 'No')
        
        prof['MiCoPDetected_Species']=prof['Species name'].apply(lambda row: 'Yes' if row in rep['TAXPATHSN'].tolist() else 'No')
        prof['MiCoPDetected_Genus']=prof['Genus name'].apply(lambda row: 'Yes' if row in rep['TAXPATHSN'].tolist() else 'No')
        prof['MiCoPDetected_Family']=prof['Family name'].apply(lambda row: 'Yes' if row in rep['TAXPATHSN'].tolist() else 'No')
        
        summary.at['TDSpecies',r]=len(prof.loc[prof['MiCoPDetected_Species']=='Yes','Species name'].unique().tolist())
        summary.at['TDGenera',r]=len(prof.loc[prof['MiCoPDetected_Genus']=='Yes','Genus name'].unique().tolist())
        summary.at['TDFamily',r]=len(prof.loc[prof['MiCoPDetected_Family']=='Yes','Family name'].unique().tolist())
        
        summary.at['FDSpecies',r]=len(rep.loc[rep['RANK']=='species'].loc[rep['TDSpecies']=='No','TAXPATHSN'].unique().tolist())
        summary.at['FDGenera',r]=len(rep.loc[rep['RANK']=='genus'].loc[rep['TDGenus']=='No','TAXPATHSN'].unique().tolist())
        summary.at['FDFamily',r]=len(rep.loc[rep['RANK']=='family'].loc[rep['TDFamily']=='No','TAXPATHSN'].unique().tolist())
    
        prof.to_csv('/'.join([profdir,r+f'_{bac}.csv']),index=False)
        prof['MockCommunity']=r
        allprofiles=pd.concat([allprofiles,prof],ignore_index=True)
    
    return summary,allprofiles

summary,allprofiles=summary_micop(replist,taxonomy)

summary=summary.T.reset_index().rename(columns={'index':'MockCommunity'})

level=['Species','Genera', 'Family']
for l in level:
    summary['Recall'+l]=summary['TD'+l]/summary['Num'+l]*100 #how many of those that were there, were detected
    summary['Precision'+l]=summary['TD'+l]/(summary['TD'+l]+summary['FD'+l])*100 #how many of those that were detected, truly were added there
    
summary.to_csv('/'.join([micdir.replace(mode,'summaries'),sumname]))

# props=pd.melt(summary[['NumSpecies','PropEukaryota','PropBacteria','PropViruses','PropUnclassified']],id_vars='NumSpecies')
# props['variable']=props['variable'].str.replace('Prop','')
# fig=sb.lineplot(data=props, y='value',x='NumSpecies',hue='variable',palette=['#1c6462','#6b519d','#ff66c4','#7ed957'],legend=False)
# sb.scatterplot(data=props, y='value',x='NumSpecies',hue='variable',palette=['#1c6462','#6b519d','#ff66c4','#7ed957'])
# fig.set(yscale='log',xlabel='Number of species',ylabel='Reads, %')
# legend = fig.legend_
# legend.set_title(None)
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

prec=pd.melt(summary[['NumSpecies','PrecisionSpecies','PrecisionGenera','PrecisionFamily']],id_vars='NumSpecies')
prec['variable']=prec['variable'].str.replace('Precision','')
prec['variable']=prec['variable'].apply(lambda row: 'Genus' if row=='Genera' else row)

recall=pd.melt(summary[['NumSpecies','RecallSpecies','RecallGenera','RecallFamily']],id_vars='NumSpecies')
recall['variable']=recall['variable'].str.replace('Recall','')
recall['variable']=recall['variable'].apply(lambda row: 'Genus' if row=='Genera' else row)

fig,ax=plt.subplots(2,1)
sb.lineplot(data=prec, y='value',x='NumSpecies',hue='variable',palette=['#6b519d','#ff66c4','#7ed957'],legend=False,ax=ax[0])
sb.scatterplot(data=prec, y='value',x='NumSpecies',hue='variable',palette=['#6b519d','#ff66c4','#7ed957'],ax=ax[0])
ax[0].set(xlabel='',ylabel='Precision, %')
legend = ax[0].legend_
legend.set_title(None)

sb.lineplot(data=recall, y='value',x='NumSpecies',hue='variable',palette=['#6b519d','#ff66c4','#7ed957'],legend=False,ax=ax[1])
sb.scatterplot(data=recall, y='value',x='NumSpecies',hue='variable',palette=['#6b519d','#ff66c4','#7ed957'],ax=ax[1])
ax[1].set(xlabel='Number of species',ylabel='Recall, %')
legend = ax[1].legend_
legend.set_title(None)




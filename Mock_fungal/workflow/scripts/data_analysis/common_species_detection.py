#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:04:10 2024

Check detection of five always included species
#Executed only after detection of species in mock communities by all tools has been added to the profile files

@author: ec-ekateria
"""

import pandas as pd


profdir='FULL_PATH_TO/Mock_fungal/mock_profiles' # NB! CHANGE TO LOCAL PATH
resfol='FULL_PATH_TO/Mock_fungal/data/all_methods_summary' # NB! CHANGE TO LOCAL PATH

modes={'profiles':'EqualReads', 'profiles_equal_cov':'EqualCoverage'}
prof=['small_1','small_2','small_3','median_1','median_2','median_3','large_1','large_2','complete_1']

bacteria=['Saccharomyces cerevisiae', 'Debaryomyces hansenii', 'Malassezia restricta', 'Candida albicans', 'Kluyveromyces lactis']

def combine_profiles(prof,modes):
    allbact=pd.DataFrame()
    global profdir, bacteria
    for m in modes:
        for p in prof:
            if m=='profiles':
                query=pd.read_csv('/'.join([profdir,m,p+'.csv']))
            else:
                query=pd.read_csv('/'.join([profdir,m,p+'_EC.csv']))
            genera=query['Genus name']
            query['NumSp_from_genus']=query['Genus name'].apply(lambda row: len([g for g in genera if g==row]))
            query=query.loc[query['Species name'].isin(bacteria)]
            query=query[['Species name','KrakenDetected_Species','KrakenDetected_Genus','KrakenDetected_Family',
                         'EukDetected_Species','EukDetected_Genus','EukDetected_Family',
                         'MetaphlanDetected_Species','MetaphlanDetected_Genus','MetaphlanDetected_Family','NumSp_from_genus']]
            
            query = pd.melt(query, id_vars=['Species name', 'NumSp_from_genus'], value_vars=['KrakenDetected_Species','KrakenDetected_Genus',
                                                                                             'KrakenDetected_Family','EukDetected_Species','EukDetected_Genus','EukDetected_Family',
                                                                                             'MetaphlanDetected_Species','MetaphlanDetected_Genus','MetaphlanDetected_Family'],
                            var_name='Val', value_name='Detection')
            query[['Tool','Level']]=query['Val'].str.split('_',expand=True)
            query['Tool']=query['Tool'].str.replace('Detected','')
            query['Tool']=query['Tool'].str.replace('Euk','EukDetect')
            query['ComType']=modes[m]
            query['Community']=p
            query=query.drop(columns='Val')
            allbact=pd.concat([allbact,query],ignore_index=True)
            
    return allbact

allbact=combine_profiles(prof,modes)

def find_detection(allbact,level):
    
    query=allbact.loc[allbact['Level']==level]
    query=query.loc[query['Detection']=='No']
    return query

#Find all cases where genus was not detected 
genus=find_detection(allbact,'Genus') #-> the only genus that was not detected is Kluyveromyces (all communities, both modes, by Metaphlan only)

#Find all cases where family was not detected 
family=find_detection(allbact,'Family') #-> empty dataframe; family was detected

#Find all cases where species was not detected
species=find_detection(allbact,'Species')
species=species.loc[species['Species name']!='Kluyveromyces lactis'] #remove the already known case
species.to_csv('/'.join([resfol,'NotDetected_CommonSpecies.csv']),index=False)

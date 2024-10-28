#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:38:57 2024

Describe detection of species in ER/EC context by different tools; compare detection of five commonly included species

@author: ec-ekateria
"""

import pandas as pd

profdir='FULL_PATH_TO/Mock_fungal/mock_profiles' # NB! CHANGE TO LOCAL PATH
metadb=pd.read_csv('/PATH_TO_DATABASES/db/biobakery/jun24/metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202307_species.txt',sep='\t',header=None) # NB! CHANGE TO LOCAL PATH
krdb=pd.read_csv('/PATH_TO_DATABASES/db/kraken-PlusPF/library_report.tsv',sep='\t') # NB! CHANGE TO LOCAL PATH
eukdb=pd.read_csv('/PATH_TO_TOOLS/EukDetect/eukdb/busco_taxid_link.txt',sep='\t',header=None) # NB! CHANGE TO LOCAL PATH
hmsdb=pd.read_csv('/PATH_TO_TOOLS/MScan2.0/HMS/var/HMS_taxonomy.txt',sep='\t') # NB! CHANGE TO LOCAL PATH
fundb=pd.read_csv('/PATH_TO_DATABASES/db/FunOMIC/FunOMIC-T/FunOMIC-Tv1/taxonomy.txt',sep='\t',header=None)
micdb=pd.read_csv('/PATH_TO_TOOLS/tools/MiCoP/data/accession2info-fungi.txt',sep='\t',header=None)

resfol='FULL_PATH_TO/Mock_fungal/data/all_methods_summary' # NB! CHANGE TO LOCAL PATH

suffix={'profiles_equal_cov':'EqualCoverage','profiles':'EqualReads'}
prof=['small_1','small_2','small_3','median_1','median_2','median_3','large_1','large_2','complete_1']

def read_proftab(profname):
    fulltab=pd.DataFrame()
    for s in suffix:
        p='' if s=='profiles' else '_EC'
        proftab=pd.read_csv('/'.join([profdir,s,profname+p+'.csv']))
        proftab['ComType']=suffix[s]
        fulltab=pd.concat([fulltab,proftab],ignore_index=False)
    return fulltab

#Find if there are any species/genera/families that are detected in either ER or EC
def find_difference_coverage(tab,profname):
    tools=['Kraken','Euk','Metaphlan','HMS','FunOMIC','MiCoP']
    level=['Species','Genus','Family']

    difference=pd.DataFrame()
    for t in tools:
        for l in level:
            qt=tab[['Taxid','Species',t+'Detected_'+l,'ComType']]
            er=qt.loc[qt['ComType']=='EqualReads']
            er=er.loc[er[t+'Detected_'+l]=='Yes','Species'].tolist()
            ec=qt.loc[qt['ComType']=='EqualCoverage']
            ec=ec.loc[ec[t+'Detected_'+l]=='Yes','Species'].tolist()
            if ec==er:
                print( f'{profname}: {t} and {l} correspond')
            else:
                ec_det=[i for i in ec if i not in er]
                er_det=[i for i in er if i not in ec]
                dif=pd.DataFrame({'Profile': [profname], 'Tool':[t],'Level':[l],'OnlyEC':[ec_det],'OnlyER':[er_det]})
                difference=pd.concat([difference,dif])
            
    return difference

tabdif=pd.DataFrame()
for p in prof:
    tab=read_proftab(p)
    td=find_difference_coverage(tab,p)
    tabdif=pd.concat([tabdif,td])

onlyEC=tabdif[['Profile','Tool','Level','OnlyEC']].explode(column='OnlyEC').dropna(subset='OnlyEC').rename(columns={'OnlyEC':'Species'})
onlyER=tabdif[['Profile','Tool','Level','OnlyER']].explode(column='OnlyER').dropna(subset='OnlyER').rename(columns={'OnlyER':'Species'})

onlyEC.to_csv('/'.join([resfol,'EukDetect_onlyEC_detected_species.csv']),sep='\t',index=False)

#Find genome sizes of the species that are detected only in EC or only in ER
com=read_proftab('complete_1')
onlyEC=onlyEC.merge(com.loc[com['ComType']=='EqualCoverage',['Species','GenomeLength','NumReads']],on='Species',how='left')
onlyER=onlyER.merge(com.loc[com['ComType']=='EqualReads',['Species','GenomeLength','AvgCoverage']],on='Species',how='left')

#Find how many species out of genera were identified in cases where there are more than one species in the genus

tools=['Kraken','Euk','Metaphlan','HMS']
Detection=pd.DataFrame()
for p in prof:
    tab=read_proftab(p)
    for s in suffix:
        query=tab.loc[tab['ComType']==suffix[s]]
        gen=query['Genus name'].value_counts().to_frame()
        #Keep genera with more than one representative
        gen=gen.query('count>1')
        if len(gen)>0:
            for ix,g in gen.iterrows():
                q=query.loc[query['Genus name']==ix]
                inclsp=q['Species name'].tolist()
                for t in tools:
                    if t+'Detected_Species' in q.columns.tolist():
                        dsp=q.loc[q[t+'Detected_Species']=='Yes']
                        numdet=len(dsp)
                        missed=[i for i in inclsp if i not in dsp['Species name'].tolist()]
                        if len(q.loc[q[t+'Detected_Genus']=='Yes'])>0:
                            detgen='Yes'
                        else:
                            detgen='No'
                        d=pd.DataFrame({'Profile':[p],'ComType':[suffix[s]],'Genus':[ix],'NumSpecies':[g['count']],
                                    'InclSpecies':[inclsp],'Tool':[t],'SpeciesDet':[numdet],'MissedSp': [missed], 'GenusDetected':[detgen]})
                        Detection=pd.concat([Detection,d],ignore_index=True)
                        

Detection.to_csv('/'.join([resfol,'Several_Species_inGenus_detection.csv']),index=False)

DetSummary=pd.DataFrame()
for s in suffix:
    det=Detection.loc[Detection['ComType']==suffix[s]] 
    det['cases']=det['Profile']+det['Genus']
    det['PropDet']=det['SpeciesDet']/det['NumSpecies']*100
    num_cases=len(det['cases'].unique().tolist())
    num_genera=len(det['Genus'].unique().tolist())
    for t in tools:
        tdet=det.loc[det['Tool']==t]
        alldet=len(tdet.query('PropDet==100'))
        nonedet=tdet.query('PropDet==0')
        gendet=len(nonedet.loc[nonedet['GenusDetected']=='Yes'])

        ts=pd.DataFrame({'ComType':[suffix[s]],'Tool':[t],'AllDetected':[alldet],'SomeDetected':[len(tdet)-alldet-len(nonedet)],
                         'NoSpeciesButGenusDetected':gendet,
                         'NoSpeciesNoGenusDetected':[len(nonedet)-gendet]})
        DetSummary=pd.concat([DetSummary,ts],ignore_index=True)

#Plot the values
DetSummary=DetSummary.loc[DetSummary['ComType']=='EqualReads']
DetSummary=DetSummary.drop(columns={'ComType'})
fig=DetSummary.plot.barh(stacked=True,color=['#6b519d','#ff66c4','#7ed957','#f4cbb0'],legend=False)
fig.set(yticklabels=['Kraken2','EukDetect','Metaphlan4','HMS'],ylabel='')

#Check which genera are not detected at all and if they are represented in a related database

#Metaphlan
missed_meta=Detection.loc[Detection['Tool']=='Metaphlan']
missed_meta=missed_meta.query('SpeciesDet==0')
missed_meta=missed_meta.loc[missed_meta['GenusDetected']=='No']
missed=pd.DataFrame({'MissedGenus':missed_meta['Genus'].unique().tolist()})

metadb.columns=['SGB','Taxon']
metadb['genus']=metadb['Taxon'].apply(lambda row: row.split('g__')[-1])
metadb['genus']=metadb['genus'].apply(lambda row: row.split('|')[0])   

missed['In_DB']=missed['MissedGenus'].apply(lambda row: 'Yes' if row in metadb['genus'].unique().tolist() else 'No')


#Kraken
missed_kraken=Detection.loc[Detection['Tool']=='Kraken']
missed_kraken=missed_kraken.query('SpeciesDet==0')
missed_kraken=missed_kraken.loc[missed_kraken['GenusDetected']=='No']
missed=pd.DataFrame({'MissedGenus':missed_kraken['Genus'].unique().tolist()})

krdb=krdb.drop(columns='URL')
krdb['genus']=krdb['Sequence Name'].apply(lambda row: ' '.join(row.split(' ')[1:2])) 
krdb=krdb.loc[krdb['#Library']=='fungi']

missed['In_DB']=missed['MissedGenus'].apply(lambda row: 'Yes' if row in krdb['genus'].unique().tolist() else 'No')


#EukDetect
missed_euk=Detection.loc[Detection['Tool']=='Euk']
missed_euk=missed_euk.query('SpeciesDet==0')
missed_euk=missed_euk.loc[missed_euk['GenusDetected']=='No']
missed=pd.DataFrame({'MissedGenus':missed_euk['Genus'].unique().tolist()})

eukdb.columns=['ID','Rec']
eukdb['group']=eukdb['ID'].apply(lambda row: row.split('-')[0])
eukdb['name']=eukdb['ID'].apply(lambda row: row.split('-')[1] if '-' in row else row)
eukdb['genus']=eukdb['name'].apply(lambda row: row.split('_')[0])

missed['In_DB']=missed['MissedGenus'].apply(lambda row: 'Yes' if row in eukdb['genus'].unique().tolist() else 'No')

kaz=eukdb.loc[eukdb['genus']=='Kazachstania']
allkaz=kaz['name'].unique().tolist()


#HMS
missed_hms=Detection.loc[Detection['Tool']=='HMS']
missed_hms=missed_hms.loc[missed_hms['ComType']=='EqualReads']
missed_hms=missed_hms.query('SpeciesDet==0')
missed_hms=missed_hms.loc[missed_hms['GenusDetected']=='No']
missed=pd.DataFrame({'MissedGenus':missed_hms['Genus'].unique().tolist()})

missed['In_DB']=missed['MissedGenus'].apply(lambda row: 'Yes' if row in hmsdb['genus'].unique().tolist() else 'No')

q=hmsdb.loc[hmsdb['genus']=='Blastomyces']
q1=q['species'].unique().tolist()

q=hmsdb.loc[hmsdb['genus']=='Trichophyton']
q1=q['species'].unique().tolist()



#FunOMIC
missed_fun=Detection.loc[Detection['Tool']=='FunOMIC']
missed_fun=missed_fun.loc[missed_fun['ComType']=='EqualReads']
missed_fun=missed_fun.query('SpeciesDet==0')
missed_fun=missed_fun.loc[missed_fun['GenusDetected']=='No']
missed=pd.DataFrame({'MissedGenus':missed_fun['Genus'].unique().tolist()})

missed['In_DB']=missed['MissedGenus'].apply(lambda row: 'Yes' if row in hmsdb['genus'].unique().tolist() else 'No')

#MiCoP
missed_mic=Detection.loc[Detection['Tool']=='MiCoP']
missed_mic=missed_mic.loc[missed_mic['ComType']=='EqualReads']
missed_mic=missed_mic.query('SpeciesDet==0')
missed_mic=missed_mic.loc[missed_mic['GenusDetected']=='No']
missed=pd.DataFrame({'MissedGenus':missed_mic['Genus'].unique().tolist()})

micdb.columns=['ID','Num','Taxid','Taxonomy']
micdb['Taxonomy']=micdb['Taxonomy'].astype(str)
micdb['species']=micdb['Taxonomy'].apply(lambda row: row.split('|')[-2] if 'nan' not in row else None)
micdb=micdb.dropna(subset='species')
missed['In_DB']=missed['MissedGenus'].apply(lambda row: 'Yes' if row in hmsdb['genus'].unique().tolist() else 'No')

q=micdb.loc[micdb['species'].str.contains('Fusarium')]
q1=q['species'].unique().tolist()
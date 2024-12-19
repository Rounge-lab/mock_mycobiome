#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 16:12:32 2024

Check if species is included in the databases

@author: ec-ekateria
"""

import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

profdir='FULL_PATH_TO/Mock_fungal/mock_profiles' # NB! CHANGE TO LOCAL PATH
krdb=pd.read_csv('/PATH_TO_DATABASES/db/kraken-PlusPF/library_report.tsv',sep='\t') # NB! CHANGE TO LOCAL PATH
metadb=pd.read_csv('/PATH_TO_DATABASES/db/biobakery/jun24/metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202307_species.txt',sep='\t',header=None) # NB! CHANGE TO LOCAL PATH
hmsdb=pd.read_csv('/PATH_TO_HMS/tools/MScan2.0/HMS/var/HMS_taxonomy.txt',sep='\t') # NB! CHANGE TO LOCAL PATH
eukdb=pd.read_csv('/PATH_TO_EUKDETECT/EukDetect/eukdb/busco_taxid_link.txt',sep='\t',header=None) # NB! CHANGE TO LOCAL PATH
fundb=pd.read_csv('/PATH_TO_DATABASES/db/FunOMIC/FunOMIC-T/FunOMIC-Tv1/taxonomy.txt',sep='\t',header=None)
micdb=pd.read_csv('/PATH_TO_MICOP/MiCoP/data/accession2info-fungi.txt',sep='\t',header=None)

taxonomy=pd.read_csv('/'.join([profdir,'metadata_ufcg.csv']),sep='\t')
taxonomy[['domain','phylum','class','order','family','genus','species']] = taxonomy['taxonomy'].str.split(';', expand=True)

dbsum=pd.DataFrame({'TotalSpecies':[],'Fungi':[],'MockSpecies':[]})

krdb=krdb.drop(columns='URL')
krdb['Species']=krdb['Sequence Name'].apply(lambda row: ' '.join(row.split(' ')[1:3]))
dbsum.at['Kraken','TotalSpecies']=len(krdb['Species'].unique().tolist())
krdb=krdb.loc[krdb['#Library']=='fungi']
krspecies=krdb['Species'].unique().tolist()
dbsum.at['Kraken','Fungi']=len(krspecies)
taxonomy['In_KrakenDb']=taxonomy['species'].apply(lambda row: 'Yes' if row in krspecies else 'No')
dbsum.at['Kraken','MockSpecies']=len(taxonomy.loc[taxonomy['In_KrakenDb']=='Yes'])

metadb.columns=['SGB','Taxon']
metadb['species']=metadb['Taxon'].apply(lambda row: row.split('s__')[-1])
metadb['species']=metadb['species'].apply(lambda row: row.replace('_',' '))
dbsum.at['Metaphlan','TotalSpecies']=len(metadb['species'].unique().tolist())
metadb=metadb.loc[metadb['Taxon'].str.contains('Euk')]
dbsum.at['Metaphlan','Fungi']=len(metadb.loc[metadb['Taxon'].str.contains('Basidiomycota')])+len(metadb.loc[metadb['Taxon'].str.contains('Ascomycota')])
metaspecies=metadb['species'].unique().tolist()
taxonomy['In_MetaDb']=taxonomy['species'].apply(lambda row: 'Yes' if row in metaspecies else 'No')
dbsum.at['Metaphlan','MockSpecies']=len(taxonomy.loc[taxonomy['In_MetaDb']=='Yes'])

hmsdb['species']=hmsdb['species'].apply(lambda row: row.replace('_',' '))
hmsdb['species']=hmsdb['species'].apply(lambda row: row.replace('[',''))
hmsdb['species']=hmsdb['species'].apply(lambda row: row.replace(']',''))
dbsum.at['HMS','TotalSpecies']=len(hmsdb['species'].unique().tolist())
dbsum.at['HMS','Fungi']=len(hmsdb['species'].unique().tolist())
hmsspecies=hmsdb['species'].unique().tolist()
taxonomy['In_HMSDb']=taxonomy['species'].apply(lambda row: 'Yes' if row in hmsspecies else 'No')
dbsum.at['HMS','MockSpecies']=len(taxonomy.loc[taxonomy['In_HMSDb']=='Yes'])

eukdb.columns=['ID','Rec']
eukdb['group']=eukdb['ID'].apply(lambda row: row.split('-')[0])
eukdb['name']=eukdb['ID'].apply(lambda row: row.split('-')[1] if '-' in row else row)
eukdb['species']=eukdb['name'].apply(lambda row: ' '.join(row.split('_')[0:2]))
dbsum.at['EukDetect','TotalSpecies']=len(eukdb['species'].unique().tolist())
eukdb=eukdb.loc[eukdb['group']=='fungi']
eukspecies=eukdb['species'].unique().tolist()
dbsum.at['EukDetect','Fungi']=len(eukspecies)
taxonomy['In_EukDb']=taxonomy['species'].apply(lambda row: 'Yes' if row in eukspecies else 'No')
dbsum.at['EukDetect','MockSpecies']=len(taxonomy.loc[taxonomy['In_EukDb']=='Yes'])

fundb.columns=['ID','Taxonomy']
fundb['species']=fundb['Taxonomy'].apply(lambda row: row.split('s__')[-1].split('|')[0])
fundb['species']=fundb['species'].apply(lambda row: row.replace('[',''))
fundb['species']=fundb['species'].apply(lambda row: row.replace(']',''))
funspecies=fundb['species'].unique().tolist()
dbsum.at['Funomic','TotalSpecies']=len(funspecies)
dbsum.at['Funomic','Fungi']=len(funspecies)
taxonomy['In_FunDb']=taxonomy['species'].apply(lambda row: 'Yes' if row in funspecies else 'No')
dbsum.at['Funomic','MockSpecies']=len(taxonomy.loc[taxonomy['In_FunDb']=='Yes'])

micdb.columns=['ID','Num','Taxid','Taxonomy']
micdb['Taxonomy']=micdb['Taxonomy'].astype(str)
micdb['species']=micdb['Taxonomy'].apply(lambda row: row.split('|')[-2] if 'nan' not in row else None)
micdb=micdb.dropna(subset='species')
micspecies=micdb['species'].unique().tolist()
dbsum.at['Micop','TotalSpecies']=len(micspecies)
dbsum.at['Micop','Fungi']=len(micspecies)
taxonomy['In_MicDb']=taxonomy['species'].apply(lambda row: 'Yes' if row in micspecies else 'No')
dbsum.at['Micop','MockSpecies']=len(taxonomy.loc[taxonomy['In_MicDb']=='Yes'])


taxonomy[['accession','label','In_KrakenDb','In_MetaDb','In_HMSDb','In_EukDb','In_FunDb','In_MicDb']].to_csv('/'.join([profdir,'metadata_Db_presence.csv']),sep='\t',index=False)
taxonomy.loc[taxonomy['genus'].str.contains('Encephalitozoon'),'class']='Microsporea'

#Make a summary plot
fig, (ax1, ax2,ax3) = plt.subplots(1, 3)
classsum=taxonomy['class'].value_counts().to_frame().reset_index()
sb.barplot(data=classsum,x='count',y='class',color='#6b519d',ax=ax1)
ax1.set(ylabel='',xlabel='Number of species')
sb.boxplot(data=taxonomy, x='GenomeLength',y='class',color='#6b519d',ax=ax2,order=classsum['class'],log_scale=True)
ax2.set(yticklabels='',ylabel='',xlabel='Genome length, bp', xlim=[10**6,10**8])
sb.boxplot(data=taxonomy, x='PercCoreGenes',y='class',color='#6b519d',ax=ax3,order=classsum['class'])
ax3.set(yticklabels='',ylabel='',xlabel='Universal core genes, %')

#Make stacked barplot of database structure
dbsum_norm = dbsum.div(dbsum.sum(axis=1), axis=0) * 100
ax = dbsum_norm.plot(kind='barh', stacked=True, figsize=(7, 4), color=['#6b519d','#ff66c4','#7ed957'])
ax.legend().remove()


####----------------------------------------------------------------------
##Perform nomenclature normalization and compare the lists
import os
from itertools import combinations

os.mkdir('/'.join([profdir,'nomenclature_normalized']))
taxonomy[['species']].to_csv('/'.join([profdir,'nomenclature_normalized','mock_fungi_original']),sep='\t',index=False,header=False)
pd.DataFrame(hmsspecies,columns=['species']).to_csv('/'.join([profdir,'nomenclature_normalized','HMSdb_original']),sep='\t',index=False,header=False)
pd.DataFrame(funspecies,columns=['species']).to_csv('/'.join([profdir,'nomenclature_normalized','FunOMICdb_original']),sep='\t',index=False,header=False)
pd.DataFrame(eukspecies,columns=['species']).to_csv('/'.join([profdir,'nomenclature_normalized','EukDetectdb_original']),sep='\t',index=False,header=False)
pd.DataFrame(krspecies,columns=['species']).to_csv('/'.join([profdir,'nomenclature_normalized','Krakendb_original']),sep='\t',index=False,header=False)
pd.DataFrame(metaspecies,columns=['species']).to_csv('/'.join([profdir,'nomenclature_normalized','MetaPhlAndb_original']),sep='\t',index=False,header=False)
pd.DataFrame(micspecies,columns=['species']).to_csv('/'.join([profdir,'nomenclature_normalized','MiCoPdb_original']),sep='\t',index=False,header=False)

def correct_list(norm):
    names=norm.columns.tolist()
    names=[n.replace('\t','') for n in names]
    norm.columns=names
    norm=norm.rename(columns={'preferred name':'preferred_name', 'name':'TaxID'})
    norm=norm.applymap(lambda x: x.replace('\t', '') if isinstance(x, str) else x)
    norm['combined_list']=norm.apply(lambda row: row.preferred_name if len(row.preferred_name)>1 else str(row.TaxID), axis=1)
    final=norm['combined_list'].tolist()
    return final, norm

norm_taxonomy=pd.read_csv('/'.join([profdir,'nomenclature_normalized','mock_fungi_normalized.txt']),sep='|')
norm_mock,mock_df=correct_list(norm_taxonomy)
mock_df.to_csv('/'.join([profdir,'nomenclature_normalized','mock_fungi_normalized.csv']),sep='\t',index=False)

#Save updated species names list
tax=pd.read_csv('/projects/ec34/katya/projects/Find_fungi_genomes/final_genomes_summary.csv')
tax=tax.merge(mock_df[['TaxID','combined_list']],left_on='Species name',right_on='TaxID')
tax=tax.drop(columns=['Species name','TaxID'])
tax=tax.rename(columns={'combined_list':'Species name'})
tax.to_csv('/projects/ec34/katya/projects/Find_fungi_genomes/final_genomes_summary.csv',sep='\t',index=False)

dbs=['Kraken','MetaPhlAn','HMS','FunOMIC','MiCoP','EukDetect']

norm_dbs=dict()
norm_frames=dict()
for d in dbs:
    dbn=pd.read_csv('/'.join([profdir,'nomenclature_normalized',f'{d}db_normalized.txt']),sep='|')
    n,nf=correct_list(dbn)
    norm_dbs[d]=n
    norm_frames[d]=nf
    nf.to_csv('/'.join([profdir,'nomenclature_normalized',f'{d}db_normalized.csv']),sep='\t',index=False)

#Find number of species which identifier was not recognized by NCBI
norm_summary=pd.DataFrame()
for d in dbs:
    nf=norm_frames[d]
    count=pd.DataFrame(nf['code'].value_counts()).T
    if d=='FunOMIC':
        count=count.reset_index().rename(columns={'index':'Database','1':'Primary names','1+':'PrimaryDuplicated','2':'Secondary names','2+':'SecondaryDuplicated','3':'Name not found'})
    else:
        count=count.reset_index().rename(columns={'index':'Database',1:'Primary names',2:'Secondary names',3:'Name not found'})
    count.at[0,'Database']=d
    count.at[0,'Total']=len(nf)
    norm_summary=pd.concat([norm_summary,count])
norm_summary.to_csv('/'.join([profdir,'nomenclature_normalized','Databases_summary.csv']),sep='\t',index=False)
    
#Find intersect between mock and databases
for d in dbs:
    common=list(set(norm_mock)&set(norm_dbs[d]))
    print(f'Number of species in common between mock community and {d}: {len(common)}')
    
#Find intersect between databases
comb=list(combinations(dbs,2))

inter_len=[]
inter_df=pd.DataFrame()
for c in comb:
    inter_=list(set(norm_dbs[c[0]])&set(norm_dbs[c[1]]))
    inter_len.append(len(inter_))
    print(f'In common between {c[0]} and {c[1]}: {len(inter_)}')
    idf=pd.DataFrame({"Database1":[c[0]],'Database2':[c[1]],'Intersect':[len(inter_)]})
    inter_df=pd.concat([inter_df,idf])
                    
inter_df.to_csv('/'.join([profdir,'nomenclature_normalized','databases_intersect_pairwise_norm.csv']),sep='\t',index=False)
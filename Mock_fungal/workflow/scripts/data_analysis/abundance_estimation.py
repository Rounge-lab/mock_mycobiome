#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 09:52:51 2024

Summarize abundance predictions

@author: ec-ekateria
"""

import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error
import seaborn as sb
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

profdir='FULL_PATH_TO/Mock_fungal/mock_profiles' # NB! CHANGE TO LOCAL PATH
preddir='FULL_PATH_TO/Mock_fungal/data' # NB! CHANGE TO LOCAL PATH
tools={'kraken':'report','metaphlan':'metarep','eukdetect':'filtered_hits_table','HMS':''}
prof=['small_1','small_2','small_3','median_1','median_2','median_3','large_1','large_2','complete_1']

def read_files(p,mod,tool):
   
    global tools,cols,level

    if '_equal_cov' in mod:
        target=pd.read_csv('/'.join([profdir,mod,p+'_EC.csv']))
        if tool !='metaphlan':
            if tool=='kraken':
                rep=pd.read_csv('/'.join([preddir,tool,'equal_coverage_metagenomes',p+'_EC_'+tools[tool]+'.txt']),sep='\t',header=None)
            elif tool=='HMS':
                if 'small' in p or 'median' in p:
                    rep_sp=pd.read_csv('/'.join([preddir,tool,'equal_coverage_metagenomes',p+'_EC',p+'_EC_Species_level_results-%.txt']),sep='\t').reset_index().rename(columns={'index':'Taxon','Unnamed: 0':'count'})
                    rep_g=pd.read_csv('/'.join([preddir,tool,'equal_coverage_metagenomes',p+'_EC',p+'_EC_Genera_level_results-%.txt']),sep='\t').reset_index().rename(columns={'index':'Taxon','V1':'count'})
                    rep_f=pd.read_csv('/'.join([preddir,tool,'equal_coverage_metagenomes',p+'_EC',p+'_EC_Family_level_results-Counts.txt']),sep='\t',header=None).rename(columns={0:'Taxon',1:'count'})
                
                    rep_f['count']=rep_f['count']/rep_f['count'].sum()*100
                    rep={'Species':[rep_sp],'Genus':[rep_g],'Family':[rep_f]}
                else:
                    rep=[]
            elif tool=='funomic':
                rep=pd.read_csv('/'.join([preddir,tool,'equal_coverage_metagenomes',p+'_EC','taxonomic_profiling',p+'_EC_buglist_stratified.txt']),sep='\t')
            elif tool=='micop':
                rep=pd.read_csv('/'.join([preddir,tool,'equal_coverage_metagenomes',p+'_EC_abundance.txt']),sep='\t')
            else: #if EukDetect
                rep=pd.read_csv('/'.join([preddir,tool,'equal_coverage_metagenomes',p+'_EC_'+tools[tool]+'.txt']),sep='\t')
        else: #if metaphlan
            rep=pd.read_csv('/'.join([preddir,tool,'equal_coverage_metagenomes',p+'_EC_'+tools[tool]+'.txt']),sep='\t',skiprows=4)

    else:
        target=pd.read_csv('/'.join([profdir,mod,p+'.csv']))
        if tool !='metaphlan':
            if tool=='kraken':
                rep=pd.read_csv('/'.join([preddir,tool,'equal_reads_metagenomes',p+'_'+tools[tool]+'.txt']),sep='\t',header=None)
            elif tool=='HMS':
                rep_sp=pd.read_csv('/'.join([preddir,tool,'equal_reads_metagenomes',p,p+'_Species_level_results-%.txt']),sep='\t').reset_index().rename(columns={'index':'Taxon','Unnamed: 0':'count'})
                rep_g=pd.read_csv('/'.join([preddir,tool,'equal_reads_metagenomes',p,p+'_Genera_level_results-%.txt']),sep='\t').reset_index().rename(columns={'index':'Taxon','V1':'count'})
                rep_f=pd.read_csv('/'.join([preddir,tool,'equal_reads_metagenomes',p,p+'_Family_level_results-Counts.txt']),sep='\t',header=None).rename(columns={0:'Taxon',1:'count'})
                
                rep_f['count']=rep_f['count']/rep_f['count'].sum()*100
                rep={'Species':[rep_sp],'Genus':[rep_g],'Family':[rep_f]}
            elif tool=='funomic':
                rep=pd.read_csv('/'.join([preddir,tool,'equal_reads_metagenomes',p,'taxonomic_profiling',p+'_buglist_stratified.txt']),sep='\t')
            elif tool=='micop':
                rep=pd.read_csv('/'.join([preddir,tool,'equal_reads_metagenomes',p+'_abundance.txt']),sep='\t')
            else: #if EukDetect
                rep=pd.read_csv('/'.join([preddir,tool,'equal_reads_metagenomes',p+'_'+tools[tool]+'.txt']),sep='\t')           
        else: #if Metaphlan
            rep=pd.read_csv('/'.join([preddir,tool,'equal_reads_metagenomes',p+'_'+tools[tool]+'.txt']), skiprows=4,sep='\t')

    if tool=='kraken':
        rep=rep.rename(columns={0:'Taxonomy',1: 'NumReads'})

        rep['Species name']=rep['Taxonomy'].apply(lambda row: row.split('s__')[-1] if 's__' in row else None)
        rep['Genus name']=rep['Taxonomy'].apply(lambda row: row.split('g__')[-1].split('|s__')[0] if 'g__' in row else None)
        rep['Family name']=rep['Taxonomy'].apply(lambda row: row.split('f__')[-1].split('|g__')[0] if 'f__' in row else None)
    
    elif tool=='eukdetect':
        rep['Species name']=rep['Lineage'].apply(lambda row: row.split('|species-')[-1] if '|species-' in row else None)
        rep['Species name']=rep['Species name'].str.replace('_',' ')
        rep['Genus name']=rep['Lineage'].apply(lambda row: row.split('|genus-')[-1].split('|species-')[0] if '|genus-' in row else None)
        rep['Family name']=rep['Lineage'].apply(lambda row: row.split('|family-')[-1].split('|genus-')[0] if '|family-' in row else None)
        
        rep['PredictedRelAb']=rep['Read_counts']/rep['Read_counts'].sum()*100
    
    elif tool=='metaphlan':
        rep=rep.rename(columns={'#clade_name':'Taxonomy','relative_abundance': 'RelAb'})
        rep=rep.drop(columns={'NCBI_tax_id','additional_species'})
        
        rep['Species name']=rep['Taxonomy'].apply(lambda row: row.split('s__')[-1].split('|t__')[0] if 's__' in row else None)
        rep['Species name']=rep['Species name'].apply(lambda row: ' '.join(row.split('_')[0:2]) if row is not None else row)
        rep['Genus name']=rep['Taxonomy'].apply(lambda row: row.split('g__')[-1].split('|s__')[0] if 'g__' in row else None)
        rep['Family name']=rep['Taxonomy'].apply(lambda row: row.split('f__')[-1].split('|g__')[0] if 'f__' in row else None)
    
    elif tool=='funomic':
        rep=rep.rename(columns={'taxa':'Taxonomy','Abundance': 'NumReads'})

        rep['Species name']=rep['Taxonomy'].apply(lambda row: row.split('s__')[-1] if 's__' in row else None)
        rep['Genus name']=rep['Taxonomy'].apply(lambda row: row.split('g__')[-1].split('|s__')[0] if 'g__' in row else None)
        rep['Family name']=rep['Taxonomy'].apply(lambda row: row.split('f__')[-1].split('|g__')[0] if 'f__' in row else None)
       
    elif tool=='micop':
        rep=rep[['RANK','TAXPATHSN','PERCENTAGE']]
        rep['TAXPATHSN']=rep['TAXPATHSN'].apply(lambda row: row.split('|')[-1])
        rep['RANK']=rep['RANK'].apply(lambda row: row.capitalize())
        rep=rep.loc[rep['RANK'].isin(level)]
        
    target['RelAbReads']=target['NumReads']/target['NumReads'].sum()*100
    target['RelAbCoverage']=target['AvgCoverage']/target['AvgCoverage'].sum()*100
  
    return target, rep


#difference between true and predicted number of reads
def kraken_error(tar,kraken):

    report=pd.DataFrame()
    all_pred=pd.DataFrame()
    for l in level:
        pq=kraken.dropna(subset=l+' name')
        if l!='Species':
            if l=='Genus':
                pq=pq.loc[pq['Species name'].isnull()]
            else:
                pq=pq.loc[pq['Genus name'].isnull()]
        pq['PredictedRelAb']=pq['NumReads']/pq['NumReads'].sum()*100
        tq=tar.loc[tar[cols['kraken']+l]=='Yes']
        pq=pq.loc[pq[l+' name'].isin(tq[l+' name'].tolist())]

        tq=tq.loc[tq[l+' name'].isin(pq[l+' name'].unique().tolist())]
        tq=tq[[l+' name','RelAbReads']]
        tq = tq.groupby(l+' name', as_index=False).sum()
        tq=tq.merge(pq[[l+' name', 'PredictedRelAb']],on=l+' name', how='left')
        
        mae=mean_absolute_error(tq['RelAbReads'],tq['PredictedRelAb'])
        rmse=mean_squared_error(tq['RelAbReads'],tq['PredictedRelAb'],squared=False)
        
        r=pd.DataFrame({'Mock_community':[p],'Tool':['kraken'],'Level':[l],'AvgRelAb':[tq['RelAbReads'].mean().astype(int)],'MAE':[mae.astype(int)],'RMSE':[rmse.astype(int)]})
        report=pd.concat([report,r],ignore_index=True)
        
        tq=tq.rename(columns={l+' name':'Taxon name'})
        tq['Mock_community']=p
        tq['Tool']='kraken'
        tq['Level']=l
        
        all_pred=pd.concat([all_pred,tq],ignore_index=True)
        
    return report,all_pred


#difference between true (based on genome coverage) and predicted relative abundance 
def eukdetect_error(tar,euk):
    
    report=pd.DataFrame()
    all_pred=pd.DataFrame()

    for l in level:
        tq=tar.loc[tar[cols['eukdetect']+l]=='Yes']
        pq=euk.loc[euk[l+' name'].isin(tq[l+' name'].tolist())]
        pq=pq[[l+' name','PredictedRelAb']]
        pq = pq.groupby(l+' name', as_index=False).sum()
        tq=tq.loc[tq[l+' name'].isin(pq[l+' name'].unique().tolist())]
        tq=tq[[l+' name','RelAbCoverage']]
        tq = tq.groupby(l+' name', as_index=False).sum()
        tq=tq.merge(pq[[l+' name', 'PredictedRelAb']],on=l+' name', how='left')
        mae=mean_absolute_error(tq['RelAbCoverage'],tq['PredictedRelAb'])
        rmse=mean_squared_error(tq['RelAbCoverage'],tq['PredictedRelAb'],squared=False)
        r=pd.DataFrame({'Mock_community':[p],'Tool':['eukdetect'],'Level':[l],'AvgRelAb':[tq['RelAbCoverage'].mean().astype(int)],'MAE':[mae.astype(int)],'RMSE':[rmse.astype(int)]})
        report=pd.concat([report,r])
        
        tq=tq.rename(columns={l+' name':'Taxon name'})
        tq['Mock_community']=p
        tq['Tool']='eukdetect'
        tq['Level']=l
        
        all_pred=pd.concat([all_pred,tq],ignore_index=True)
        
    return report,all_pred

     
def metaphlan_error(tar,meta):
    
    report=pd.DataFrame()
    all_pred=pd.DataFrame()

    for l in level:
        tq=tar.loc[tar[cols['metaphlan']+l]=='Yes']
        pq=meta.loc[meta[l+' name'].isin(tq[l+' name'].tolist())]       
        if l!='Species':
            if l=='Genus':
                pq=pq.loc[pq['Species name'].isnull()]
            else:
                pq=pq.loc[pq['Genus name'].isnull()]
        else:
            pq=pq.drop_duplicates(subset='Species name')
        pq=pq.rename(columns={'RelAb':'PredictedRelAb'})
        pq=pq[[l+' name','PredictedRelAb']]
        tq=tq.loc[tq[l+' name'].isin(pq[l+' name'].unique().tolist())]
        tq=tq[[l+' name','RelAbCoverage']]
        tq = tq.groupby(l+' name', as_index=False).sum()
        tq=tq.merge(pq[[l+' name', 'PredictedRelAb']],on=l+' name', how='left')
        mae=mean_absolute_error(tq['RelAbCoverage'],tq['PredictedRelAb'])
        rmse=mean_squared_error(tq['RelAbCoverage'],tq['PredictedRelAb'],squared=False)
        r=pd.DataFrame({'Mock_community':[p],'Tool':['metaphlan'],'Level':[l],'AvgRelAb':[tq['RelAbCoverage'].mean().astype(int)],'MAE':[mae.astype(int)],'RMSE':[rmse.astype(int)]})
        report=pd.concat([report,r])
        
        tq=tq.rename(columns={l+' name':'Taxon name'})
        tq['Mock_community']=p
        tq['Tool']='metaphlan'
        tq['Level']=l
        
        all_pred=pd.concat([all_pred,tq],ignore_index=True)
        
    return report,all_pred

def hms_error(tar,hms):
    
    report=pd.DataFrame()
    all_pred=pd.DataFrame()

    for l in level:
        tq=tar.loc[tar[cols['HMS']+l]=='Yes']
        pq=hms[l][0]     
        if l=='Species':
            pq['Taxon']=pq['Taxon'].apply(lambda row: row.replace('_',' '))
       
        pq=pq.rename(columns={'count':'PredictedRelAb'})
       
        tq=tq.loc[tq[l+' name'].isin(pq['Taxon'].unique().tolist())]
        tq=tq[[l+' name','RelAbReads']]
        tq = tq.groupby(l+' name', as_index=False).sum()
        tq=tq.merge(pq[['Taxon', 'PredictedRelAb']],left_on=l+' name',right_on='Taxon', how='left')
        tq=tq.drop(columns='Taxon')
        mae=mean_absolute_error(tq['RelAbReads'],tq['PredictedRelAb'])
        rmse=mean_squared_error(tq['RelAbReads'],tq['PredictedRelAb'],squared=False)
        r=pd.DataFrame({'Mock_community':[p],'Tool':['HMS'],'Level':[l],'AvgRelAb':[tq['RelAbReads'].mean().astype(int)],'MAE':[mae.astype(int)],'RMSE':[rmse.astype(int)]})
        report=pd.concat([report,r])
        
        tq=tq.rename(columns={l+' name':'Taxon name'})
        tq['Mock_community']=p
        tq['Tool']='HMS'
        tq['Level']=l
        
        all_pred=pd.concat([all_pred,tq],ignore_index=True)
        
    return report,all_pred

def funomic_error(tar,fun):

    report=pd.DataFrame()
    all_pred=pd.DataFrame()
    for l in level:
        pq=fun.dropna(subset=l+' name')
        if l!='Species':
            if l=='Genus':
                pq=pq.loc[pq['Species name'].isnull()]
            else:
                pq=pq.loc[pq['Genus name'].isnull()]
        pq['PredictedRelAb']=pq['NumReads']/pq['NumReads'].sum()*100
        tq=tar.loc[tar[cols['funomic']+l]=='Yes']
        pq=pq.loc[pq[l+' name'].isin(tq[l+' name'].tolist())]

        tq=tq.loc[tq[l+' name'].isin(pq[l+' name'].unique().tolist())]
        tq=tq[[l+' name','RelAbCoverage']]
        tq = tq.groupby(l+' name', as_index=False).sum()
        tq=tq.merge(pq[[l+' name', 'PredictedRelAb']],on=l+' name', how='left')
        
        mae=mean_absolute_error(tq['RelAbCoverage'],tq['PredictedRelAb'])
        rmse=mean_squared_error(tq['RelAbCoverage'],tq['PredictedRelAb'],squared=False)
        
        r=pd.DataFrame({'Mock_community':[p],'Tool':['funomic'],'Level':[l],'AvgRelAb':[tq['RelAbCoverage'].mean().astype(int)],'MAE':[mae.astype(int)],'RMSE':[rmse.astype(int)]})
        report=pd.concat([report,r],ignore_index=True)
        
        tq=tq.rename(columns={l+' name':'Taxon name'})
        tq['Mock_community']=p
        tq['Tool']='funomic'
        tq['Level']=l
        
        all_pred=pd.concat([all_pred,tq],ignore_index=True)
        
    return report,all_pred

def micop_error(tar,mic):
    
    report=pd.DataFrame()
    all_pred=pd.DataFrame()

    for l in level:
        tq=tar.loc[tar[cols['micop']+l]=='Yes']
        
        pq=mic.loc[mic['RANK']==l]
        tq=tq.loc[tq[l+' name'].isin(pq['TAXPATHSN'].unique().tolist())]
        tq=tq[[l+' name','RelAbReads']]
        tq = tq.groupby(l+' name', as_index=False).sum()
        tq=tq.merge(pq[['TAXPATHSN', 'PERCENTAGE']],left_on=l+' name',right_on='TAXPATHSN', how='left')
        tq=tq.drop(columns='TAXPATHSN')
        mae=mean_absolute_error(tq['RelAbReads'],tq['PERCENTAGE'])
        rmse=mean_squared_error(tq['RelAbReads'],tq['PERCENTAGE'],squared=False)
        r=pd.DataFrame({'Mock_community':[p],'Tool':['micop'],'Level':[l],'AvgRelAb':[tq['RelAbReads'].mean().astype(int)],'MAE':[mae.astype(int)],'RMSE':[rmse.astype(int)]})
        report=pd.concat([report,r])
        
        tq=tq.rename(columns={l+' name':'Taxon name'})
        tq['Mock_community']=p
        tq['Tool']='micop'
        tq['Level']=l
        
        all_pred=pd.concat([all_pred,tq],ignore_index=True)
        
    return report,all_pred
#Read a profile and detected reads

mod={'profiles':'EqualReads','profiles_equal_cov':'EqualCoverage'}

KrakenReport=pd.DataFrame()
EukReport=pd.DataFrame()
MetaReport=pd.DataFrame()
HMSReport=pd.DataFrame()
FunomicReport=pd.DataFrame()
MicopReport=pd.DataFrame()


KrakenPred=pd.DataFrame()
EukPred=pd.DataFrame()
MetaPred=pd.DataFrame()
HMSPred=pd.DataFrame()
FunomicPred=pd.DataFrame()
MicopPred=pd.DataFrame()



for m in mod:
    print(m)
    for p in prof:
        tar,euk=read_files(p,m,'eukdetect')
        _,meta=read_files(p,m,'metaphlan')
        _,kraken=read_files(p,m,'kraken')
        _,hms=read_files(p,m,'HMS')
        _,fun=read_files(p,m,'funomic')
        _,mic=read_files(p,m,'micop')
    
        #Assess predictions errors
        krep,krpredictions=kraken_error(tar,kraken)
        eukrep,eukpredictions=eukdetect_error(tar,euk)
        metarep,metapredictions=metaphlan_error(tar,meta) 
        if len(hms)>0:
            hmsrep,hmspredictions=hms_error(tar,hms)
        funrep,funpredictions=funomic_error(tar,fun)
        micrep,micpredictions=micop_error(tar,mic)


        krep['ComType']=mod[m]
        eukrep['ComType']=mod[m]
        metarep['ComType']=mod[m]
        if len(hms)>0:
            hmsrep['ComType']=mod[m]
        funrep['ComType']=mod[m]
        micrep['ComType']=mod[m]

    
        krpredictions['ComType']=mod[m]
        eukpredictions['ComType']=mod[m]
        metapredictions['ComType']=mod[m]
        if len(hms)>0:
            hmspredictions['ComType']=mod[m]
        funpredictions['ComType']=mod[m]
        micpredictions['ComType']=mod[m]


        KrakenReport=pd.concat([KrakenReport,krep])
        EukReport=pd.concat([EukReport,eukrep])
        MetaReport=pd.concat([MetaReport,metarep])
        if len(hms)>0:
            HMSReport=pd.concat([HMSReport,hmsrep])
        FunomicReport=pd.concat([FunomicReport,funrep])
        MicopReport=pd.concat([MicopReport,micrep])

        
        KrakenPred=pd.concat([KrakenPred,krpredictions])
        EukPred=pd.concat([EukPred,eukpredictions])
        MetaPred=pd.concat([MetaPred,metapredictions])
        if len(hms)>0:
            HMSPred=pd.concat([HMSPred,hmspredictions])
        FunomicPred=pd.concat([FunomicPred,funpredictions])
        MicopPred=pd.concat([MicopPred,micpredictions])

##Add Bracken classification for Kraken RelAb calculations
profdir='/fp/homes01/u01/ec-ekateria/ec34/katya/projects/Find_fungi_genomes/fungal_genomes'
preddir='/fp/homes01/u01/ec-ekateria/ec34/katya/projects/Mock_fungal/data'

KrakenPred=pd.read_csv('/'.join([preddir,'all_methods_summary/KrakenRelAbPred.csv']))
BrackenReport=pd.DataFrame()
mod={'profiles':['EqualReads','equal_reads_metagenomes'],'profiles_equal_cov':['EqualCoverage','equal_coverage_metagenomes']}
for m in mod:
    for p in prof:
        if m=='profiles':
            pred=pd.read_csv('/'.join([preddir,'kraken',mod[m][1],p+'_bracken_report.txt']),sep='\t',header=None)
        else:
            pred=pd.read_csv('/'.join([preddir,'kraken',mod[m][1],p+'_EC_bracken_report.txt']),sep='\t',header=None)
        
        pred.columns=['PredictedRelAb','NumReads','NRThisLevel','Rank','TaxID','Taxon name']
        pred['Taxon name']=pred['Taxon name'].str.lstrip()
        
        target=KrakenPred.loc[KrakenPred['ComType']==mod[m][0]]
        target=target.loc[target['Mock_community']==p]
        
        for l in level:
            pr=pred.loc[pred['Rank']==l[0]]
            tar=target.loc[target['Level']==l]
            
            comb=tar[['Taxon name','RelAbReads']].merge(pr[['PredictedRelAb','Taxon name']],on='Taxon name',how='left')
            mae=mean_absolute_error(comb['RelAbReads'],comb['PredictedRelAb'])
            rmse=mean_squared_error(comb['RelAbReads'],comb['PredictedRelAb'],squared=False)
            
            r=pd.DataFrame({'Mock_community':[p],'Tool':['kraken'],'Level':[l],'AvgRelAb':[comb['RelAbReads'].mean().astype(int)],'MAE':[mae.astype(int)],'RMSE':[rmse.astype(int)],'ComType':[mod[m][0]]})
            BrackenReport=pd.concat([BrackenReport,r])
        
BrackenReport.to_csv('/'.join([preddir,'all_methods_summary/KrakenRMSE_bracken.csv']),index=False)
         
        
KrakenReport.to_csv('/'.join([preddir,'all_methods_summary/KrakenRMSE.csv']),index=False)
EukReport.to_csv('/'.join([preddir,'all_methods_summary/EukdetectRMSE.csv']),index=False)
MetaReport.to_csv('/'.join([preddir,'all_methods_summary/MetaphlanRMSE.csv']),index=False)
HMSReport.to_csv('/'.join([preddir,'all_methods_summary/HMSRMSE.csv']),index=False)
FunomicReport.to_csv('/'.join([preddir,'all_methods_summary/FunomicRMSE.csv']),index=False)
MicopReport.to_csv('/'.join([preddir,'all_methods_summary/MicopRMSE.csv']),index=False)

KrakenPred.to_csv('/'.join([preddir,'all_methods_summary/KrakenRelAbPred.csv']),index=False)
EukPred.to_csv('/'.join([preddir,'all_methods_summary/EukdetectRelAbPred.csv']),index=False)
MetaPred.to_csv('/'.join([preddir,'all_methods_summary/MetaphlanRelAbPred.csv']),index=False)
HMSPred.to_csv('/'.join([preddir,'all_methods_summary/HMSRelAbPred.csv']),index=False)
FunomicPred.to_csv('/'.join([preddir,'all_methods_summary/FunomicRelAbPred.csv']),index=False)
MicopPred.to_csv('/'.join([preddir,'all_methods_summary/MicopRelAbPred.csv']),index=False)

#Visualize data
KrakenReport=pd.read_csv('/'.join([preddir,'all_methods_summary/KrakenRMSE.csv']))
EukReport=pd.read_csv('/'.join([preddir,'all_methods_summary/EukdetectRMSE.csv']))
MetaReport=pd.read_csv('/'.join([preddir,'all_methods_summary/MetaphlanRMSE.csv']))
HMSReport=pd.read_csv('/'.join([preddir,'all_methods_summary/HMSRMSE.csv']))
FunomicReport=pd.read_csv('/'.join([preddir,'all_methods_summary/FunomicRMSE.csv']))
MicopReport=pd.read_csv('/'.join([preddir,'all_methods_summary/MicopRMSE.csv']))

# KrakenReport=KrakenReport.rename(columns={'AvgRelAb':'AvgAbunEst'})
# EukReport=EukReport.rename(columns={'AvgRelAb':'AvgAbunEst'})
# MetaReport=MetaReport.rename(columns={'AvgRelAb':'AvgAbunEst'})

Report=pd.concat([KrakenReport,EukReport,MetaReport,HMSReport, FunomicReport,MicopReport],ignore_index=True)

sb.boxplot(data=Report,x='Tool',y='RMSE',palette=['#6b519d','#ff66c4','#7ed957','#f4cbb0'],order=['kraken','metaphlan','HMS','funomic','micop','eukdetect'],hue='Level')


def kruskal_test(data,grcol,val):
    groups = data.groupby(grcol)[val].apply(list)
    hst, pval = kruskal(*groups)
    
    return hst, pval

def kruskal_erec(data):
    KruskalEREC=pd.DataFrame()
    level=['Species','Genus','Family']
    tools=data['Tool'].unique().tolist()
    for l in level:
        query=data.loc[data['Level']==l]
        for t in tools:
            if len(query.loc[query['Tool']==t,'RMSE'].unique())>1:
                hst,pval=kruskal_test(query.loc[query['Tool']==t],'ComType','RMSE')
            else:
                hst=None
                pval=1
            resp=pd.DataFrame({'Tool':[t],'Level':[l],'Hst':[hst],'Pval':[pval]})
            KruskalEREC=pd.concat([KruskalEREC,resp],ignore_index=True)
    _, fdrp, _, _ = multipletests(KruskalEREC['Pval'].values, alpha=0.05, method='fdr_bh')
    KruskalEREC['FDRp']=fdrp
    
    return KruskalEREC

KruskalEREC=kruskal_erec(Report)
KruskalEREC.to_csv('/'.join([preddir,'all_methods_summary/Kruskal_EREC_RMSE.csv']),index=False)

def kruskal_tool(data,col):
    KruskalTool=pd.DataFrame()
    level=['Species','Genus','Family']
    tools=data['Tool'].unique().tolist()
    for l in level:
        query=data.loc[data['Level']==l]
        comb=list(combinations(tools, 2))
        for t in comb:
            hst,pval=kruskal_test(query.loc[query['Tool'].isin(t)],'Tool',col)
            resp=pd.DataFrame({'Tool':['_vs_'.join(list(t))],'Level':[l],'Hst':[hst],'Pval':[pval]})
            KruskalTool=pd.concat([KruskalTool,resp],ignore_index=True)
    _, fdrp, _, _ = multipletests(KruskalTool['Pval'].values, alpha=0.05, method='fdr_bh')
    KruskalTool['FDRp']=fdrp
    
    return KruskalTool

KruskalTool=kruskal_tool(Report,'RMSE')
KruskalTool.to_csv('/'.join([preddir,'all_methods_summary/Kruskal_Tools_RMSE.csv']),index=False)
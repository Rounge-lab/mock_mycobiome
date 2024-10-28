#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 09:46:42 2024

Summarize results

@author: ec-ekateria
"""

import pandas as pd
import seaborn as sb
from scipy.stats import kruskal
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests


wdir='/FULL_PATH_TO/Mock_fungal/data/' # NB! CHANGE TO LOCAL PATH
tool=['kraken','metaphlan','eukdetect','HMS','funomic','micop']
com=['EqualReads','EqualCoverage']
resdir='/FULL_PATH_TO/data/all_methods_summary' # NB! CHANGE TO LOCAL PATH



def get_summaries():
    allsum=pd.DataFrame()
    global wdir,tool,com
    for t in tool:
        for c in com:
            q=pd.read_csv('/'.join([wdir,t,'summaries',c+'_Summary.csv'])).drop(columns='Unnamed: 0')
            q['Tool']=t
            q['ComType']=c
            allsum=pd.concat([allsum,q],ignore_index=True)
            
    return allsum

allsum=get_summaries()

prec=pd.melt(allsum[['NumSpecies','PrecisionSpecies','PrecisionGenera','PrecisionFamily','Tool','ComType']],id_vars=['NumSpecies','Tool','ComType'])
prec['variable']=prec['variable'].str.replace('Precision','')
prec['variable']=prec['variable'].apply(lambda row: 'Genus' if row=='Genera' else row)

recall=pd.melt(allsum[['NumSpecies','RecallSpecies','RecallGenera','RecallFamily','Tool','ComType']],id_vars=['NumSpecies','Tool','ComType'])
recall['variable']=recall['variable'].str.replace('Recall','')
recall['variable']=recall['variable'].apply(lambda row: 'Genus' if row=='Genera' else row)

sb.catplot(data=prec,x='Tool',y='value',hue='ComType',col='variable',kind='box',palette=['#9697b1','#1c6462'],aspect=.6)
sb.catplot(data=recall,x='Tool',y='value',hue='ComType',col='variable',kind='box',palette=['#9697b1','#1c6462'],aspect=.6)


def kruskal_test(data,grcol,val):
    groups = data.groupby(grcol)[val].apply(list)
    hst, pval = kruskal(*groups)
    
    return hst, pval

def kruskal_erec(data,col):
    KruskalEREC=pd.DataFrame()
    level=['Species','Genus','Family']
    for l in level:
        query=data.loc[data['variable']==l]
        for t in tool:
            if len(query.loc[query['Tool']==t,col].unique())>1:
                hst,pval=kruskal_test(query.loc[query['Tool']==t],'ComType',col)
            else:
                hst=None
                pval=1
            resp=pd.DataFrame({'Tool':[t],'Level':[l],'Hst':[hst],'Pval':[pval]})
            KruskalEREC=pd.concat([KruskalEREC,resp],ignore_index=True)
    _, fdrp, _, _ = multipletests(KruskalEREC['Pval'].values, alpha=0.05, method='fdr_bh')
    KruskalEREC['FDRp']=fdrp
    
    return KruskalEREC

PrecKruskalEREC=kruskal_erec(prec,'value')
PrecKruskalEREC.to_csv('/'.join([resdir,'Precision_ERvsEC_Kruskal.csv']),index=False)
    
RecallKruskalEREC=kruskal_erec(recall,'value')
RecallKruskalEREC.to_csv('/'.join([resdir,'Recall_ERvsEC_Kruskal.csv']),index=False)


def kruskal_tool(data,col):
    KruskalTool=pd.DataFrame()
    level=['Species','Genus','Family']
    for l in level:
        query=data.loc[data['variable']==l]
        for t in tool:
            hst,pval=kruskal_test(query.loc[query['Tool']!=t],'Tool',col)
            resp=pd.DataFrame({'Tool':['_vs_'.join([tl for tl in tool if tl!=t])],'Level':[l],'Hst':[hst],'Pval':[pval]})
            KruskalTool=pd.concat([KruskalTool,resp],ignore_index=True)
    _, fdrp, _, _ = multipletests(KruskalTool['Pval'].values, alpha=0.05, method='fdr_bh')
    KruskalTool['FDRp']=fdrp
    
    return KruskalTool
  
PrecKruskalTool=kruskal_tool(prec,'value')
PrecKruskalTool.to_csv('/'.join([resdir,'Precision_Tools_Kruskal.csv']),index=False)

RecallKruskalTool=kruskal_tool(recall,'value')
RecallKruskalTool.to_csv('/'.join([resdir,'Recall_Tools_Kruskal.csv']),index=False)


#Check if there is correlation between number of species in the community and precision/recall of tools

def pearson_cor(data,col):
    Pearson=pd.DataFrame()
    level=['Species','Genus','Family']
    for l in level:
        query=data.loc[data['variable']==l]
        for t in tool:
            if len(query.loc[query['Tool']==t,col].unique())>1:
                r,pval=pearsonr(query.loc[query['Tool']==t,'NumSpecies'],query.loc[query['Tool']==t,col])
            else:
                r=0
                pval=1
            resp=pd.DataFrame({'Tool':[t],'Level':[l],'PearsonR':[r],'PearsonR2':[r**2],'Pval':[pval]})
            Pearson=pd.concat([Pearson,resp],ignore_index=True)
    _, fdrp, _, _ = multipletests(Pearson['Pval'].values, alpha=0.05, method='fdr_bh')
    Pearson['FDRp']=fdrp
    
    return Pearson

PearsonPrec=pearson_cor(prec,'value')
PearsonPrec.to_csv('/'.join([resdir,'Precision_vs_NumSpecies_Pearson.csv']),index=False)

PearsonRecall=pearson_cor(recall,'value')
PearsonRecall.to_csv('/'.join([resdir,'Recall_vs_NumSpecies_Pearson.csv']),index=False)

fig=sb.barplot(data=PearsonPrec,x='PearsonR',y='Tool',hue='Level',palette=['#9697b1','#1c6462','#eac0bf'],edgecolor='gray',width=.6)
fig=sb.barplot(data=PearsonRecall,x='PearsonR',y='Tool',hue='Level',palette=['#9697b1','#1c6462','#eac0bf'],edgecolor='gray',width=.6).legend_.remove()

PearsonPrec['Metric']='Precision'
PearsonRecall['Metric']='Recall'

Pearson=pd.concat([PearsonPrec,PearsonRecall],ignore_index=True)
fig=sb.boxplot(data=Pearson,y='Tool',x='PearsonR',hue='Metric',palette=['#6b519d','#ff66c4'])


#calculate F1 score for the accuracy
prec['group']=prec.apply(lambda row: '_'.join([str(row.NumSpecies), row.Tool, row.ComType, row.variable]), axis=1)
prec=prec.rename(columns={'value':'Precision'})

recall['group']=recall.apply(lambda row: '_'.join([str(row.NumSpecies), row.Tool, row.ComType, row.variable]), axis=1)
recall=recall.rename(columns={'value':'Recall'})

compl=prec.merge(recall[['Recall']], left_index=True, right_index=True, how='left')
compl['F1']=2*(compl['Precision']*compl['Recall'])/(compl['Precision']+compl['Recall'])

fig=sb.boxplot(data=compl,x='Tool',y='F1',hue='variable',palette=['#6b519d','#ff66c4','#7ed957','#f4cbb0'],order=['kraken','metaphlan','HMS','eukdetect','funomic','micop'])
fig.set(ylabel='F1 score',xlabel='')

#Calulate median and IQR
ER=compl.loc[compl['ComType']=='EqualReads',['Tool','variable','F1']].groupby(['Tool', 'variable'])['F1'].median().reset_index()
ER.columns=['Tool','Level','F1_ER']
EC=compl.loc[compl['ComType']=='EqualCoverage',['Tool','variable','F1']].groupby(['Tool', 'variable'])['F1'].median().reset_index()
EC.columns=['Tool','Level','F1_EC']

def calc_iqr(series):
    Q1 = series.quantile(0.25)
    Q3 = series.quantile(0.75)
    IQR = Q3 - Q1
    return IQR

ERIQR=compl.loc[compl['ComType']=='EqualReads',['Tool','variable','F1']].groupby(['Tool', 'variable'])['F1'].agg(calc_iqr).reset_index()
ERIQR.columns=['Tool','Level','IQR_ER']

ECIQR=compl.loc[compl['ComType']=='EqualCoverage',['Tool','variable','F1']].groupby(['Tool', 'variable'])['F1'].agg(calc_iqr).reset_index()
ECIQR.columns=['Tool','Level','IQR_EC']

MedianIQR=ER.merge(EC,on=['Tool','Level'])
MedianIQR=MedianIQR.merge(ERIQR,on=['Tool','Level'])
MedianIQR=MedianIQR.merge(ECIQR,on=['Tool','Level'])

AccKruskalEREC=kruskal_erec(compl,'F1')
AccKruskalEREC=AccKruskalEREC.merge(MedianIQR,on=['Tool','Level'])
AccKruskalEREC=AccKruskalEREC[['Tool','Level','F1_ER','IQR_ER','F1_EC','IQR_EC','Hst','Pval','FDRp']]

AccKruskalEREC.to_csv('/'.join([resdir,'AccuracyF1_ERvsEC_Kruskal.csv']),index=False)

AccKruskalTool=kruskal_tool(compl,'F1')
AccKruskalTool.to_csv('/'.join([resdir,'AccuracyF1_Tools_Kruskal.csv']),index=False)

PearsonAcc=pearson_cor(compl,'F1')
PearsonAcc.to_csv('/'.join([resdir,'F1_vs_NumSpecies_Pearson.csv']),index=False)

#Add data on relative abundance
KrakenReport=pd.read_csv('/'.join([resdir,'KrakenRMSE.csv']))
EukReport=pd.read_csv('/'.join([resdir,'EukdetectRMSE.csv']))
MetaReport=pd.read_csv('/'.join([resdir,'MetaphlanRMSE.csv']))
HMSReport=pd.read_csv('/'.join([resdir,'HMSRMSE.csv']))

#KrakenReport=KrakenReport.rename(columns={'AvgRelAb':'AvgAbunEst'})
#EukReport=EukReport.rename(columns={'AvgRelAb':'AvgAbunEst'})
#MetaReport=MetaReport.rename(columns={'AvgRelAb':'AvgAbunEst'})
#HMSReport=HMSReport.rename(columns={'AvgRelAb':'AvgAbunEst'})

Report=pd.concat([KrakenReport,EukReport,MetaReport,HMSReport],ignore_index=True)
Report=Report.rename(columns={'Level':'variable'})
Report['NumSpecies']=Report['Mock_community'].apply(lambda row: int(10) if 'small' in row else (int(50) if 'median' in row else (int(100) if 'large' in row else int(165))))
PearsonRMSE=pearson_cor(Report,'RMSE')
PearsonRMSE.to_csv('/'.join([resdir,'RMSE_vs_NumSpecies_Pearson.csv']),index=False)

PearsonAcc['Metric']='F1'
PearsonRMSE['Metric']='RMSE'

Pearson=pd.concat([PearsonAcc,PearsonRMSE],ignore_index=True)
fig=sb.boxplot(data=Pearson,y='Tool',x='PearsonR',hue='Metric',palette=['#6b519d','#ff66c4','#7ed957','#f4cbb0'],order=['kraken','metaphlan','HMS','eukdetect','funomic','micop'])

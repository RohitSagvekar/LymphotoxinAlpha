import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import scipy.stats as stats
from scipy import stats
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator
#matplotlib inline
import os


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 150)

#cancers = ['BRCA','CESC','CHOL','COAD','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','PAAD','PRAD','PCPG','READ','SARC','SKCM','STAD','THCA','THYM','UCEC']
#cancers = ['CHOL','KICH','SKCM']
cancers = ['CHOL']
gene_number = 12130 
normalization_technique = 'fpkm_unstranded'

df3 = pd.DataFrame()

for i in cancers:
    print(i)

    df = pd.DataFrame(pd.read_csv('data/mRNA_exp_Normal_RAW_txt/'+i+'_normal.txt',delim_whitespace=True,header=0 ).loc[gene_number])
    df = df.reset_index()
    print(df)
    df = df.loc[df['index'].str.startswith(normalization_technique,na=False)] #This deletes entries not starting with 'fpkm_unstranded'
    df['tumor_status']='no'
    df['cancer_name']=i
    df[gene_number] = df[gene_number].apply(pd.to_numeric)
    print(df)

    protien = df.iloc[1,1]

    df2 = pd.DataFrame(pd.read_csv('data/mRNA_exp_TUMOR_RAW_TXT/'+i+'_tumor.txt',delim_whitespace=True).loc[gene_number])
    df2 = df2.reset_index()
    df2 = df2.loc[df2['index'].str.startswith(normalization_technique,na=False)] #This deletes entries not starting with 'fpkm_unstranded'
    df2['tumor_status']='yes'
    df2['cancer_name']=i
    df2[gene_number] = df2[gene_number].apply(pd.to_numeric)
    print(df2)
   # print(df2[gene_number])
   # df2.loc[gene_number].to_csv(i+'_normal_LTA_csv.csv')

    print(stats.mannwhitneyu(df[gene_number],df2[gene_number]).pvalue)
    df['pvalue']=stats.mannwhitneyu(df[gene_number],df2[gene_number]).pvalue
    df2['pvalue']=stats.mannwhitneyu(df[gene_number],df2[gene_number]).pvalue


    df3 = df3.append([df,df2]) 
    print(df3)
    df3[gene_number] = df3[gene_number].apply(pd.to_numeric)
   # df3['serial_number'] = np.arange(df3.shape[0])
   # df3 = df3.set_index('serial_number')
   # df5.to_csv('CHOL_all_csv.csv')
    
df3[gene_number] = df3[gene_number] + 1

'''
print(df3)
print(list(df3))
print(df3.dtypes)
'''

sns.set()
#sns.set_theme(palette = 'deep')
ax = sns.boxplot(x='cancer_name', y = gene_number , data=df3, hue='tumor_status',whis=200, palette = 'Set1')
ax.set_yscale('log')

pairs = [] #declares array of columns to be compared
for j in cancers:
    pair = [((j, "no"),(j, "yes"))]
    pairs = pairs + pair


annot=Annotator(ax, pairs, data=df3, x='cancer_name', y=gene_number ,  hue='tumor_status') #library to generate significance stars. 
annot.configure(test='Mann-Whitney', verbose=2)
annot.apply_test()
annot.annotate()

plt.title("Pan Cancer Plot " + protien)
plt.ylabel('Log( '+ normalization_technique +' +1)')
ax.set_xticklabels(ax.get_xticklabels(),rotation=40, ha="right")
plt.tight_layout()
plt.show()
sns.plt.savefig('pancancer1.png', dpi=300)
#print(df.loc[14020])

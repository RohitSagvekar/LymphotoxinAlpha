import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 150)

#cancers = ['BRCA','CESC','CHOL','COAD','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','PAAD','PRAD','PCPG','READ','SARC','SKCM','STAD','THCA','THYM','UCEC']

cancers = ['CHOL']

#category = [5, 9, 18]

category = 24

for i in cancers:
   # print(i)
    
    df1 = pd.read_csv('nationwidechildrens.org_clinical_patient_chol.txt',sep='\t',header=None, dtype={'69': 'Int64'})
    df = pd.read_csv('nationwidechildrens.org_clinical_patient_chol.txt',sep='\t',header=None, usecols=[1,category], dtype={'69': 'Int64'})
    df_kirc= pd.read_csv('nationwidechildrens.org_clinical_patient_kirc.txt',sep='\t',header=None,  dtype={'69': 'Int64'})
    df_brca= pd.read_csv('nationwidechildrens.org_clinical_patient_brca.txt',sep='\t',header=None,  dtype={'69': 'Int64'})
    print(df1.T[[0,1]])
    print(df_kirc.T[[0,1]])
    print(df_brca.T[[0,1]])
    

    df2 = pd.read_csv('data/normal_LTA/'+i+'_normal_LTA_csv.csv')
    df3 = pd.read_csv('data/tumor_LTA/'+i+'_normal_LTA_csv.csv')

    
    df2 = df2.loc[df2["Unnamed: 0"].str.startswith('fpkm_unstranded',na=False)] #This deletes entries not starting with 'fpkm_unstranded'
    df2['Unnamed: 0'] = df2['Unnamed: 0'].map(lambda x: x.lstrip('-11A')) #This section deletes parts of string other than patient barcode.
    df2['Unnamed: 0'] = df2['Unnamed: 0'].map(lambda x: x.lstrip('fpkm_unstranded'))
    df2['Unnamed: 0'] = df2['Unnamed: 0'].map(lambda x: str(x)[:-16])

    df3 = df3.loc[df3["Unnamed: 0"].str.startswith('fpkm_unstranded',na=False)] #This deletes entries not starting with 'fpkm_unstranded'
    df3['Unnamed: 0'] = df3['Unnamed: 0'].map(lambda x: x.lstrip('-11A')) #This section deletes parts of string other than patient barcode.
    df3['Unnamed: 0'] = df3['Unnamed: 0'].map(lambda x: x.lstrip('fpkm_unstranded'))
    df3['Unnamed: 0'] = df3['Unnamed: 0'].map(lambda x: str(x)[:-16])

    df = df.add_prefix('col_') #This makes the heading of df a string and removes 1st 3 entries which are gene name, gene number etc
    df = df.drop(df.index[[0,1,2]])   

   # print(df1)
   # print(df)
   # print(df2) 
   # print(df3)


    df2=df2.merge(df, left_on='Unnamed: 0', right_on='col_1' ) #This does the pivot table thingy
    df3=df3.merge(df, left_on='Unnamed: 0', right_on='col_1' ) #This does the pivot table thingy

    
    df2['status'] = 'normal'
    df3['status'] = 'tumor'
    both = pd.concat((df2, df3))

#    print(df)
#    print(df2)
#    print(df3)
    both['27105'] = pd.to_numeric(both['27105'])
   # print(both)

   
''' Wow, finally here concludes the data formatting, now I can peacefully make graphs'''

sns.set()
ax = sns.boxplot(x='status', y='27105', data=both, hue='col_'+str(category),whis=200)
ax.set_yscale('log')
plt.title(i+" "+df1.iat[0,category]+ " boxplot")
plt.ylabel('Log(fpkm+1)')
#ax.set_xticklabels(ax.get_xticklabels(),rotation=40, ha="right")
plt.tight_layout()
plt.show()
sns.plt.savefig('pancancer1.png', dpi=300)
#print(df.loc[14020])

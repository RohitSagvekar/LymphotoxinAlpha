import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import scipy.stats as stats
from scipy.stats import mannwhitneyu

cancers = ['BRCA','CESC','CHOL','COAD','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','PAAD','PRAD','PCPG','READ','SARC','SKCM','STAD','THCA','THYM','UCEC']

#cancers = ['CHOL']

df3 = pd.DataFrame()

for i in cancers:
   # print(i)
    
    df = pd.read_csv('data/normal_LTA/'+i+'_normal_LTA_csv.csv')
    df['tumor_status']='no'
    df = df.drop(df.index[[0,1,2]])   
    df['cancer_name']=i
    df['27105'] = df['27105'].apply(pd.to_numeric)
    print(df) 
    print(df.dtypes)
   
    df2 = pd.read_csv('data/tumor_LTA/'+i+'_normal_LTA_csv.csv')
    df2['tumor_status']='yes'
    df2 = df2.drop(df2.index[[0,1,2]])   
    df2['cancer_name']=i
    df2['27105'] = df2['27105'].apply(pd.to_numeric)
    print(df2)

   # U1, p = mannwhitneyu(df['27105'], df2['27105'], method='exact')
   # print(p)

    df3=df3.append([df,df2]) 
    #df3['27105'] = df3['27105'].apply(pd.to_numeric)
    df3 = df3.loc[df3["Unnamed: 0"].str.startswith('fpkm_unstranded',na=False)]
    df3.to_csv('CHOL_all_csv.csv')
    
df3["27105"] = df3["27105"] + 1
print(df3)
print(list(df3))
print(df3.dtypes)

sns.set()
ax = sns.boxplot(x='cancer_name', y='27105', data=df3, hue='tumor_status',whis=200)
ax.set_yscale('log')
plt.title("Pan Cancer Plot")
plt.ylabel('Log(fpkm+1)')
ax.set_xticklabels(ax.get_xticklabels(),rotation=40, ha="right")
plt.tight_layout()
plt.show()
sns.plt.savefig('pancancer1.png', dpi=300)
#print(df.loc[14020])



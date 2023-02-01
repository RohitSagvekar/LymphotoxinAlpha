import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

#cancers = ['BRCA','CESC','CHOL','COAD','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','PAAD','PRAD','PCPG','READ','SARC','SKCM','STAD','THCA','THYM','UCEC']

cancers = ['CHOL']


for i in cancers:
   # print(i)
    
    df = pd.read_csv('nationwidechildrens.org_clinical_patient_chol.txt',sep='\t',header=None,usecols=[0,1,13,69,76], dtype={'69': 'Int64'})
    #df = pd.read_csv('nationwidechildrens.org_clinical_patient_kirc.txt',sep='\t',header=None, dtype={'69': 'Int64'})
    df = df.drop(df.index[[0,1,2]])   
    df = df.add_prefix('col_')
    
    print(df)
    
    #A little insertion from hispanicBoxplot here now

    df2 = pd.read_csv('data/normal_LTA/'+i+'_normal_LTA_csv.csv')
    df3 = pd.read_csv('data/tumor_LTA/'+i+'_normal_LTA_csv.csv')

    
    df2 = df2.loc[df2["Unnamed: 0"].str.startswith('fpkm_unstranded',na=False)] #This deletes entries not starting with 'fpkm_unstranded'
    df2['Unnamed: 0'] = df2['Unnamed: 0'].map(lambda x: x.lstrip('-11A')) #This section deletes parts of string other than patient barcode.
    df2['Unnamed: 0'] = df2['Unnamed: 0'].map(lambda x: x.lstrip('fpkm_unstranded'))
    df2['Unnamed: 0'] = df2['Unnamed: 0'].map(lambda x: str(x)[:-16])
    print(df2)

    df3 = df3.loc[df3["Unnamed: 0"].str.startswith('fpkm_unstranded',na=False)] #This deletes entries not starting with 'fpkm_unstranded'
    df3['Unnamed: 0'] = df3['Unnamed: 0'].map(lambda x: x.lstrip('-11A')) #This section deletes parts of string other than patient barcode.
    df3['Unnamed: 0'] = df3['Unnamed: 0'].map(lambda x: x.lstrip('fpkm_unstranded'))
    df3['Unnamed: 0'] = df3['Unnamed: 0'].map(lambda x: str(x)[:-16])
    print(df3)

    #df3['27105'] = pd.to_numeric(df3['27105'], errors='coerce')
    #df4 = df3[df3['27105'].between(0,df3.median())]
   # print(df3['27105'].between(0,df3.median()))
    print(df3.median())

   # df = df.add_prefix('col_') #This makes the heading of df a string and removes 1st 3 entries which are gene name, gene number etc
   # df = df.drop(df.index[[0,1,2]])   

    #df2 = 0

    df2=df2.merge(df, left_on='Unnamed: 0', right_on='col_1' ) #This does the pivot table thingy
    df3=df3.merge(df, left_on='Unnamed: 0', right_on='col_1' ) #This does the pivot table thingy

#    print(df2)
#    print(df3)
#    print(df)

    

    df = df.sort_values(by=['col_69'], ignore_index=True)
    df['col_69'] = pd.to_numeric(df['col_69'], errors='coerce')
    df['col_76'] = pd.to_numeric(df['col_76'], errors='coerce')
    df['col_69'] = df['col_69'].fillna(0) + df['col_76'].fillna(0)
    df = df.drop(columns=['col_76'])   
    region_dictionary = {'Dead': 1, 'Alive' : 0}
    df['col_13'] = df['col_13'].apply(lambda x: region_dictionary[x])
    #df.to_csv('CHOL_clinical_csv.csv')
    
    # Prepare unique durations in ascending order
    durations = df.sort_values('col_69')['col_69'].unique()

    # Initialise the table
    columns = ['duration', 'n_at_risk', 'n_events','survival_probability']
    km = pd.DataFrame(columns=columns, dtype=np.number)
    km = km.append(pd.DataFrame([[0, df.shape[0], 0, 1]],columns=columns))

    # Calculate survival probability for each duration
    for i, t in enumerate(durations):
        n = np.sum(df['col_69']>=t)
        d = np.sum((df['col_69']==t) & (df['col_13']==1))
        s = (1 - d / n) * km.loc[i, 'survival_probability']
        km = km.append(pd.DataFrame([[t, n, d, s]],index=[i+1],columns=columns))

    print(km)

'''
#yaha thoda bakchodi karenge 

    df1 = pd.read_csv('nationwidechildrens.org_clinical_patient_chol.txt',sep='\t',header=None, usecols=[1,9], dtype={'69': 'Int64'})
    print(df1.T[[0,1]])

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

    df1 = df.add_prefix('col_') #This makes the heading of df a string and removes 1st 3 entries which are gene name, gene number etc
    df1 = df.drop(df.index[[0,1,2]])   

    print(df1)
    print(df)
    print(df2) 
    print(df3)


    df2=df2.merge(df1, left_on='Unnamed: 0', right_on='col_1' ) #This does the pivot table thingy
    df3=df3.merge(df1, left_on='Unnamed: 0', right_on='col_1' ) #This does the pivot table thingy

    
    df2['status'] = 'normal'
    df3['status'] = 'tumor'
    both = pd.concat((df2, df3))

#    print(df)
#    print(df2)
#    print(df3)
    both['27105'] = pd.to_numeric(both['27105'])
    print(both)




'''

print(km)
print(list(df))
print(df.dtypes)

sns.set()
ax = sns.lineplot(x='duration', y='survival_probability', data=km,drawstyle='steps-post', palette = 'Set1')
plt.title("CHOL Kaplan-Meier Curve")
plt.ylim(0,1.1)
plt.show()
#plt.savefig('Kaplan-Meier_CHOL', dpi=300)

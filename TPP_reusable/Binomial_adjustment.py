import pandas as pd
from scipy.stats import binom
import matplotlib.pyplot as plt
import re
import numpy as np

def calulate_decoy_FLR(input,decoy,targets):
    df = pd.read_csv(input,dtype={'PTM_positions': str})
    STY_count=0
    for target in list(targets):
        T_count=df['Peptide'].str.count(target).sum()
        STY_count+=T_count
    A_count=df['Peptide'].str.count(decoy).sum()
    STY_A_ratio=STY_count/A_count
    pA_count = 0
    df.fillna('-')
    for i in range(len(df)):
        for a in (str(df.loc[i,'PTM positions']).split(";")):
            if a!="-" and a!="nan":
                peptide=re.sub('[()+0-9. "]', '',df.loc[i,'Peptide'])
                if peptide[int(float(a))-1]==decoy:
                    pA_count+=1
                    df.loc[i,'DecoyP'] = 1
                else:
                    df.loc[i,'DecoyP'] = 0
        df.loc[i, 'p'+decoy+'_count'] = pA_count
        #decoy pX_FLR = STY:X ratio * pX_count * 2 / Count
        df.loc[i,decoy+'_FLR']=(STY_A_ratio*df.loc[i,'p'+decoy+'_count']*2)/(i+1)

    df[decoy+'_q_value'] = df[decoy+'_FLR']
    df[decoy+'_q_value'] = df.iloc[::-1][decoy+'_FLR'].cummin()
    df.to_csv(input,index=False)

def model_FLR_binomial(file,FLR_output):
    df=pd.read_csv(file)
    df = df.sort_values(['Peptide_pos', 'Binomial_final_score'], ascending=[True, True])
    USI_list = []
    PSM_threshold = []
    sig1=[]
    sig2=[]
    sig3=[]
    sig4=[]
    sig5=[]
    sig6=[]
    for i in df['USI'].groupby(df['Peptide_pos']).apply(';'.join):
        USI_list.append(i)
    for i in df.groupby(['Peptide_pos'])['0.05FLR_threshold'].sum():
        PSM_threshold.append(i)
    for i in df.groupby(['Peptide_pos'])['0.01<P<=0.05'].sum():
        sig1.append(i)
    for i in df.groupby(['Peptide_pos'])['0.05<P<=0.19'].sum():
        sig2.append(i)
    for i in df.groupby(['Peptide_pos'])['0.19<P<=0.81'].sum():
        sig3.append(i)
    for i in df.groupby(['Peptide_pos'])['0.81<P<=0.95'].sum():
        sig4.append(i)
    for i in df.groupby(['Peptide_pos'])['0.95<P<0.99'].sum():
        sig5.append(i)
    for i in df.groupby(['Peptide_pos'])['P>=0.99'].sum():
        sig6.append(i)
    df = df.drop_duplicates(subset=('Peptide_pos'), keep='last', inplace=False)
    df['All_USI'] = USI_list
    df['0.05FLR_threshold_count'] = PSM_threshold
    df['0.01<P<=0.05_count']=sig1
    df['0.05<P<=0.19_count']=sig2
    df['0.19<P<=0.81_count']=sig3
    df['0.81<P<=0.95_count']=sig4
    df['0.95<P<0.99_count']=sig5
    df['P>=0.99_count']=sig6
    df = df.reset_index(drop=True)
    df = df.sort_values(by=(['Binomial_final_score']), ascending=[False])
    df = df.reset_index(drop=True)
    df['Count_binomial'] = (df.index) + 1
    df['final_temp'] = 1 - df['Binomial_final_score']
    df['Binomial_final_prob_FLR'] = df['final_temp'].cumsum() / df['Count_binomial']
    df['Binomial_final_prob_q_value'] = df['Binomial_final_prob_FLR']
    df['Binomial_final_prob_q_value'] = df.iloc[::-1]['Binomial_final_prob_FLR'].cummin()
    df = df.drop(['final_temp'], axis=1)
    df = df.reset_index(drop=True)
    df.to_csv(FLR_output,index=False)


def Binomial(file,decoy, targets):
    d="/".join(file.split("/")[:-1])
    output = d+"/binomial.csv"
    FLR_output = d+"/binomial_collapsed_FLR.csv"
    df = pd.read_csv(file)
    df['0.05FLR_threshold']=np.where(df['p'+decoy+'_FLR']<=0.05,1,0)
    df['Protein-pos'] = df['Protein']+"-"+df['Protein position'].astype(str)
    df['Peptide_pos']  = df['Peptide']+"-"+df['PTM positions'].astype(str)
    df['Peptide_start_Protein'] = df['Protein position'] - df['PTM positions'] +1
    df['Peptide_end_Protein'] = df['Peptide_start_Protein']+(df['Peptide'].str.len()) -1
    df = df.sort_values(by=(['Protein','Protein position']), ascending=[False,False])
    df['Protein_loc'] = df.groupby('Protein-pos')['Protein-pos'].transform('count')
    counter=1
    for i in range(len(df)):
        protein=df.loc[i,'Protein-pos'].rsplit("-",1)[0]
        pos=int(df.loc[i,'Protein-pos'].rsplit("-",1)[-1])
        df_temp=df.loc[(df['Protein']==protein) & (df['Peptide_start_Protein']<=pos) & (df['Peptide_end_Protein']>=pos)]
        df.loc[i,'Protein_occ']=len(df_temp)
        counter+=1
    prob_chance=df['pA_count'].max()/len(df)
    df.loc[:, 'pA_prob_chance_binomial'] = prob_chance
    df['Binomial_p'] = binom.pmf(df['Protein_loc'],df['Protein_occ'],prob_chance)
    df['Binomial_i'] = df['PTM Score'] * (1-df['Binomial_p'])
    df['Binomial_final_score'] = df['Binomial_i'] * df['Score']
    df['0.01<P<=0.05']=np.where(((df['Binomial_final_score']>0.01)&(df['Binomial_final_score']<=0.05)),1,0)
    df['0.05<P<=0.19']=np.where(((df['Binomial_final_score']>0.05)&(df['Binomial_final_score']<=0.19)),1,0)
    df['0.19<P<=0.81']=np.where(((df['Binomial_final_score']>0.19)&(df['Binomial_final_score']<=0.81)),1,0)
    df['0.81<P<=0.95']=np.where(((df['Binomial_final_score']>0.81)&(df['Binomial_final_score']<=0.95)),1,0)
    df['0.95<P<0.99']=np.where(((df['Binomial_final_score']>0.95)&(df['Binomial_final_score']<0.99)),1,0)
    df['P>=0.99']=np.where(df['Binomial_final_score']>=0.99,1,0)
    df.to_csv(output,index=False)

    model_FLR_binomial(output,FLR_output)
    calulate_decoy_FLR(FLR_output,decoy,targets)
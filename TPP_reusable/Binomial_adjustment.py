import pandas as pd
from scipy.stats import binom
import numpy as np
import time

def calulate_decoy_FLR(input,decoy,targets):
    print("Running: Bimonial calulate_decoy_FLR")
    start_time = time.time()
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
        if decoy+"[" in df.loc[i,'Peptide_mod']:
            pA_count+=1
            df.loc[i,'DecoyP'] = 1
            if df.loc[i,'Peptide_mod'].startswith(decoy+"[Acetyl]") and df.loc[i,'Protein position']==2:
                df.loc[i,'DecoyP'] = 0
                pA_count-=1
        else:
            df.loc[i,'DecoyP'] = 0

        df.loc[i, 'p'+decoy+'_count'] = pA_count
        #decoy pX_FLR = STY:X ratio * pX_count * 2 / Count
        #df.loc[i,decoy+'_FLR']=(STY_A_ratio*df.loc[i,'p'+decoy+'_count']*2)/(i+1)
        #df.loc[i,'p'+decoy+'_FLR']=((STY_A_ratio*df.loc[i,'p'+decoy+'_count'])+df.loc[i,'p'+decoy+'_count'])/(i+1)
        df.loc[i,'p'+decoy+'_FLR_BA']=(STY_A_ratio*df.loc[i,'p'+decoy+'_count'])/(i+1-df.loc[i,'p'+decoy+'_count'])

    df['p'+decoy+'_q_value_BA'] = df['p'+decoy+'_FLR_BA']
    df['p'+decoy+'_q_value_BA'] = df.iloc[::-1]['p'+decoy+'_FLR_BA'].cummin()
    df.to_csv(input,index=False)
    print("Complete --- %s seconds ---" % (time.time() - start_time))

def model_FLR_binomial(file,FLR_output):
    print("Running: model_FLR_binomial")
    start_time = time.time()
    df=pd.read_csv(file)
    df = df.sort_values(['Binomial_final_score','Peptide','Peptide_pos'], ascending=[True, True,True])
    USI_list = []
    PSM_threshold = []
    sig1=[]
    sig2=[]
    sig3=[]
    sig4=[]
    sig5=[]
    sig6=[]
    df['Peptide_mod_pos']=df['Peptide_mod']+"-"+df['PTM positions'].astype(str)
    for i in df['USI'].groupby(df['Peptide_mod_pos']).apply(';'.join):
        USI_list.append(i)
    for i in df.groupby(['Peptide_mod_pos'])['0.05FLR_threshold'].sum():
        PSM_threshold.append(i)
    for i in df.groupby(['Peptide_mod_pos'])['0.01<P<=0.05'].sum():
        sig1.append(i)
    for i in df.groupby(['Peptide_mod_pos'])['0.05<P<=0.19'].sum():
        sig2.append(i)
    for i in df.groupby(['Peptide_mod_pos'])['0.19<P<=0.81'].sum():
        sig3.append(i)
    for i in df.groupby(['Peptide_mod_pos'])['0.81<P<=0.95'].sum():
        sig4.append(i)
    for i in df.groupby(['Peptide_mod_pos'])['0.95<P<0.99'].sum():
        sig5.append(i)
    for i in df.groupby(['Peptide_mod_pos'])['P>=0.99'].sum():
        sig6.append(i)

    #original collapse method
    #df = df.drop_duplicates(subset=('Peptide_pos'), keep='last', inplace=False)
    #updated collapse - peptidoform-site based
    df=df.drop_duplicates(subset=(['Peptide_mod_pos']),keep='last',inplace=False)

    df['All_USI'] = USI_list
    df['0.05FLR_threshold_count'] = PSM_threshold
    df['0.01<P<=0.05_count']=sig1
    df['0.05<P<=0.19_count']=sig2
    df['0.19<P<=0.81_count']=sig3
    df['0.81<P<=0.95_count']=sig4
    df['0.95<P<0.99_count']=sig5
    df['P>=0.99_count']=sig6
    df = df.reset_index(drop=True)
    df['Peptide']=df['Peptide'].astype(str)
    df = df.sort_values(by=(['Binomial_final_score','Peptide','Peptide_pos']), ascending=[False,True,True])
    df = df.reset_index(drop=True)
    df['Count_binomial'] = (df.index) + 1
    df['final_temp'] = 1 - df['Binomial_final_score']
    df['Binomial_final_prob_FLR'] = df['final_temp'].cumsum() / df['Count_binomial']
    df['Binomial_final_prob_q_value'] = df['Binomial_final_prob_FLR']
    df['Binomial_final_prob_q_value'] = df.iloc[::-1]['Binomial_final_prob_FLR'].cummin()
    df = df.drop(['final_temp'], axis=1)
    df = df.reset_index(drop=True)
    df.to_csv(FLR_output,index=False)
    print("Complete --- %s seconds ---" % (time.time() - start_time))

def prot_occ(df,protein,pos):
    df_temp=df.loc[(df['Protein']==protein) & (df['Peptide_start_Protein']<=pos) & (df['Peptide_end_Protein']>=pos)]
    df_temp=df_temp.drop_duplicates(subset=(['USI']),keep='last',inplace=False)
    return(len(df_temp))

def Binomial(file,decoy, targets):
    print("Running: Binomial")
    start_time = time.time()
    d="/".join(file.split("/")[:-1])
    output = d+"/binomial.csv"
    FLR_output = d+"/binomial_peptidoform_collapsed_FLR.csv"
    df = pd.read_csv(file)
    df['0.05FLR_threshold']=np.where(df['p'+decoy+'_FLR']<=0.05,1,0)
    df['Protein-pos'] = df['Protein']+"-"+df['Protein position'].astype(str)
    df['Peptide_pos']  = df['Peptide']+"-"+df['PTM positions'].astype(str)
    df['Peptide_start_Protein'] = df['Protein position'] - df['PTM positions'] +1
    df['Peptide_end_Protein'] = df['Peptide_start_Protein']+(df['Peptide'].str.len()) -1
    df = df.sort_values(by=(['Protein','Protein position']), ascending=[False,True])
    df['Protein_loc'] = df.groupby('Protein-pos')['Protein-pos'].transform('count')
    counter=1

    #dict of occ for protein-pos
    occ_dict = {}
    for p in df['Protein-pos'].unique():
        protein=p.rsplit("-",1)[0]
        pos=int(p.rsplit("-",1)[-1])
        df_temp=df.loc[(df['Protein']==protein) & (df['Peptide_start_Protein']<=pos) & (df['Peptide_end_Protein']>=pos)]
        df_temp=df_temp.drop_duplicates(subset=(['USI']),keep='last',inplace=False)
        occ=len(df_temp)
        occ_dict[p]=occ

    # for i in range(len(df)):
    #     protein=df.loc[i,'Protein-pos'].rsplit("-",1)[0]
    #     pos=int(df.loc[i,'Protein-pos'].rsplit("-",1)[-1])
    #     df_temp=df.loc[(df['Protein']==protein) & (df['Peptide_start_Protein']<=pos) & (df['Peptide_end_Protein']>=pos)]
    #     df_temp=df_temp.drop_duplicates(subset=(['USI']),keep='last',inplace=False)
    #     df.loc[i,'Protein_occ']=len(df_temp)
    #df['Protein_occ']=df.apply(lambda row: prot_occ(df,row['Protein'],row['Protein position']),axis=1)
    df['Protein_occ']=df['Protein-pos'].map(occ_dict)
    counter+=1
    #prob_chance=df['pA_count'].max()/len(df) #/number of PSMs not sites
    unique_USIs=df['USI'].unique()
    prob_chance=sum(decoy+"[" in s for s in unique_USIs)/len(unique_USIs)
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
    df['Peptide']=df['Peptide'].astype(str)
    df = df.sort_values(by=(['Binomial_final_score','Peptide','Peptide_pos']), ascending=[False,True,True])
    df.to_csv(output,index=False)
    print("Complete --- %s seconds ---" % (time.time() - start_time))

    model_FLR_binomial(output,FLR_output)
    calulate_decoy_FLR(FLR_output,decoy,targets)
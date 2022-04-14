import pandas as pd
from scipy.stats import binom
import matplotlib.pyplot as plt
import re
import numpy as np

def calulate_decoy_FLR(input,decoy):
    if decoy=="pAla":
        x="A"
    if decoy=="pLeu":
        x="L"
    if decoy=="pGly":
        x="G"
    if decoy=="pAsp":
        x="D"
    if decoy=="pGlu":
        x="E"
    if decoy=="pPro":
        x="P"
    df = pd.read_csv(input,dtype={'PTM positions': str})
    S_count=df['Peptide'].str.count('S').sum()
    T_count=df['Peptide'].str.count('T').sum()
    Y_count=df['Peptide'].str.count('Y').sum()
    A_count=df['Peptide'].str.count(x).sum()
    STY_count=S_count+T_count+Y_count
    STY_A_ratio=STY_count/A_count

    pA_count = 0
    df.fillna('-')
    for i in range(len(df)):
        for a in (str(df.loc[i,'PTM positions']).split(";")):
            if a!="-" and a!="nan":
                peptide=re.sub('[()+0-9. "]', '',df.loc[i,'Peptide'])
                if peptide[int(float(a))-1]==x:
                    pA_count+=1
                    df.loc[i,'DecoyP'] = 1
                else:
                    df.loc[i,'DecoyP'] = 0
        df.loc[i, 'p'+x+'_count'] = pA_count
        #decoy pX_FLR = STY:X ratio * pX_count * 2 / Count
        df.loc[i,decoy+'_FLR']=(STY_A_ratio*df.loc[i,'p'+x+'_count']*2)/(i+1)

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

def plot_FLR_comparisons(file_list,file_name_list,working):
    import seaborn as sns;sns.set_palette("colorblind")
    legend=[]
    for j,k in zip(file_list,file_name_list):
        print(j,k)
        if "binomial" in j:
            legend.append(k+"_Binomial pAla Decoy FLR")
        elif "pL" in j or "pLeu" in j:
            legend.append(k+"_pLeu Decoy FLR")
        elif "pGlu" not in j and ("pG" in j or "pGly" in j):
            legend.append(k+"_pGly Decoy FLR")
        elif "pD" in j or "pAsp" in j:
            legend.append(k+"_pAsp Decoy FLR")
        elif "pE" in j or "pGlu" in j:
            legend.append(k+"_pGlu Decoy FLR")
        elif "pP" in j or "pPro" in j:
            legend.append(k+"_pPro Decoy FLR")
        else:
            legend.append(k+"_pAla Decoy FLR")
    print(legend)
    output = working + "/Comparisons/"
    #Combined prob FLR plot
    x_list=['Count']
    for x in x_list:
        counter = 0
        max = 0
        for a in file_list:
            if "pGly" in a:
                p="pGly_q_value"
            elif "pLeu" in a:
                p="pLeu_q_value"
            elif "pAsp" in a:
                p="pAsp_q_value"
            elif "pGlu" in a:
                p="pGlu_q_value"
            elif "pPro" in a:
                p="pPro_q_value"
            else:
                p="pAla_q_value"
            df=pd.read_csv(a)
            if "binomial_collapsed_FLR.csv" in a:
                if df['Count_binomial'].max() > max:
                    max = df['Count_binomial'].max()
            else:
                if df[x].max() > max:
                    max = df[x].max()
            if x == "PTM_final_prob":
                lim = (1, 0)
            else:
                print(max)
                lim=(0,max)
            c = "C" + str(counter + 1)
            if "binomial_collapsed_FLR.csv" not in a:
                if counter == 0:
                    df = df.sort_values(by=(['Count']), ascending=[False])
                    ax = df.plot.line(x=x, y=p, linestyle="-", color=c, xlim=lim,figsize=(14, 7))
                else:
                    df = df.sort_values(by=(['Count']), ascending=[False])
                    df.plot.line(x=x, y=p, linestyle="-", ax=ax, color=c, xlim=lim)
                counter += 1
            else:
                if counter == 0:
                    df = df.sort_values(by=(['Count_binomial']), ascending=[False])
                    df.plot.line(x="Count_binomial", y="pAla_q_value", linestyle="-", color=c, xlim=lim, ax=ax)
                else:
                    df = df.sort_values(by=(['Count_binomial']), ascending=[False])
                    df.plot.line(x="Count_binomial", y="pAla_q_value", linestyle="-", color=c, xlim=lim, ax=ax)
                counter += 1
        ax.legend(legend)
        ax.set_ylabel("FLR")
        ax.set_xlabel("Count of Sites")
        output_file=output+"Binomial.jpg"
        print(output_file)
        plt.savefig(output_file, dpi=300)

def Binomial(file,decoy):
    #file = "D:/Dropbox/PTMExchange/PXD008355/RAPA_tryptic/TPP_binomial_test/FDR_0.01_PTM_score_0/ALL_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
    collapsed_output = file.replace(".csv", "_FLR_collapse_pep_pos_final_score.csv")
    d="/".join(file.split("/")[:-1])
    print(d)
    output = d+"/binomial.csv"
    FLR_output = d+"/binomial_collapsed_FLR.csv"
    output_decoy = d+"/binomial_no_decoy.csv"
    FLR_output_decoy = d+"/binomial_no_decoy_collapsed_FLR.csv"
    print(output)
    df = pd.read_csv(file)
    df['0.05FLR_threshold']=np.where(df['pAla_FLR']<=0.05,1,0)
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
        #print(str(counter)+"/"+str(len(df)))
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

    #remove decoys before collapse and scoring
    df = pd.read_csv(file)
    if df.loc[i,'DecoyP'] == 1:
        #Replace pA with highest scoring STY
        PTM_info=df.loc[i,'PTM_info'].split(";")
        USI_temp=df.loc[i,'USI']

        mods=[]
        mods_prob=[]
        mods_pos=[]

        max_score=-1
        max_score_pos=-1
        for z in PTM_info:
            pos=int(z.split(':')[0])
            mod=z.split(':')[1]
            score=float(z.split(':')[2])
            if mod=="Phosphorylation" and score>=max_score and max_score<=df.loc[i,'PTM Score'] and str(pos)!=df.loc[i,'PTM positions'] and pos not in mods_pos and df.loc[i,'Peptide'][pos-1]!="A":
                max_score=score
                max_score_pos=pos
        mods.append("Phosphorylation")
        mods_prob.append(float(max_score))
        mods_pos.append(int(max_score_pos))
        phos_pos=max_score_pos

        #reorder by position
        zipped=sorted(zip(mods_pos,mods,mods_prob))
        mods_pos,mods,mods_prob=zip(*zipped)

        mod_list=""
        for a in list(mods):
            mod_list+=a+";"

        pos_list=""
        for a in list(mods_pos):
            pos_list+=str(a)+";"

        prob_list=""
        final_probs_list=""
        for a in list(mods_prob):
            prob_list+=str(a)+";"
            final_probs_list+=str(a*df.loc[i,'Score'])+";"

        df.loc[i,'PTM']=mod_list[:-1]
        df.loc[i,'PTM Score']=prob_list[:-1]
        df.loc[i,'PTM_final_prob']=final_probs_list[:-1]
        df.loc[i,'PTM positions']=phos_pos

        mod_pep=USI_temp.split(':')[-1].split("/")[0]
        seq=df.loc[i,'Peptide']
        for x,y in zip(mods_pos[::-1],mods[::-1]):
            if y=="Acetyl":
                seq=seq[:x]+"["+y+"]-"+seq[x:]
            else:
                if "Acetyl" in seq:
                    seq = seq[:x + 9] + "[" + y + "]" + seq[x + 9:]
                else:
                    seq=seq[:x]+"["+y+"]"+seq[x:]
        df.loc[i,'USI']=df.loc[i,'USI'].replace(mod_pep,seq)
        df.loc[i,'Peptide_mod']=seq

    df['0.05FLR_threshold'] = np.where(df['pAla_FLR'] <= 0.05, 1, 0)
    df['Protein-pos'] = df['Protein']+"-"+df['Protein position'].astype(str)
    df['Peptide_pos']  = df['Peptide']+"-"+df['PTM positions'].astype(str)
    df['Peptide_start_Protein'] = df['Protein position'] - df['PTM positions'] +1
    df['Peptide_end_Protein'] = df['Peptide_start_Protein']+(df['Peptide'].str.len()) -1
    df = df.sort_values(by=(['Protein','Protein position']), ascending=[False,False])
    df['Protein_loc'] = df.groupby('Protein-pos')['Protein-pos'].transform('count')
    counter=1
    for i in range(len(df)):
        protein = df.loc[i, 'Protein-pos'].rsplit("-", 1)[0]
        pos = int(df.loc[i, 'Protein-pos'].rsplit("-", 1)[-1])
        df_temp=df.loc[(df['Protein']==protein) & (df['Peptide_start_Protein']<=pos) & (df['Peptide_end_Protein']>=pos)]
        df.loc[i,'Protein_occ']=len(df_temp)
        #print(str(counter)+"/"+str(len(df)))
        counter+=1
    prob_chance=df['pA_count'].max()/len(df)
    df.loc[:,'pA_prob_chance_binomial']=prob_chance
    df['Binomial_p'] = binom.pmf(df['Protein_loc'],df['Protein_occ'],prob_chance)
    df['Binomial_i'] = df['PTM Score'].astype(float) * (1-df['Binomial_p'])
    df['Binomial_final_score'] = df['Binomial_i'] * df['Score']
    df['0.01<P<=0.05']=np.where(((df['Binomial_final_score']>0.01)&(df['Binomial_final_score']<=0.05)),1,0)
    df['0.05<P<=0.19']=np.where(((df['Binomial_final_score']>0.05)&(df['Binomial_final_score']<=0.19)),1,0)
    df['0.19<P<=0.81']=np.where(((df['Binomial_final_score']>0.19)&(df['Binomial_final_score']<=0.81)),1,0)
    df['0.81<P<=0.95']=np.where(((df['Binomial_final_score']>0.81)&(df['Binomial_final_score']<0.95)),1,0)
    df['0.95<P<0.99']=np.where(((df['Binomial_final_score']>0.95)&(df['Binomial_final_score']<=0.99)),1,0)
    df['P>=0.99']=np.where(df['Binomial_final_score']>=0.99,1,0)
    df.to_csv(output_decoy,index=False)

    df=pd.read_csv(file)
    df['Peptide_pos']  = df['Peptide']+"-"+df['PTM positions'].astype(str)
    df=df.sort_values(['Peptide_pos','PTM_final_prob'],ascending=[True,True])
    USI_list = []
    for i in df['USI'].groupby(df['Peptide_pos']).apply(';'.join):
        USI_list.append(i)

    df=df.drop_duplicates(subset=('Peptide_pos'),keep='last', inplace=False)
    df['All_USI'] = USI_list
    df=df.sort_values(['PTM_final_prob'],ascending=[False])
    df = df.reset_index(drop=True)
    df['Count'] = (df.index) + 1
    df['final_temp'] = 1 - df['PTM_final_prob']
    df['PTM_final_prob_FLR'] = df['final_temp'].cumsum() / df['Count']
    df['PTM_final_prob_q_value'] = df['PTM_final_prob_FLR']
    df['PTM_final_prob_q_value'] = df.iloc[::-1]['PTM_final_prob_FLR'].cummin()
    df = df.drop(['final_temp'], axis=1)
    df = df.reset_index(drop=True)
    df.to_csv(collapsed_output,index=False)

    model_FLR_binomial(output,FLR_output)
    model_FLR_binomial(output_decoy,FLR_output_decoy)
    calulate_decoy_FLR(FLR_output,decoy)
    calulate_decoy_FLR(FLR_output_decoy,decoy)
    calulate_decoy_FLR(collapsed_output,decoy)




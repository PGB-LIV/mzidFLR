#Calculates decoy amino acid FLR 
#Takes input from "FLR files" generated from Post_analysis and calculates Decoy amino acid FLR. 
#output=[input]_p[Decoy amino acid].csv


import pandas as pd
import re

#calculate decoy FLR from decoy input 
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
                #print(a, ":" , peptide)
                if peptide[int(float(a))-1]==decoy:
                    pA_count+=1
                    df.loc[i,'DecoyP'] = 1
                    if mod.lower()=="acetyl" and a==0:
                        df.loc[i,'DecoyP'] = 0
                else:
                    df.loc[i,'DecoyP'] = 0
        df.loc[i, 'p'+decoy+'_count'] = pA_count
        #decoy pX_FLR = STY:X ratio * pX_count * 2 / Count
        df.loc[i,'p'+decoy+'_FLR']=(STY_A_ratio*df.loc[i,'p'+decoy+'_count']*2)/(i+1)

    df['p'+decoy+'_q_value'] = df['p'+decoy+'_FLR']
    df['p'+decoy+'_q_value'] = df.iloc[::-1]['p'+decoy+'_FLR'].cummin()
    output=input.replace(".csv","_p"+decoy+".csv")
    df.to_csv(output,index=False)
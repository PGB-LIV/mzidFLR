#Calculates decoy amino acid FLR 
#Takes input from "FLR files" generated from Post_analysis and calculates Decoy amino acid FLR. 
#output=[input]_p[Decoy amino acid].csv


import pandas as pd
import re

#calculate decoy FLR from decoy input 
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
    df = pd.read_csv(input,dtype={'PTM_positions': str})
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
                #print(a, ":" , peptide)
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
    output=input.replace(".csv","_"+decoy+".csv")
    df.to_csv(output,index=False)
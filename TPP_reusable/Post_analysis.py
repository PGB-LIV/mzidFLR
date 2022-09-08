#Expands PSMs to site-based format, where one row is one site on a peptide. Calculates FLR calculations.


import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

def explode(df, lst_cols, fill_value='', preserve_index=False):
    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
        and len(lst_cols) > 0
        and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
                col:np.repeat(df[col].values, lens)
                for col in idx_cols},
                index=idx)
             .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
                            for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (res.append(df.loc[lens==0, idx_cols], sort=False)
                  .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:
        res = res.reset_index(drop=True)
    return res

#Explode rows for each PSM postion - for each mod in peptide, a seperate row will show PSM score and corresponding PTM position/score for each site on the peptide
def site_based(input, FDR_cutoff):
    df = pd.read_csv(input, dtype={'Protein position': str, 'PTM Score':str, 'PTM positions':str})
    df = df.loc[df['FDR'] <= FDR_cutoff]
    df = df.sort_values(['Peptide', 'Score','PTM Score'], ascending=[True, True, True])
    df['All_Proteins']=df['Protein']
    df['All_USI']=df['USI']

    df.dropna(subset=['PTM'], inplace=True)
    df=df.fillna('')
    df['All_PTMs'] = df['PTM']
    df['All_PTM_scores'] = df['PTM Score']
    df['All_PTM_positions'] = df['PTM positions']
    df['All_PTM_protein_positions'] = df['Protein position']
    df['PTM'] = df['PTM'].str.split(';')
    df['PTM Score'] = df['PTM Score'].str.split(';')
    df['PTM positions'] = df['PTM positions'].str.split(';')
    df3 = df.loc[~df['Protein position'].str.contains(':')]
    df3['Protein'] = df3['All_Proteins']
    df2 = df.loc[df['Protein position'].str.contains(':')]
    df2['Protein position'] = df2['Protein position'].str.split(':').str[0]
    df2['Protein'] = df2['All_Proteins'].str.split(':').str[0]
    df = pd.concat([df3,df2])
    df['Protein position'] = df['Protein position'].str.split(';')
    df = explode(df, ['PTM', 'PTM Score', 'PTM positions', 'Protein position'], fill_value='')
    df = df[df.PTM == "Phospho"]
    df = df[~df.Protein.str.contains("DECOY")]
    df = df[~df.Protein.str.contains("CONTAM")]
    df = df.reset_index(drop=True)
    output = input.replace("FDR_output.csv", "Site-based.csv")
    df.to_csv(output,index=False)

#model FLR from combined probability: 1-(PSM prob*PTM prob)/Count
def model_FLR(file):
    df=pd.read_csv(file, dtype={'Score': float,'PTM Score':float})
    df = df[df['PTM positions'].notna()]
    df = df[df.PTM == "Phospho"]
    df = df[~df.Protein.str.contains("DECOY",na=False)]
    df = df[~df.Protein.str.contains("CONTAM",na=False)]
    df = df.reset_index(drop=True)
    df['PTM_final_prob'] = df['Score'] * df['PTM Score']
    df = df.sort_values(by=(['PTM_final_prob']), ascending=[False])
    df = df.reset_index(drop=True)
    df['Count'] = (df.index) + 1
    df['final_temp'] = 1 - df['PTM_final_prob']
    df['PTM_final_prob_FLR'] = df['final_temp'].cumsum() / df['Count']
    df['PTM_final_prob_q_value'] = df['PTM_final_prob_FLR']
    df['PTM_final_prob_q_value'] = df.iloc[::-1]['PTM_final_prob_FLR'].cummin()
    df = df.drop(['final_temp'], axis=1)
    df = df.reset_index(drop=True)
    output=file.replace(".csv","_FLR.csv")
    df.to_csv(output,index=False)

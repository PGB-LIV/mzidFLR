#Gives collapsed outputs of unique hits
#	a) collapse by mass shift - collapsing PSMs by peptide+mass shift, keeping the best scoring: "Peptide_mass_confident_PTM_unique.csv"
#	b) collapse by peptide  - collapsing PSMs by peptide, keeping best scoring: "Peptide_confident_PTM_unique.csv"
#	c) collapse by protein position - collapsing modification by position on protein, keeping best scoring: "Site_confident_PTM_unique.csv"
#Also gives none-collapsed in same format: "All_confident_no_collapse.csv

import pandas as pd
import re
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

#Collapse for best scoring for each peptide, protein site, mass shift on protein or no collapse
def unique(output,FDR_cutoff, peptide_based,site_based,mass_based,non_collapse):
    #filter for FDR here
    df = pd.read_csv(output)
    df = df.loc[df['FDR']<=FDR_cutoff]
    #Collapse by mass shift
    df=df.sort_values(['Peptide','Mass shift','Score','PTM Score'],ascending=[True,True,True,True])
    df['Pep-shift'] = df['Peptide'] +";"+df['Mass shift'].astype(str)
    df['PSM_count'] = df.groupby(['Pep-shift'], sort=False).cumcount() + 1
    USI_list=[]
    for i in df['USI'].groupby(df['Pep-shift']).apply(';'.join):
        USI_list.append(i)
    protein_list=[]
    for i in df['Protein'].groupby(df['Pep-shift']).apply(','.join):
        protein_list_temp = []
        all_protein = ""
        for p in re.split(",|:",i):
            if p not in protein_list_temp:
                protein_list_temp.append(p)
                all_protein+=";"+p
        protein_list.append(all_protein[1:])
    df['Unique'] = np.where(df['PSM_count'] == 1, "TRUE", "FALSE")
    df=df.drop_duplicates(subset=('Pep-shift'),keep='last', inplace=False)
    df['All_Proteins']=protein_list
    df['All_USI']=USI_list
    df.to_csv(mass_based,index=False)
	
	#Collapse by peptide
    df2 = pd.read_csv(output)
    df2 = df2.loc[df2['FDR'] <= FDR_cutoff]
    df2=df2.sort_values(['Peptide','Score','PTM Score'],ascending=[True,True,True])
    df2['Peptide_count'] = df2.groupby(['Peptide'], sort=False).cumcount()+1
    USI_list2 = []
    for i in df2['USI'].groupby(df2['Peptide']).apply(';'.join):
        USI_list2.append(i)
    protein_list2 = []
    for i in df2['Protein'].groupby(df2['Peptide']).apply(','.join):
        protein_list_temp2 = []
        all_protein2 = ""
        for p in re.split(",|:",i):
            if p not in protein_list_temp2:
                protein_list_temp2.append(p)
                all_protein2 += ";" + p
        protein_list2.append(all_protein2[1:])
    df2['Unique'] = np.where(df2['Peptide_count'] == 1, "TRUE", "FALSE")
    df2=df2.drop_duplicates(subset=('Peptide'),keep='last', inplace=False)
    df2['All_Proteins'] = protein_list2
    df2['All_USI']=USI_list2
    df2.to_csv(peptide_based,index=False)
	
	#Collapse by protein position
    df3 = pd.read_csv(output, dtype={'Protein position': str, 'PTM Score':str, 'PTM positions':str})
    df3 = df3.loc[df3['FDR'] <= FDR_cutoff]
    df3.dropna(subset=['PTM'], inplace=True)
    df3.reset_index(drop=True)
    df3['All_PTMs'] = df3['PTM']
    df3['All_PTM_scores'] = df3['PTM Score']
    df3['All_PTM_protein_positions'] = df3['Protein position']
    df3['PTM'] = df3['PTM'].str.split(';')
    df3['PTM Score'] = df3['PTM Score'].str.split(';')
    df3['Protein position'] = df3['Protein position'].str.split(';')
    df3 = explode(df3, ['Protein position','PTM Score', 'PTM'], fill_value='')
    df3 = df3.sort_values(['Protein', 'Protein position', 'Score', 'PTM Score'], ascending=[True, True, True, True])
    df3 = df3.reset_index(drop=True)
    df3['Pro-pos'] = df3['Protein'] + "-" + df3['Protein position'].astype(str)
    df3['Protein_count'] = df3.groupby(['Pro-pos'], sort=False).cumcount()+1
    USI_list3 = []
    for i in df3['USI'].groupby(df3['Pro-pos']).apply(';'.join):
        USI_list3.append(i)
    df3['Unique'] = np.where(df3['Protein_count'] == 1, "TRUE", "FALSE")
    df3=df3.drop_duplicates(subset=('Pro-pos'),keep='last', inplace=False)
    df3['All_Proteins']=df3['Protein']
    df3['All_USI']=USI_list3
    df3.to_csv(site_based,index=False)

	#No collapse
    df4 = pd.read_csv(output)
    df4 = df4.loc[df4['FDR'] <= FDR_cutoff]
    df4 = df4.sort_values(['Peptide', 'Score','PTM Score'], ascending=[True, True, True])
    df4['All_Proteins']=df4['Protein']
    df4['All_USI']=df4['USI']
    df4.to_csv(non_collapse,index=False)



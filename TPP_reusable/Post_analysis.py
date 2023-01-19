#Expands PSMs to site-based format, where one row is one site on a peptide. Calculates FLR calculations.


import pandas as pd
import numpy as np
import sys
import time

def explode(df, lst_cols, fill_value='', preserve_index=False):
    print("\t Running: Explode")
    start_time=time.time()
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
    print("\t Complete --- %s seconds ---" % (time.time() - start_time))

    return res

#Explode rows for each PSM postion - for each mod in peptide, a seperate row will show PSM score and corresponding PTM position/score for each site on the peptide
def site_based(input,FDR_cutoff,mod,verbose):
    print("Running: site_based")
    start_time = time.time()
    if verbose:
        output=input.replace("FDR_output", "Site-based_verbose")
    else:
        output = input.replace("FDR_output", "Site-based")
    df = pd.read_csv(input, dtype={'Protein position': str, 'PTM Score':str, 'PTM positions':str,'FDR':float})
    df = df.loc[df['FDR'] <= float(FDR_cutoff)]
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
    #df = df[df.PTM == "Phospho"]
    df = df.loc[df['PTM'].str.lower()==mod.lower()]
    df = df[~df.Protein.str.contains("DECOY")]
    df = df[~df.Protein.str.contains("CONTAM")]
    df = df.reset_index(drop=True)
    df['Peptide_pos']  = df['Peptide']+"-"+df['PTM positions'].astype(str)
    df['PTM_final_prob'] = pd.to_numeric(df['Score']) * pd.to_numeric(df['PTM Score'])
    df['Peptide']=df['Peptide'].astype(str)
    df = df.sort_values(by=(['PTM_final_prob','Peptide','Peptide_pos']), ascending=[False,True,True])
    df.to_csv(output,index=False)
    print("Complete --- %s seconds ---" % (time.time() - start_time))

#model FLR from combined probability: 1-(PSM prob*PTM prob)/Count
def model_FLR(file,mod,verbose):
    print("Running: Model_FLR")
    start_time = time.time()
    if verbose:
        file=file.replace(".csv","_verbose.csv")
    df=pd.read_csv(file, dtype={'Score': float,'PTM Score':float})
    df = df[df['PTM positions'].notna()]
    #df = df[df.PTM == "Phospho"]
    df = df.loc[df['PTM'].str.lower()==mod.lower()]
    if len(df)==0:
        sys.exit("Modification not found, check spelling or try specifiying modification name and mass to 2dp (eg. Phospho:79.97) \n TPP_comparison.py [mzid_file] [PXD] [modification:target:decoy] [optional: modification:mass(2dp)] [optional: PSM FDR_cutoff]")
    df = df[~df.Protein.str.contains("DECOY",na=False)]
    df = df[~df.Protein.str.contains("CONTAM",na=False)]
    df = df.reset_index(drop=True)
    df['Peptide']=df['Peptide'].astype(str)
    df = df.sort_values(by=(['PTM_final_prob','Peptide','Peptide_pos']), ascending=[False,True,True])
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
    print("Complete --- %s seconds ---" % (time.time() - start_time))

def peptidoform_to_peptide(file,mod,verbose):
    print("Running: peptidoform_to_peptide")
    start_time = time.time()
    if verbose:
        file=file.replace(".csv","_verbose.csv")
    df=pd.read_csv(file)
    df['Peptide_mod_count'] = df['Peptide']+"_"+df['Peptide_mod'].str.count(mod).astype(str)
    df['Prob'] = 1-df['Binomial_final_score']
    df['Prod_prob'] = df.groupby('Protein-pos')['Prob'].transform('prod')
    df['BA_score_new'] = 1-df['Prod_prob']
    df['PepMeanScore'] = df.groupby('Peptide_mod')['BA_score_new'].transform('mean')
    df = df.sort_values(by=(['Peptide_mod_count','PepMeanScore']), ascending=[False,True])
    df_temp = df[df.PepMeanScore == df.PepMeanScore.groupby(df['Peptide_mod_count']).transform('max')]
    #print(len(df_temp))
    peptidoform_list=df_temp['Peptide_mod'].to_list()
    df=df[df['Peptide_mod'].isin(peptidoform_list)]
    df = df.sort_values(by=(['Binomial_final_score','Peptide','Peptide_pos']), ascending=[False,True,True])
    output=file.replace("binomial_peptidoform_collapsed","binomial_peptide_collapsed")

    #return as csv
    if verbose:
        df.to_csv(output,index=False)
    df=df.drop(["Peptide_mod_count","Prob","Prod_prob","BA_score_new","PepMeanScore"], axis=1)
    df.to_csv(output,index=False)

    print("Complete--- %s seconds ---" % (time.time() - start_time))

#calculate decoy FLR from decoy input
def calculate_decoy_FLR(input,decoy,targets, verbose):
    print("Running: calculate_decoy_FLR")
    start_time = time.time()
    if verbose:
        input=input.replace("Site-based_FLR.csv","Site-based_verbose_FLR.csv")
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
        #df.loc[i,'p'+decoy+'_FLR']=(STY_A_ratio*df.loc[i,'p'+decoy+'_count']*2)/(i+1)
        #df.loc[i,'p'+decoy+'_FLR']=((STY_A_ratio*df.loc[i,'p'+decoy+'_count'])+df.loc[i,'p'+decoy+'_count'])/(i+1)
        df.loc[i,'p'+decoy+'_FLR']=(STY_A_ratio*df.loc[i,'p'+decoy+'_count'])/(i+1-df.loc[i,'p'+decoy+'_count'])

    df['p'+decoy+'_q_value'] = df['p'+decoy+'_FLR']
    df['p'+decoy+'_q_value'] = df.iloc[::-1]['p'+decoy+'_FLR'].cummin()

    #return as csv
    if verbose:
        df.to_csv(input,index=False)
    df=df.drop(["Count","DecoyP","pA_count",], axis=1)
    df.to_csv(input,index=False)
    print("Complete --- %s seconds ---" % (time.time() - start_time))
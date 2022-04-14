#Expands PSMs to site-based format, where one row is one site on a peptide. Calculates and generates figures for FLR calculations.
#Output="All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR.csv", can also be used for filtering for no-choice peptides ("All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_filtered.csv") or sorting by PTM score rather than combined probability ("All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_PTM_sort.csv") or collapsing by peptide+mod ("All_confident_PTM_no_collapse_Site-based_spectrum_match_new_FLR_collapse.csv") or collapsed by protein position ("Site_confident_PTM_unique_FLR.csv") 
#These output files are then processed by "pAla" (see below) to provide the decoy amino acid FLR calculations
#Final FLR calculations are then passed through either "plot_FLR_comparisons" or "plot_FLR_comparisons_PTM_prob" functions 
#Output folder is generated "Comparisons/[software names, seperated by "_"]/" and multiple comparison plots, from each of the ordering/filtering options seen above, are generated here


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

#Explode rows for each protein PTM position - for each mod in peptide, a seperate row will show PSM score and corresponding PTM position/score for each site on each protein
def site_based_all(input):
    df = pd.read_csv(input, dtype={'Protein position': str, 'PTM Score':str, 'PTM positions':str})
    df.dropna(subset=['PTM'], inplace=True)
    df['Protein'] = df['Protein'].str.split(":")
    df['Protein position'] = df['Protein position'].str.split(':')
    df = explode(df, ['Protein', 'Protein position'], fill_value='')
    df['PTM'] = df['PTM'].str.split(';')
    df['PTM Score'] = df['PTM Score'].str.split(';')
    df['PTM positions'] = df['PTM positions'].str.split(';')
    df['Protein position'] = df['Protein position'].str.split(';')
    df = explode(df, ['PTM', 'PTM Score', 'PTM positions', 'Protein position'], fill_value='')
    df = df[df.PTM == "Phospho"]
    df=df[~df.Protein.str.contains("DECOY")]
    df=df[~df.Protein.str.contains("CONTAM")]
    df = df.reset_index(drop=True)
    output=input.replace(".csv","_Site-based_all_proteins.csv")
    df.to_csv(output,index=False)

#Explode rows for each PSM postion - for each mod in peptide, a seperate row will show PSM score and corresponding PTM position/score for each site on the peptide
def site_based(input):
    df = pd.read_csv(input, dtype={'Protein position': str, 'PTM Score':str, 'PTM positions':str})
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
    output = input.replace(".csv", "_Site-based.csv")
    df.to_csv(output,index=False)

#Add column for if spectrum matches throughout comparisons
def spectrum_comparison(input,spectrum_match_list):
    df=pd.read_csv(input).fillna('N/A')
    match_counts=[]
    for i in range(len(df)):
        match_count=0
        spectrum=df.loc[i,'Spectrum']
        spectrum=spectrum.replace(".0",".")
        if spectrum in spectrum_match_list:
            match_count=1
        match_counts.append(match_count)
    df['Spectrum_matches']=match_counts
    output=input.replace(".csv","_spectrum_match.csv")
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

#calculate FLR collapsing by peptide+mod and ranking by both final probablity and supporting PSM count 
def model_FLR_new(file):
    df=pd.read_csv(file)
    df = df[df['PTM positions'].notna()]
    df = df[df.PTM == "Phospho"]
    df = df[~df.Protein.str.contains("DECOY",na=False)]
    df = df[~df.Protein.str.contains("CONTAM",na=False)]
    df = df.reset_index(drop=True)
    df['PSM_count'] = df.groupby('Peptide_mod')['Peptide_mod'].transform('count')
    df['PTM_final_prob'] = round(df['Score'] * df['PTM Score'],2)
    df = df.sort_values(by=(['PTM_final_prob','PSM_count']), ascending=[False,False])
    df=df.drop_duplicates(subset=('Peptide_mod'), keep='first', inplace=False)
    df = df.reset_index(drop=True)
    df['Count'] = (df.index) + 1
    df['final_temp'] = 1 - df['PTM_final_prob']
    df['PTM_final_prob_FLR'] = df['final_temp'].cumsum() / df['Count']
    df['PTM_final_prob_q_value'] = df['PTM_final_prob_FLR']
    df['PTM_final_prob_q_value'] = df.iloc[::-1]['PTM_final_prob_FLR'].cummin()
    df = df.drop(['final_temp'], axis=1)
    df = df.reset_index(drop=True)
    output=file.replace(".csv","_new_FLR_collapse.csv")
    df.to_csv(output,index=False)
	
#calculate FLR, removing all "no choice" peptides where the number of phosphosites is equal to the number of potential phosphosites
def model_FLR_filter(file):
    df=pd.read_csv(file)
    df = df[df['PTM positions'].notna()]
    df = df[df.PTM == "Phospho"]
    df = df[~df.Protein.str.contains("DECOY",na=False)]
    df = df[~df.Protein.str.contains("CONTAM",na=False)]
    df = df.reset_index(drop=True)
    before_filter=len(df.index)
    if "pG" in file:
        df = df[df.Peptide.str.count('S|T|Y|G') > df.Peptide_mod.str.count("Phospho")]
    elif "pL" in file:
        df = df[df.Peptide.str.count('S|T|Y|L') > df.Peptide_mod.str.count("Phospho")]
    elif "pD" in file:
        df = df[df.Peptide.str.count('S|T|Y|D') > df.Peptide_mod.str.count("Phospho")]
    elif "pE" in file:
        df = df[df.Peptide.str.count('S|T|Y|E') > df.Peptide_mod.str.count("Phospho")]
    elif "pP" in file:
        df = df[df.Peptide.str.count('S|T|Y|P') > df.Peptide_mod.str.count("Phospho")]
    else:
        df = df[df.Peptide.str.count('S|T|Y|A')>df.Peptide_mod.str.count("Phospho")]
    df = df.reset_index(drop=True)
    after_filter=len(df.index)
    print("Rows removed by no choice filtering: "+str(before_filter-after_filter))
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
    output=file.replace(".csv","_FLR_filtered.csv")
    df.to_csv(output,index=False)

#Plot FLR comparisons - pX_FLR and model FLR for each software
def plot_FLR_comparisons(filter,file_list,software_list,working):
    import seaborn as sns;sns.set_palette("colorblind")
    folder=""
    legend=[]
    for i,j in zip(software_list,file_list):
        folder+="_"+i
        if "Site_confident" not in j:
            legend.append(i + " Model FLR")
        if "pL" in j or "pLeu" in j:
            legend.append(i + " pLeu Decoy FLR")
        elif "pGlu" not in j and ("pG" in j or "pGly" in j):
            legend.append(i + " pGly Decoy FLR")
        elif "pD" in j or "pAsp" in j:
            legend.append(i + " pAsp Decoy FLR")
        elif "pE" in j or "pGlu" in j:
            legend.append(i + " pGlu Decoy FLR")
        elif "pP" in j or "pPro" in j:
            legend.append(i + " pPro Decoy FLR")
        else:
            legend.append(i + " pAla Decoy FLR")
    output=working+"/Comparisons/"+folder[1:]
    if not os.path.exists(output):
        os.mkdir(output)
    print("OUTPUT= "+output)
    #Combined prob FLR plot
    x_list=['PTM_final_prob','Count']
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
            if df[x].max()>max:
                max=df[x].max()
            if x == "PTM_final_prob":
                lim = (1, 0)
            else:
                lim=(0,max)
            c = "C" + str(counter + 1)
            if "Site_confident" not in file_list[0]:
                if counter == 0:
                    ax = df.plot.line(x=x, y="PTM_final_prob_q_value", linestyle=":", color=c, xlim=lim,figsize=(14, 7))
                    df.plot.line(x=x, y=p, linestyle="-", color=c, xlim=lim, ax=ax)
                else:
                    df.plot.line(x=x, y="PTM_final_prob_q_value", linestyle=":", color=c, ax=ax, xlim=lim)
                    df.plot.line(x=x, y=p, linestyle="-", ax=ax, color=c, xlim=lim)
                counter += 1
            else:
                if counter == 0:
                    ax = df.plot.line(x=x, y=p, linestyle="-", color=c, xlim=lim, figsize=(14, 7))
                else:
                    df.plot.line(x=x, y=p, linestyle="-", ax=ax, color=c, xlim=lim)
                counter += 1
        ax.legend(legend)
        ax.set_ylabel("FLR")
        ax.set_xlabel("Count of Sites")
        if x=="PTM_final_prob":
            ax.set_xlabel("Combined probability")
        if "Site_confident_PTM_unique" in file_list[0]:
            if "new_FLR_collapse" in file_list[0]:
                plt.savefig(output + "/" + filter + "_" + x + "_All_software_Site_Confident_new_collapse" + p + ".jpg", dpi=300)
            else:
                plt.savefig(output+"/"+filter+"_"+x+"_All_software_Site_Confident_"+p+".jpg",dpi=300)
        else:
            if "new_FLR_collapse" in file_list[0]:
                plt.savefig(output + "/" + filter + "_" + x + "_All_software_new_collapse_"+p+".jpg", dpi=300)
            elif "filtered" in file_list[0]:
                plt.savefig(output + "/" + filter + "_" + x + "_All_software_filtered_no_choice" + p + ".jpg", dpi=300)
            else:
                plt.savefig(output + "/" + filter + "_" + x + "_All_software_" + p + ".jpg", dpi=300)

#Sort on PTM score
def PTM_sort_FLR(input):
    df=pd.read_csv(input)
    df = df[df['PTM positions'].notna()]
    df = df[df.PTM == "Phospho"]
    df = df[~df.Protein.str.contains("DECOY",na=False)]
    df = df[~df.Protein.str.contains("CONTAM",na=False)]
    df = df.reset_index(drop=True)
    df = df.sort_values(by=(['PTM Score']), ascending=[False])
    df = df.reset_index(drop=True)
    df['Count'] = (df.index) + 1
    output=input.replace(".csv","_FLR_PTM_sort.csv")
    df.to_csv(output,index=False)

def plot_FLR_comparisons_PTM_prob(filter,file_list,software_list,output_folder):
    import seaborn as sns; sns.set_palette("colorblind")
    folder=""
    legend=[]
    for i,j in zip(software_list,file_list):
        folder+="_"+i
        if "pL" in j or "pLeu" in j:
            legend.append(i + " pLeu Decoy FLR")
        elif "pG" in j or "pGly" in j:
            legend.append(i + " pGly Decoy FLR")
        elif "pD" in j or "pAsp" in j:
            legend.append(i + " pAsp Decoy FLR")
        elif "pE" in j or "pGlu" in j:
            legend.append(i + " pGlu Decoy FLR")
        elif "pP" in j or "pPro" in j:
            legend.append(i + " pPro Decoy FLR")
        else:
            legend.append(i + " pAla Decoy FLR")
    output=output_folder+"/Comparisons/"+folder[1:]
    if not os.path.exists(output):
        os.mkdir(output)
    x_list = ['PTM Score','Count']
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

            if df[x].max()>max:
                max=df[x].max()
            if x == "PTM Score":
                lim = (1, 0)
            else:
                lim=(0,max)
            c = "C" + str(counter + 1)
            if counter==0:
                ax=df.plot.line(x=x, y=p, linestyle="-", color=c, xlim=lim, figsize=(14,7))
            else:
                df.plot.line(x=x, y=p, linestyle="-", color=c, ax=ax, xlim=lim)
            counter += 1
        ax.legend(legend)
        ax.set_ylabel("FLR")
        ax.set_xlabel("Count of Sites")
        if x=="PTM Score":
            ax.set_xlabel("PTM probability")
        plt.savefig(output+"/"+filter+"_"+x+"_All_software_"+p+"_PTM.jpg",dpi=300)
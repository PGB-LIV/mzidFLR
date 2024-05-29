#Calculate false detection rate statistics and plots
#	Loads CSV file to dataframe for FDR calculations
#Writes folder for downstream analysis eg. "FDR_0.01_PTM_score_0" FDR filter 0.01, no PTM score filter
#output="FDR_output.csv","FDR.jpg","FDR_filter.jpg","FDR_score.jpg"


import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import time
import sys

f = lambda x, y : "" if y=="NA" else round(float([z if z != '' else '0' for z in x][y]),2)
r = lambda x, y : x.replace("unknown_mod",str(y),1)

#Extract PTM prophet output to df
def extract_PTMprophet_IDent_df(input,PXD,mod, mod_id, mod_mass_id):
    print("Running: extract_PTMprophet_IDent_df")
    start_time = time.time()
    #start_time=time.time()
    all_peptide_mods=[]
    PTMs=[]
    all_positions=[]
    PTM_scores=[]
    mass_shift=[]
    protein_positions=[]
    df = pd.read_csv(input)
    df['PTM info']=df['PTM info'].fillna("")
    df['Modifications']=df['Modifications'].astype(str)
    df['Positions']=df['Positions'].astype(str)
    df['Modification mass']=df['Modification mass'].astype(str)
    df['Protein position']=df['Protein position'].astype(str)

    while df['Modifications'].str.contains("unknown_mod").any():
        df['Mods_temp']=df['Modifications'].str.split(";")
        df['Mod_mass_temp']=df['Modification mass'].str.split(";")
        df['Index_unknown']=df['Mods_temp'].apply(lambda x: x.index("unknown_mod") if "unknown_mod" in x else "NA")
        df['Mods_replace']=df.apply(lambda x: f(x.Mod_mass_temp, x.Index_unknown), axis=1)
        #replace specified mod mass with specified mod, round to 2dp
        df['Modifications']=df.apply(lambda x: r(x.Modifications, x.Mods_replace), axis=1)

    if mod_id!="NA" and mod_mass_id!="NA":
        df['Modifications']=df['Modifications'].str.replace(str(round(float(mod_mass_id),2)),mod_id, regex=True)

    for i in range(len(df)):
        peptide_temp=df.loc[i,'Peptide']
        if df.loc[i,'PTM info']!="":
            PTMs.append(df.loc[i,'Modifications'][:-1])
            all_positions.append(df.loc[i,'Positions'][:-1])
            for a,b in zip(df.loc[i,'Positions'].split(";")[-2::-1],df.loc[i,'Modifications'].split(";")[-2::-1]):
                peptide_temp=peptide_temp[:int(a)]+"["+b+"]"+peptide_temp[int(a):]
            protein_start=df.loc[i,'Protein position']
            protein_position_list=""
            for prot in protein_start.split(":"):
                for p in df.loc[i,'Positions'].split(";")[:-1]:
                    protein_position_list+=str(int(prot)+int(p)-1)+";"
                protein_position_list=protein_position_list[:-1]
                protein_position_list+=":"
            protein_positions.append(protein_position_list[:-1])
        else:
            all_positions.append("")
            PTMs.append("")
            protein_positions.append("")

        #mass_shift.append(mod_mass)
        all_peptide_mods.append(peptide_temp)
        PTMscore_list=""
        if df.loc[i, 'PTM info'] != "":
            for x,z in zip(df.loc[i,'Positions'].split(";")[:-1],df.loc[i,'Modifications'].split(";")[:-1]):
                if z.lower()==mod.lower():
                    score_found="No"
                    for y in df.loc[i,'PTM info'].split(";"):
                        if y.split(":")[2]==x:
                            PTMscore_list+=y.split(":")[1]+";"
                            score_found="Yes"
                    if score_found=="No":
                        PTMscore_list+="0;"
                        if z!="Acetyl":
                            print("CHECK: "+ df.loc[i,"Spectrum"])
                else:
                    PTMscore_list += "0;"
        PTM_scores.append(PTMscore_list[:-1])

        #df.loc[i,'USI']="mzspec:" + PXD + ":" + df.loc[i,'Spectrum'].split(".")[0]  + ":scan:" + df.loc[i,'Spectrum'].split(".")[1] + ":" + peptide_temp + "/" + str(df.loc[i,'Charge'])
        #df.loc[i,'Sources'] = df.loc[i,'Spectrum'].split(".")[0]
    #print("--- %s seconds ---" % (time.time() - start_time))
    df['Modification mass']=df['Modification mass'].str.replace(";;;",";")
    df['Modification mass']=df['Modification mass'].str.replace(";;",";")
    df['Modification mass']=df['Modification mass'].str.lstrip(";")

    df['Mass_shift']=df['Modification mass'].apply(lambda x: sum(list(map(float,x.split(";")[:-1]))) if x!=";" and x!="" else "")
    df['mass_diff'] = df['Calculated mass']-df['Experimental mass']
    df['ppm_error'] = (((df['mass_diff']) / df['Calculated mass']) * 1e6)
    df2 = pd.DataFrame({"Peptide_mod": all_peptide_mods, "Peptide":df['Peptide'].values , "Protein":df['Protein'].values, "Score":df['PSM probability'].values, "PTM":PTMs ,
                        "PTM positions":all_positions,"PTM Score": PTM_scores, "PTM_info":df['PTM info'].values,"Spectrum":df['Spectrum'].values,"mass_diff":df['mass_diff'].values
                           ,"ppm_error":df['ppm_error'].values, "Protein position":protein_positions,"Mass shift":df['Mass_shift'].values,"Charge":df['Charge'].values})

    df2["USI"]="mzspec:"+PXD+":"+df2['Spectrum'].str.rsplit(".",n=3).str[0] +":scan:"+df2['Spectrum'].str.rsplit(".",n=2).str[1]+":"+df2['Peptide_mod']+"/"+df2['Charge'].astype(str)
    df2['Source'] = df2['Spectrum'].str.split(".").str[0]
    #replace weird mods which cause issues with USI validation
    df2=df2.replace("pryo","pyro", regex=True)
    df2=df2.replace("Q\[Pyro_glu\]","Q[Gln->pyro-Glu]", regex=True)
    df2=df2.replace("E\[Pyro_glu\]","E[Glu->pyro-Glu]", regex=True)
    df2=df2.replace("C\[Pyro_glu\]","C[Ammonia-loss]", regex=True)
    df2=df2.replace("Carbamidomethylation","Carbamidomethyl", regex=True)
    df2=df2.replace("346\.21","iTRAQ8plex", regex=True)
    df2=df2.replace("304\.2","iTRAQ8plex", regex=True)
    df2=df2.replace("\[229\.16\]","[TMT6plex]", regex=True)
    df2=df2.replace("42\.01","Acetyl", regex=True)


    return(df2)
    print("complete --- %s seconds ---" % (time.time() - start_time))

#calculate FDR
def calculateFDR(results_file,output,PXD, decoy_prefix, mod, mod_id, mod_mass, verbose):
    print("Running: CalculateFDR")
    start_time = time.time()
    #extract results to df
    df=extract_PTMprophet_IDent_df(results_file,PXD,mod, mod_id, mod_mass)
    df['Score']=df['Score'].astype(float)
    df=df.sort_values(by='Score',ascending=False)
    df=df.reset_index(drop=True)
    #set decoy and target
    df['Protein_count'] = df['Protein'].str.count(":")+1
    df['Decoy_count'] = df['Protein'].str.count(decoy_prefix)
    df['Decoy'] = np.where(df['Protein_count']==df['Decoy_count'],1,0)
    df=df.drop(columns=['Protein_count', 'Decoy_count'])
    #target_count = row_count-decoy_count
    df['decoy_count']=df['Decoy'].cumsum()
    df['row_count']=(df.index)+1
    df['target_count']=df['row_count']-df['decoy_count']
    #FDR = decoy_count/target_count
    df['FDR']= df['decoy_count']/df['target_count']
    df['FDR'] = df['FDR'].astype(float)

    #warning for no decoys
    if df['FDR'].max==0:
        sys.exit('No decoys found - check decoy protein prefix, if not given, default is "DECOY"')

    #q_val = min FDR at this position or lower
    df['q_value']=df['FDR']
    df['q_value']=df.iloc[::-1]['FDR'].cummin()

    #FDR plots
    ax1=df.plot.scatter(x='FDR',y='row_count')
    ax2=df.plot.line(x='q_value',y='row_count')
    plt.ylabel("PSM count")
    plt.savefig('FDR.jpg',dpi=300)

    df.plot.line(x='Score', y='FDR',ylim=(0,1),xlim=(1,0))
    plt.xlabel("Peptide Probabilty")
    plt.ylabel("Global FDR")
    plt.savefig('FDR_score.jpg',dpi=300)
    plt.close('all')

    #return as csv
    if verbose:
        verbose_out=output.replace(".csv","_verbose.csv")
        df.to_csv(verbose_out,index=False)
    df=df.drop(["Decoy","decoy_count","row_count","target_count"], axis=1)
    df.to_csv(output,index=False)

    print("Complete --- %s seconds ---" % (time.time() - start_time))


#precursor mass tolerance plot
def ppm_error(results_file):
    print("Running: ppm_error")
    start_time = time.time()
    df = pd.read_csv(results_file)

    df['Error_round']=round(df['ppm_error'],0)
    df['ppm_error']=df['ppm_error']-(df['Error_round']*1.00866491595)

    df=df.loc[df['mass_diff'] <=0.9]
    df.plot.scatter(x='Score', y='ppm_error', s=0.5)
    plt.ylabel("PPM error")
    plt.xlabel("PSM probability")
    plt.savefig('PPM_error.jpg', dpi=300)
    df.to_csv("error.csv", index=False)

    df = df.loc[df['q_value'] <= 0.01]
    df=df.loc[df['mass_diff'] <=0.9]
    df.plot.scatter(x='Score', y='ppm_error', s=0.5)
    plt.ylabel("PPM error")
    plt.xlabel("PSM probability")
    plt.savefig('PPM_error_FDR.jpg', dpi=300)
    df.to_csv("error_FDR.csv", index=False)
    print("Complete --- %s seconds ---" % (time.time() - start_time))

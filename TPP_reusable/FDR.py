#Calculate false detection rate statistics and plots
#-Calls Extract_DF
#	Loads CSV file to dataframe for FDR calculations
#Writes folder for downstream analysis eg. "FDR_0.01_PTM_score_0" FDR filter 0.01, no PTM score filter
#output="FDR_output.csv","FDR.jpg","FDR_filter.jpg","FDR_score.jpg"

	
import TPP_reusable.Extract_DF as Extract_DF
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import time

#Extract PTM prophet output to df
def extract_PTMprophet_IDent_df(input,PXD,mod):
    start_time=time.time()
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
    counter=1
    for i in range(len(df)):
        if "unknown" in df.loc[i,'Modifications']:
            mods_temp=df.loc[i,'Modifications'].split(";")
            index=mods_temp.index('unknown_mod')
            mods_temp[index]=str(df.loc[i, 'Modification mass'].split(";")[index])
            df.loc[i,'Modifications']=";".join([str(x)for x in mods_temp])
        peptide_temp=df.loc[i,'Peptide']
        mod_mass=0
        if df.loc[i,'PTM info']!="":
            PTMs.append(df.loc[i,'Modifications'][:-1])
            all_positions.append(df.loc[i,'Positions'][:-1])
            for a,b in zip(df.loc[i,'Positions'].split(";")[-2::-1],df.loc[i,'Modifications'].split(";")[-2::-1]):
                peptide_temp=peptide_temp[:int(a)]+"["+b+"]"+peptide_temp[int(a):]
            for m in df.loc[i,'Modification mass'].split(";")[:-1]:
                if m!="":
                    mod_mass+=float(m)
            protein_start=df.loc[i,'Protein position']
            protein_position_list=""
            for p in df.loc[i,'Positions'].split(";")[:-1]:
                protein_position_list+=str(protein_start+int(p)-1)+";"
            protein_positions.append(protein_position_list[:-1])
        else:
            all_positions.append("")
            PTMs.append("")
            protein_positions.append("")

        mass_shift.append(mod_mass)
        all_peptide_mods.append(peptide_temp)
        PTMscore_list=""
        if df.loc[i, 'PTM info'] != "":
            for x,z in zip(df.loc[i,'Positions'].split(";")[:-1],df.loc[i,'Modifications'].split(";")[:-1]):
                #if z=="Carbamidomethylation" or z=="Acetyl" or z=="Glu->pyro-Glu" or z=="unknown_mod" or z=="Deamidated" or z=="Gln->pryo-Glu" or "iTRAQ" in z or "." in z:
                #PTMscore_list+="0;"
                #else:
                #for y in df.loc[i,'PTM info'].split(";"):
                #if y.split(":")[2]==x:
                #PTMscore_list+=y.split(":")[1]+";"
                if z.lower()==mod.lower():
                    score_found="No"
                    for y in df.loc[i,'PTM info'].split(";"):
                        if y.split(":")[2]==x:
                            PTMscore_list+=y.split(":")[1]+";"
                            score_found="Yes"
                    if score_found=="No":
                        PTMscore_list+="0;"
                        print("CHECK: "+ df.loc[i,"Spectrum"])
                else:
                    PTMscore_list += "0;"
        PTM_scores.append(PTMscore_list[:-1])

        df.loc[i,'USI']="mzspec:" + PXD + ":" + df.loc[i,'Spectrum'].split(".")[0]  + ":scan:" + df.loc[i,'Spectrum'].split(".")[1] + ":" + peptide_temp + "/" + str(df.loc[i,'Charge'])
        df.loc[i,'Sources'] = df.loc[i,'Spectrum'].split(".")[0]
        #print(str(counter)+"/"+str(len(df)))
        counter+=1
    #print("--- %s seconds ---" % (time.time() - start_time))

    df['mass_diff'] = df['Calculated mass']-df['Experimental mass']
    df['ppm_error'] = (((df['mass_diff']) / df['Calculated mass']) * 1e6)
    df2 = pd.DataFrame({"Peptide_mod": all_peptide_mods, "Peptide":df['Peptide'].values , "Protein":df['Protein'].values, "Source":df['Sources'].values, "Score":df['PSM probability'].values, "PTM":PTMs ,
                        "PTM positions":all_positions,"PTM Score": PTM_scores, "PTM_info":df['PTM info'].values,"Spectrum":df['Spectrum'].values, "USI": df['USI'].values,"mass_diff":df['mass_diff'].values
                           ,"ppm_error":df['ppm_error'].values, "Protein position":protein_positions,"Mass shift":mass_shift})
    return(df2)

#calculate FDR
def calculateFDR(results_file,output,PXD,mod):
    #extract results to df
    df=extract_PTMprophet_IDent_df(results_file,PXD,mod)
    df['Score']=df['Score'].astype(float)
    df=df.sort_values(by='Score',ascending=False)
    df=df.reset_index(drop=True)
    #set decoy and target
    df['Protein_count'] = df['Protein'].str.count(":")+1
    df['Decoy_count'] = df['Protein'].str.count("DECOY")
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
        sys.exit('No decoys found - check decoy protein prefix contains "DECOY"')

    #q_val = min FDR at this position or lower
    df['q_value']=df['FDR']
    df['q_value']=df.iloc[::-1]['FDR'].cummin()

    #return as csv
    df.to_csv(output,index=False)

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


#precursor mass tolerance plot
def ppm_error(results_file):
    df = pd.read_csv(results_file)
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
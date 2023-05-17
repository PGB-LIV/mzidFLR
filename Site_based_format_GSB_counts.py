import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import os

#read in csv from FLR pipeline
#DECOY and CONTAM removed, filtered for contains Phospho, exploded for site based, binomial adjustment and collapsed by peptidoform

flr_filter=0.05

wd="C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/"
folder_list=["PXD002222","PXD004939","PXD005241","PXD004705","PXD002756","PXD012764","PXD000923","PXD019291"]
meta="C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/rice_ptm_metadata.csv"

with open(meta,'r') as f:
    reader = csv.reader(f)
    dict = {rows[0]:rows[1:] for rows in reader}
f.close()

for i in folder_list:
    file=wd+i+"/FDR_updated_0.01/binomial_peptidoform_collapsed_FLR.csv"
    df_temp = pd.read_csv(file)
    search_mod="Phosphorylation"
    #Filter for decoy FLR? - not needed due to "Site Passes Threshold" column
    df = df_temp[['All_Proteins','Peptide','USI','pA_q_value_BA','PTM Score','All_PTM_positions','Score','PTM positions','PTM_info','Binomial_final_score','Peptide_mod','All_PTMs','All_PTM_positions',
                  'All_PTM_scores','All_USI','Protein position','0.05FLR_threshold_count','0.01<P<=0.05_count','0.05<P<=0.19_count','0.19<P<=0.81_count', '0.81<P<=0.95_count', '0.95<P<0.99_count', 'P>=0.99_count']].copy()
    df = df.rename(columns={'All_Proteins':'Proteins','pA_q_value_BA':'Binomial adjusted pA q value','Score':'PSM Score','PTM Score':'Site probability','All_USI':'Supporting PSM count',
                            'PTM positions':'Site position','Binomial_final_score':'Final site probability','All_PTMs':'Modifications','All_PTM_scores':'Modification probabilities',
                            'All_PTM_positions':'Modification positions','0.05FLR_threshold_count':'Supporting PSM Count (Thresholded 0.05 Global FLR)','0.01<P<=0.05_count':'PSM count 0.01<P<=0.05','0.05<P<=0.19_count':'PSM count 0.05<P<=0.19',
                            '0.19<P<=0.81_count':'PSM count 0.19<P<=0.81', '0.81<P<=0.95_count':'PSM count 0.81<P<=0.95', '0.95<P<0.99_count':'PSM count 0.95<P<0.99', 'P>=0.99_count':'PSM count P>=0.99'})

    df['Proteins'] = df['Proteins'].str.replace(":",";")
    dataset_ID=[]
    #If no choice peptide, delete row
    df = df[df.Peptide.str.count('S|T|Y|A')>df.Peptide_mod.str.count("Phospho")]
    #Remove peptidoforms with pA
    df=df[~df['Peptide_mod'].str.contains("A\[Phosphorylation\]")==True]
    df=df[~df['Peptide_mod'].str.contains("A\[Phospho\]")==True]
    df['Modification']=search_mod

    df = df.rename(columns={'Peptide_mod':'Modified Peptide','Peptide':'Unmodified Sequence'})
    df = df.reset_index(drop=True)

    ######
    df['PRIDE dataset URL']=""
    df['PubMedID']=""
    df['Sample ID']=""
    df['PRIDE dataset URL (including the USI)']=""
    df['Supporting PSM count']=df['Supporting PSM count'].str.count(";")+1

    df['Pep_pos'] = df['Unmodified Sequence'] + "-" + df['Site position'].astype(str)
    organism=[]
    PMID=[]
    for i in range(len(df)):
        USI_temp=df.loc[i,'USI']
        PXD=USI_temp.split(":")[1]
        dataset_ID.append(PXD)
        organism.append(dict[PXD][1])
        PMID.append(dict[PXD][0])
    df['Site Passes Threshold']=np.where(df['Binomial adjusted pA q value'] <= flr_filter, 1, 0)
    df['Dataset ID']=dataset_ID
    df['Site ID']=df.index
    df['Organism']=organism
    df['PubMedID']=PMID

    df = df[['Site ID', 'Proteins', 'Unmodified Sequence', 'Modified Peptide','Modification','Site position', 'Site probability', 'PSM Score',
             'Final site probability','Binomial adjusted pA q value','Site Passes Threshold', 'Supporting PSM count','Supporting PSM Count (Thresholded 0.05 Global FLR)', 'Dataset ID', 'PRIDE dataset URL', 'PubMedID',
             'Sample ID', 'Organism', 'USI','PRIDE dataset URL (including the USI)','PSM count 0.01<P<=0.05','PSM count 0.05<P<=0.19','PSM count 0.19<P<=0.81','PSM count 0.81<P<=0.95','PSM count 0.95<P<0.99','PSM count P>=0.99']]
    if not os.path.exists(wd+"Site_format_files"):
        os.mkdir(wd+"Site_format_files")
    df.to_csv(wd+"Site_format_files/"+dataset_ID[0]+"_Site_peptidoform_centric_binomial_FLR_no_pA_bins_"+str(flr_filter)+".tsv", sep="\t", index=False)


for loc in folder_list:
    loc_full=wd+"/"+loc+"/FDR_updated_0.01/binomial_peptidoform_collapsed_FLR.csv"
    df=pd.read_csv(loc_full)
    df=df.loc[df['pA_q_value_BA']<=flr_filter]
    df=df.sort_values(['Protein-pos','Binomial_final_score','pA_q_value_BA'],ascending=[True,True,False])
    df['PTM_residue']=df.apply(lambda x: x['Peptide'][x['PTM positions']-1],axis=1)
    df['Protein_pos_res']=df['Protein-pos']+"_"+df['PTM_residue']
    df=df.drop_duplicates(subset=('Protein_pos_res'),keep="last",inplace=False)
    df=df[['Peptide_mod_pos', 'pA_q_value_BA','Binomial_final_score', 'Protein_pos_res']]
    df.rename(columns = {'pA_q_value_BA':loc+"_FLR",'Binomial_final_score':loc+'_BinomialScore','Peptide_mod_pos':loc+"_peptide_mod_pos"}, inplace = True)
    df=df.set_index('Protein_pos_res')

    if loc==folder_list[0]:
        df_final=df
    else:
        df_final=pd.concat([df_final,df],axis=1)

for i in df_final.index.values.tolist():
    df_final.loc[i,'PTM_residue']=i.rsplit("_")[-1]
    flr1_count=0
    for dataset in folder_list:
        FLR_count_all=0
        FLR_cols = [x for x in folder_list if dataset in x]
        for FLR_col in FLR_cols:
            if df_final.loc[i,FLR_col+"_FLR"]!="N/A":
                if float(df_final.loc[i,FLR_col+"_FLR"])<=0.01:
                    FLR_count_all+=1
        if FLR_count_all!=0:
            flr1_count+=1
    if flr1_count>1:
        df_final.loc[i,'PTM_FLR_category']="Gold"
    elif flr1_count==1:
        df_final.loc[i,'PTM_FLR_category']="Silver"
    else:
        df_final.loc[i,'PTM_FLR_category']="Bronze"

df_final.to_csv(wd+"Site_format_files/GSB_"+str(flr_filter)+"_protein_pos.csv")
counts=pd.crosstab(df_final.PTM_residue,df_final.PTM_FLR_category).replace(0,np.nan).stack().reset_index().rename(columns={0:'Count of Sites'})
fig, axes = plt.subplots(ncols=3)
counts['PTM_FLR_category'] = pd.Categorical(counts['PTM_FLR_category'], ["Bronze", "Silver", "Gold"])
counts.sort_values("PTM_FLR_category")
for i, (name, group) in enumerate(counts.groupby("PTM_FLR_category")):
    axes[i].set_title(name)
    group.plot(kind="bar", x = "PTM_residue", y="Count of Sites", ax=axes[i], legend=False)
    axes[i].set_ylabel("count")
    axes[i].set_xlabel("")
plt.tight_layout()
#plt.show()
plt.savefig(wd+"Site_format_files/GSB_"+str(flr_filter)+"_protein_pos_categories.png",dpi=300)

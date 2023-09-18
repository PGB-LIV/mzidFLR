import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import requests
import json

#read in csv from FLR pipeline
#DECOY and CONTAM removed, filtered for contains Phospho, exploded for site based, binomial adjustment and collapsed by peptidoform

#bins=[-1,0.01,0.05,1]
#labels=["less than 1%FLR","Between 1% and 5%FLR","Does not pass 5% FLR threshold"]

wd="C:/users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/"
folder_list=["PXD002222","PXD004939","PXD005241","PXD004705","PXD002756","PXD012764","PXD000923","PXD019291"]

meta_all="C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/rice_ptm_metadata.csv"
with open(meta_all,'r') as f:
    reader = csv.reader(f)
    dict_all = {rows[0]:rows[1:] for rows in reader}
f.close()


if not os.path.exists(wd+"/All_site_formats"):
    os.mkdir(sub)

r = lambda x, y : 0 if x[int(y)-1]!="A" else 1

### Peptide site format
for i in folder_list:
    print(i)
    file=wd+i+"/FDR_0.01/binomial.csv"
    df_temp = pd.read_csv(file)
    search_mod="Phospho"
    #Filter for decoy FLR? - not needed due to "Site Passes Threshold" column
    df = df_temp[['All_Proteins','Peptide','USI','pA_q_value','PTM Score','All_PTM_positions','Score','PTM positions','PTM_info','Binomial_final_score','Peptide_mod','All_PTMs','All_PTM_positions',
                  'All_PTM_scores','All_USI','Protein position','Source']].copy()
    df = df.rename(columns={'All_Proteins':'Proteins','pA_q_value':'Site Q-Value','Score':'PSM Probability','PTM Score':'PTM Probability',
                            'PTM positions':'Peptide Modification Position','Binomial_final_score':'Final Probability','All_PTMs':'Modifications','All_PTM_scores':'Modification probabilities',
                            'All_PTM_positions':'Modification positions','Peptide_mod':'Peptidoform','Peptide':'Unmodified Sequence','USI':'Universal Spectrum Identifier'}) #'All_USI':'Supporting PSM count',
    df = df.reset_index(drop=True)
    df['Proteins'] = df['Proteins'].str.replace(":",";")
    df['Modification']=search_mod

    ######
    meta_file=wd+"PRIDE/SDRFs/"+i+".sdrf.tsv"
    meta=pd.read_csv(meta_file,sep="\t")
    meta['comment[data file]']=meta['comment[data file]'].str.split(".").str[0]
    meta['comment[data file]']=meta['comment[data file]'].str.replace(" ","_")
    meta.set_index('comment[data file]', inplace=True)
    dict=meta.to_dict()


    PMID=[]
    Sample=[]
    dataset_ID=[]
    organism=[]
    organism_part=[]
    cell_line=[]
    disease=[]
    for x in range(len(df)):
        source=df.loc[x,"Source"]
        source=source.replace("_raw","")
        PXD=dict["comment[proteomexchange accession number]"][source]
        PMID.append(dict_all[PXD][0])
        Sample.append(dict["source name"][source])
        dataset_ID.append(PXD)
        organism.append(dict["characteristics[organism]"][source])
        organism_part.append(dict["characteristics[organism part]"][source])
        cell_line.append(dict["characteristics[cell type]"][source])
        disease.append(dict["characteristics[disease]"][source])
    df['PubMedIDs']=PMID
    df['Sample ID']=Sample
    df['Source Dataset Identifier']=dataset_ID
    df['Reanalysis Dataset Identifier']="NA"
    df['Organism']=organism
    df['Organism Part']=organism_part
    df['Cell Line']=cell_line
    df['Disease Information']=disease
    #####


    df['Pep_pos'] = df['Unmodified Sequence'] + "-" + df['Peptide Modification Position'].astype(str)
    df['Site Passes Threshold [0.05]'] = np.where(df['Site Q-Value'] <= 0.05, 1, 0)
    df['Site Passes Threshold [0.01]'] = np.where(df['Site Q-Value'] <= 0.01, 1, 0)
    df['Decoy Peptide']=np.where(df['Proteins'].str.contains("DECOY")==True, 1, 0)
    df['Decoy Modification Site']=df.apply(lambda x: r(x.Pep_pos.split("-")[0], x.Pep_pos.split("-")[1]), axis=1)
    df['PSM Site ID']=df.index
    df['PSM Count Passing Threshold 0.05']=1


    #peptidoform without reagent labels?
    df['Peptidoform']=df['Peptidoform'].str.replace(r'\[TMT(.*?)\]','',regex="True")
    df['Peptidoform']=df['Peptidoform'].str.replace(r'\[iTRAQ.*?\]','',regex="True")


    df = df[['PSM Site ID', 'Proteins', 'Unmodified Sequence', 'Peptidoform','Modification','Peptide Modification Position', 'PSM Probability','PTM Probability',
             'Final Probability','Site Q-Value','Site Passes Threshold [0.05]','Site Passes Threshold [0.01]', 'Decoy Peptide','Decoy Modification Site',
             'PSM Count Passing Threshold 0.05','Source Dataset Identifier', 'Reanalysis Dataset Identifier', 'PubMedIDs','Sample ID', 'Organism', 'Organism Part',
             'Cell Line','Disease Information', 'Universal Spectrum Identifier']]
    df = df.replace("not applicable","NA")

    df.to_csv(wd+"All_site_formats/"+dataset_ID[0]+"_Site_PSM_centric.tsv", sep="\t", index=False)


    df=df.head(100)
    #Check USI valid
    url = 'https://proteomecentral.proteomexchange.org/api/proxi/v0.1/usi_validator'
    x = requests.post(url, json=df['Universal Spectrum Identifier'].to_list())
    res=json.loads(x.text)

    for i in range(len(df)):
        #df.loc[i,'Valid_USI']=res['validation_results'][df.loc[i,'USI']]['is_valid']
        if res['validation_results'][df.loc[i,'Universal Spectrum Identifier']]['is_valid']==False:
            print(df.loc[i,'Universal Spectrum Identifier'])

### Peptidoform site format

print("PSM format complete, creating Peptidoform format")
for i in folder_list:
    print(i)
    file=wd+i+"/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
    df_temp = pd.read_csv(file)
    search_mod="Phospho"
    #Filter for decoy FLR? - not needed due to "Site Passes Threshold" column
    df = df_temp[['All_Proteins','Peptide','USI','pA_q_value_BA','PTM Score','All_PTM_positions','Score','PTM positions','PTM_info','Binomial_final_score','Peptide_mod','All_PTMs','All_PTM_positions',
                  'All_PTM_scores','All_USI','Protein position','0.05FLR_threshold_count','0.01FLR_threshold_count','0.01<P<=0.05_count','0.05<P<=0.19_count','0.19<P<=0.81_count', '0.81<P<=0.95_count',
                  '0.95<P<0.99_count', 'P>=0.99_count','All_Source']].copy()
    df = df.rename(columns={'All_Proteins':'Proteins','Peptide_mod':'Peptidoform','Peptide':'Unmodified Sequence','pA_q_value_BA':'Site Q-Value','Score':'PSM Probability','PTM Score':'PTM Probability',
                            'PTM positions':'Peptide Modification Position','Binomial_final_score':'Final Probability','All_PTMs':'Modifications','All_PTM_scores':'Modification probabilities',
                            'All_PTM_positions':'Modification positions','0.05FLR_threshold_count':'PSM Count Passing Threshold [0.05]','0.01FLR_threshold_count':'PSM Count Passing Threshold [0.01]',
                            '0.01<P<=0.05_count':'opt_PSM count 0.01<P<=0.05','0.05<P<=0.19_count':'opt_PSM count 0.05<P<=0.19','0.19<P<=0.81_count':'opt_PSM count 0.19<P<=0.81',
                            '0.81<P<=0.95_count':'opt_PSM count 0.81<P<=0.95', '0.95<P<0.99_count':'opt_PSM count 0.95<P<0.99', 'P>=0.99_count':'opt_PSM count P>=0.99',
                            'USI':'Universal Spectrum Identifier'}) #'All_USI':'Supporting PSM count'
    df = df.reset_index(drop=True)
    df['Proteins'] = df['Proteins'].str.replace(":",";")
    dataset_ID=[]
    df['Modification']=search_mod
    #df['Supporting PSM count']=df['Supporting PSM count'].str.count(";")+1

    df['Pep_pos'] = df['Unmodified Sequence'] + "-" + df['Peptide Modification Position'].astype(str)
    df['Site Passes Threshold [0.05]']=np.where(df['Site Q-Value'] <= 0.05, 1, 0)
    df['Site Passes Threshold [0.01]'] = np.where(df['Site Q-Value'] <= 0.01, 1, 0)
    df['Decoy Peptide']=np.where(df['Proteins'].str.contains("DECOY")==True, 1, 0)
    df['Decoy Modification Site']=df.apply(lambda x: r(x.Pep_pos.split("-")[0], x.Pep_pos.split("-")[1]), axis=1)
    df['Peptidoform Site ID']=df.index

    ######
    meta_file=wd+"PRIDE/SDRFs/"+i+".sdrf.tsv"
    meta=pd.read_csv(meta_file,sep="\t")
    meta['comment[data file]']=meta['comment[data file]'].str.split(".").str[0]
    meta['comment[data file]']=meta['comment[data file]'].str.replace(" ","_")
    meta.set_index('comment[data file]', inplace=True)
    dict=meta.to_dict()

    PMID=[]
    Sample=[]
    dataset_ID=[]
    organism=[]
    organism_part=[]
    cell_line=[]
    disease=[]
    for x in range(len(df)):
        source_all=df.loc[x,"All_Source"]
        for source in source_all.split(";"):
            source=source.replace("_raw","")
            PMID_temp=[]
            Sample_temp=[]
            dataset_ID_temp=[]
            organism_temp=[]
            organism_part_temp=[]
            cell_line_temp=[]
            disease_temp=[]
            Sample_temp.append(dict["source name"][source])
            PXD=dict["comment[proteomexchange accession number]"][source]
            PMID_temp.append(dict_all[PXD][0])
            dataset_ID_temp.append(PXD)
            organism_temp.append(dict["characteristics[organism]"][source])
            organism_part_temp.append(dict["characteristics[organism part]"][source])
            cell_line_temp.append(dict["characteristics[cell type]"][source])
            disease_temp.append(dict["characteristics[disease]"][source])
        PMID.append(";".join(map(str,list(set(PMID_temp)))))
        Sample.append(";".join(map(str,list(set(Sample_temp)))))
        dataset_ID.append(";".join(map(str,list(set(dataset_ID_temp)))))
        organism.append(";".join(map(str,list(set(organism_temp)))))
        organism_part.append(";".join(map(str,list(set(organism_part_temp)))))
        cell_line.append(";".join(map(str,list(set(cell_line_temp)))))
        disease.append(";".join(map(str,list(set(disease_temp)))))

    df['PubMedIDs']=PMID
    df['Sample ID']=Sample
    df['Source Dataset Identifier']=dataset_ID
    df['Reanalysis Dataset Identifier']="NA"
    df['Organism']=organism
    df['Organism Part']=organism_part
    df['Cell Line']=cell_line
    df['Disease Information']=disease
    #####

    df = df[['Peptidoform Site ID', 'Proteins', 'Unmodified Sequence', 'Peptidoform','Modification','Peptide Modification Position', 'PSM Probability', 'PTM Probability',
    'Final Probability','Site Q-Value','Site Passes Threshold [0.05]','Site Passes Threshold [0.01]', 'Decoy Peptide','Decoy Modification Site',
    'PSM Count Passing Threshold [0.05]','PSM Count Passing Threshold [0.01]', 'Source Dataset Identifier', 'Reanalysis Dataset Identifier','PubMedIDs',
    'Sample ID', 'Organism', 'Organism Part', 'Cell Line', 'Disease Information', 'Universal Spectrum Identifier','opt_PSM count 0.01<P<=0.05','opt_PSM count 0.05<P<=0.19',
    'opt_PSM count 0.19<P<=0.81','opt_PSM count 0.81<P<=0.95','opt_PSM count 0.95<P<0.99','opt_PSM count P>=0.99']]
    df = df.replace("not applicable","NA")
    df.to_csv(wd+"All_site_formats/"+dataset_ID[0]+"_Site_Peptidoform_centric.tsv", sep="\t", index=False)


print("Peptidoform formats done, creating GSB counts")
flr_filter=0.05
### GSB counts
for loc in folder_list:
    loc_full=wd+"/"+loc+"/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
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
    df_final.loc[i,'Protein']=i.rsplit("-",1)[0]
    df_final.loc[i,'Protein_pos']=i.rsplit("-",1)[1].split("_")[0]
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

cols = list(df_final.columns.values)
cols.pop(cols.index('Protein'))
cols.pop(cols.index('Protein_pos'))
cols.pop(cols.index('PTM_residue'))
cols.pop(cols.index('PTM_FLR_category'))
df_final = df_final[['Protein','Protein_pos','PTM_residue','PTM_FLR_category']+cols]

df_final.to_csv(wd+"All_site_formats/GSB_"+str(flr_filter)+"_protein_pos.csv", index=False)
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
plt.savefig(wd+"All_site_formats/GSB_"+str(flr_filter)+"_protein_pos_categories.png",dpi=300)

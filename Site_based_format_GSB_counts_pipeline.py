import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
#import requests
#import json
import sys

# requires csv file (meta.csv) with meta data, sdrf location (eg. "PRIDE/SDRFs/"), gold count and silver count - if meta or sdrf files not available, use "NA"

# usage: python Site_based_format_GSB_counts_pipeline.py file_list.txt meta.csv PRIDE/SDRFs/ 5 2
# OR
# python Site_based_format_GSB_counts_pipeline.py file_list.txt NA NA 5 2

#read in csv from FLR pipeline
#DECOY and CONTAM removed, filtered for contains Phospho, exploded for site based, binomial adjustment and collapsed by peptidoform

folder_list_file = open(sys.argv[1],"r")
folder_list_file = folder_list_file.read()
folder_list = folder_list_file.replace('\n', '.').split(".")

dataset_list=(list(set([a.split("/",1)[0] for a in folder_list])))

meta_all = sys.argv[2]
#meta_all=wd+"/meta.csv"
if meta_all!="NA":
    with open(meta_all,'r') as f:
        reader = csv.reader(f)
        dict_all = {rows[0]:rows[1:] for rows in reader}
    f.close()

gold_count = int(sys.argv[4])
silver_count = int(sys.argv[5])

if not os.path.exists("All_site_formats"):
    os.mkdir("All_site_formats")

r = lambda x, y : 0 if x[int(float(y))-1]!="A" else 1

## Peptide site format
for i in folder_list:
    print(i)
    file=i+"/FDR_0.01/binomial.csv"
    df_temp = pd.read_csv(file)
    search_mod="Phospho"
    #Filter for decoy FLR? - not needed due to "Site Passes Threshold" column
    df = df_temp[['All_Proteins','All_PTM_protein_positions','Peptide','USI','pA_q_value','PTM Score','All_PTM_positions','Score','PTM positions','PTM_info','Binomial_final_score','Peptide_mod','All_PTMs',
                  'All_PTM_scores','All_USI','Protein position','Source']].copy()
    df = df.rename(columns={'All_Proteins':'Proteins','All_PTM_protein_positions':'Protein Modification Positions','pA_q_value':'Site Q-Value','Score':'PSM Probability','PTM Score':'PTM Probability',
                            'PTM positions':'Peptide Modification Position','Binomial_final_score':'Final Probability','All_PTMs':'Modifications','All_PTM_scores':'Modification probabilities',
                            'All_PTM_positions':'Modification positions','Peptide_mod':'Peptidoform','Peptide':'Unmodified Sequence','USI':'Universal Spectrum Identifier'}) #'All_USI':'Supporting PSM count',
    df = df.reset_index(drop=True)
    df['Proteins'] = df['Proteins'].str.replace(":",";")
    df['Modification']=search_mod
    dataset_ID=[]

    PMID = []
    Sample = []
    dataset_ID = []
    organism = []
    organism_part = []
    cell_line = []
    disease = []
    for x in range(len(df)):
        protein_PTM = ""
        for z, y in zip(df.loc[x, 'Proteins'].split(";"), df.loc[x, 'Protein Modification Positions'].split(":")):
            for a, b in zip(df.loc[x, 'Modification positions'].split(";"), y.split(";")):
                if int(a) == int(df.loc[x, 'Peptide Modification Position']):
                    protein_PTM += b + ";"
        df.loc[x, 'Protein Modification Positions'] = protein_PTM[:-1]

        source = str(df.loc[x, "Source"])
        source = source.replace("_raw", "")
        source = source.replace("(","")
        source = source.replace(")","")

        if sys.argv[3]!="NA":
            meta_file = sys.argv[3]+"/"+i.split("/")[0]+".sdrf.tsv"
            meta=pd.read_csv(meta_file,sep="\t")
            meta['comment[data file]']=meta['comment[data file]'].str.split(".").str[0]
            meta['comment[data file]']=meta['comment[data file]'].str.replace(" ","_")
            meta.set_index('comment[data file]', inplace=True)
            dict=meta.to_dict()
            Sample.append(dict["source name"][source])
            organism.append(dict["characteristics[organism]"][source])
            cell_line.append(dict["characteristics[cell type]"][source])
        else:
            Sample.append("NA")
            organism.append("NA")
            cell_line.append("NA")

        try:
            PXD = dict["comment[proteomexchange accession number]"][source]
        except:
            PXD="NA"
        try:
            organism_part.append(dict["characteristics[organism part]"][source])
        except:
            organism_part.append("NA")
        try:
            disease.append(dict["characteristics[disease]"][source])
        except:
            disease.append("NA")

        if meta_all!="NA":
            PMID.append(dict_all[PXD][0])
        else:
            PMID.append("NA")
        dataset_ID.append(PXD)

    df['PubMedIDs'] = PMID
    df['Sample ID'] = Sample
    df['Source Dataset Identifier'] = dataset_ID
    df['Reanalysis Dataset Identifier'] = "NA"
    df['Organism'] = organism
    df['Organism Part'] = organism_part
    df['Cell Line'] = cell_line
    df['Disease Information'] = disease

    #####

    df['Pep_pos'] = df['Unmodified Sequence'] + "-" + df['Peptide Modification Position'].astype(str)
    df['Site Passes Threshold [0.05]'] = np.where(df['Site Q-Value'] <= 0.05, 1, 0)
    df['Site Passes Threshold [0.01]'] = np.where(df['Site Q-Value'] <= 0.01, 1, 0)
    df['Decoy Peptide']=np.where(df['Proteins'].str.contains("DECOY")==True, 1, 0)
    df['Decoy Modification Site']=df.apply(lambda x: r(x.Pep_pos.split("-")[0], x.Pep_pos.split("-")[1]), axis=1)
    df['PSM Site ID']=df.index
    df['PSM Count Passing Threshold [0.05]']="1"

    #peptidoform without reagent labels?
    df['Peptidoform']=df['Peptidoform'].str.replace(r'\[TMT(.*?)\]','',regex="True")
    df['Peptidoform']=df['Peptidoform'].str.replace(r'\[iTRAQ.*?\]','',regex="True")

    df = df[['PSM Site ID', 'Proteins', 'Unmodified Sequence', 'Peptidoform','Modification','Peptide Modification Position', 'Protein Modification Positions', 'PSM Probability','PTM Probability',
             'Final Probability','Site Q-Value','Site Passes Threshold [0.05]','Site Passes Threshold [0.01]', 'Decoy Peptide','Decoy Modification Site',
             'PSM Count Passing Threshold [0.05]','Source Dataset Identifier', 'Reanalysis Dataset Identifier', 'PubMedIDs','Sample ID', 'Organism', 'Organism Part',
             'Cell Line','Disease Information', 'Universal Spectrum Identifier']]
    df = df.replace("not applicable","NA")

    df.to_csv("All_site_formats/"+i.replace("/","_")+"_Site_PSM_centric.tsv", sep="\t", index=False)

    # #Test for valid USI
    # #Won't work on cluster!
    # df=df.head(100)
    # #Check USI valid
    # url = 'https://proteomecentral.proteomexchange.org/api/proxi/v0.1/usi_validator'
    # x = requests.post(url, json=df['Universal Spectrum Identifier'].to_list())
    # res=json.loads(x.text)
    #
    # for i in range(len(df)):
    #     #df.loc[i,'Valid_USI']=res['validation_results'][df.loc[i,'USI']]['is_valid']
    #     if res['validation_results'][df.loc[i,'Universal Spectrum Identifier']]['is_valid']==False:
    #         print(df.loc[i,'Universal Spectrum Identifier'])

# Merged files
for dataset in dataset_list:
    folder_list_temp=[x for x in folder_list if dataset in x]
    for loc in folder_list_temp:
        loc_full="All_site_formats//"+loc.replace("/","_")+"_Site_PSM_centric.tsv"
        df=pd.read_csv(loc_full, sep="\t", dtype={'Protein Modification Positions': str})
        if loc==folder_list_temp[0]:
            df_temp=df
        else:
            df_temp=pd.concat([df_temp,df])
    df_temp=df_temp.reset_index(drop=True)
    df_temp['PSM Site ID']=df_temp.index

    df_temp.to_csv("All_site_formats/"+dataset+"_merged_Site_PSM_centric.tsv", sep="\t", index=False)

    if dataset==dataset_list[0]:
        df_all=df_temp
    else:
        df_all=pd.concat([df_temp,df_all])
df_all.to_csv("All_site_formats//all_datasets_merged_Site_PSM_centric.tsv", sep="\t", index=False)


### Peptidoform site format
print("PSM format complete, creating Peptidoform format")
for i in folder_list:
    print(i)
    file=i+"/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
    df_temp = pd.read_csv(file)
    search_mod="Phospho"
    #Filter for decoy FLR? - not needed due to "Site Passes Threshold" column
    df = df_temp[['All_Proteins','All_PTM_protein_positions','Peptide','USI','pA_q_value_BA','PTM Score','All_PTM_positions','Score','PTM positions','PTM_info','Binomial_final_score','Peptide_mod','All_PTMs',
                  'All_PTM_scores','All_USI','Protein position','0.05FLR_threshold_count','0.01FLR_threshold_count','0.01<P<=0.05_count','0.05<P<=0.19_count','0.19<P<=0.81_count', '0.81<P<=0.95_count',
                  '0.95<P<0.99_count', 'P>=0.99_count','All_Source']].copy()
    df = df.rename(columns={'All_Proteins':'Proteins','All_PTM_protein_positions':'Protein Modification Positions','Peptide_mod':'Peptidoform','Peptide':'Unmodified Sequence','pA_q_value_BA':'Site Q-Value','Score':'PSM Probability','PTM Score':'PTM Probability',
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

    if sys.argv[3]!="NA":
        meta_file = sys.argv[3]+"/"+i.split("/")[0]+".sdrf.tsv"
        meta=pd.read_csv(meta_file,sep="\t")
        meta['comment[data file]']=meta['comment[data file]'].str.split(".").str[0]
        meta['comment[data file]']=meta['comment[data file]'].str.replace(" ","_")
        meta.set_index('comment[data file]', inplace=True)
        dict=meta.to_dict()

    PMID = []
    Sample = []
    dataset_ID = []
    organism = []
    organism_part = []
    cell_line = []
    disease = []

    for x in range(len(df)):
        protein_PTM = ""
        for z, y in zip(df.loc[x, 'Proteins'].split(";"), df.loc[x, 'Protein Modification Positions'].split(":")):
            for a, b in zip(df.loc[x, 'Modification positions'].split(";"), y.split(";")):
                if int(a) == df.loc[x, 'Peptide Modification Position']:
                    protein_PTM += b + ";"
        df.loc[x, 'Protein Modification Positions'] = protein_PTM[:-1]

        source_all = df.loc[x, "All_Source"]
        for source in source_all.split(";"):
            source = source.replace("_raw", "")
            source = source.replace("(","")
            source = source.replace(")","")
            PMID_temp = []
            Sample_temp = []
            dataset_ID_temp = []
            organism_temp = []
            organism_part_temp = []
            cell_line_temp = []
            disease_temp = []

            if sys.argv[3]!="NA":
                Sample_temp.append(dict["source name"][source])
                PXD = dict["comment[proteomexchange accession number]"][source]
                dataset_ID_temp.append(PXD)
                organism_temp.append(dict["characteristics[organism]"][source])
                try:
                    organism_part_temp.append(dict["characteristics[organism part]"][source])
                except:
                    organism_part_temp.append("NA")
                cell_line_temp.append(dict["characteristics[cell type]"][source])
                try:
                    disease_temp.append(dict["characteristics[disease]"][source])
                except:
                    disease_temp.append("NA")
            else:
                Sample_temp.append("NA")
                dataset_ID_temp.append("NA")
                organism_temp.append("NA")
                organism_part_temp.append("NA")
                cell_line_temp.append("NA")
                disease_temp.append("NA")

            if meta_all!="NA":
                PMID_temp.append(dict_all[PXD][0])
            else:
                PMID_temp.append("NA")

        PMID.append(";".join(map(str, list(set(PMID_temp)))))
        Sample.append(";".join(map(str, list(set(Sample_temp)))))
        dataset_ID.append(";".join(map(str, list(set(dataset_ID_temp)))))
        organism.append(";".join(map(str, list(set(organism_temp)))))
        organism_part.append(";".join(map(str, list(set(organism_part_temp)))))
        cell_line.append(";".join(map(str, list(set(cell_line_temp)))))
        disease.append(";".join(map(str, list(set(disease_temp)))))

    df['PubMedIDs'] = PMID
    df['Sample ID'] = Sample
    df['Source Dataset Identifier'] = dataset_ID
    df['Reanalysis Dataset Identifier'] = "NA"
    df['Organism'] = organism
    df['Organism Part'] = organism_part
    df['Cell Line'] = cell_line
    df['Disease Information'] = disease
    #####

    df = df[['Peptidoform Site ID', 'Proteins', 'Unmodified Sequence', 'Peptidoform','Modification','Peptide Modification Position','Protein Modification Positions', 'PSM Probability', 'PTM Probability',
    'Final Probability','Site Q-Value','Site Passes Threshold [0.05]','Site Passes Threshold [0.01]', 'Decoy Peptide','Decoy Modification Site',
    'PSM Count Passing Threshold [0.05]','PSM Count Passing Threshold [0.01]', 'Source Dataset Identifier', 'Reanalysis Dataset Identifier','PubMedIDs',
    'Sample ID', 'Organism', 'Organism Part', 'Cell Line', 'Disease Information', 'Universal Spectrum Identifier','opt_PSM count 0.01<P<=0.05','opt_PSM count 0.05<P<=0.19',
    'opt_PSM count 0.19<P<=0.81','opt_PSM count 0.81<P<=0.95','opt_PSM count 0.95<P<0.99','opt_PSM count P>=0.99']]
    df = df.replace("not applicable","NA")
    df.to_csv("All_site_formats/"+i.replace("/","_")+"_Site_Peptidoform_centric.tsv", sep="\t", index=False)

#Merged files
for dataset in dataset_list:
    folder_list_temp=[x for x in folder_list if dataset in x]
    for loc in folder_list_temp:
        loc_full="All_site_formats/"+loc.replace("/","_")+"_Site_Peptidoform_centric.tsv"
        df=pd.read_csv(loc_full, sep="\t", dtype={'Protein Modification Positions': str})

        if loc==folder_list_temp[0]:
            df_temp=df
        else:
            df_temp=pd.concat([df_temp,df])

    df_temp['Peptidoform_pos']=df_temp['Peptidoform']+"-"+df_temp['Peptide Modification Position'].astype(str)

    df_temp=df_temp.sort_values(['Peptidoform_pos','Final Probability','Site Q-Value'],ascending=[True,True,False])
    df_temp=df_temp.drop_duplicates(subset=('Peptidoform_pos'),keep="last",inplace=False)
    df_temp=df_temp.reset_index()
    df_temp['Peptidoform Site ID']=df_temp.index

    df_merge=df_temp.drop(['index','Peptidoform_pos'], axis=1)
    df_merge.to_csv("All_site_formats/"+dataset+"_merged_Site_Peptidoform_centric.tsv", sep="\t", index=False)

    if dataset==dataset_list[0]:
        df_all=df_merge
    else:
        df_all=pd.concat([df_merge,df_all])
df_all.to_csv("All_site_formats//all_datasets_merged_Site_Peptidoform_centric.tsv", sep="\t", index=False)


print("Peptidoform formats done, creating GSB counts")
flr_filter=0.05
### GSB counts
#filter for sites seen multiple times across dataset studies
for dataset in dataset_list:
    folder_list_temp=[x for x in folder_list if dataset in x]
    for loc in folder_list_temp:
        loc_full=loc+"/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
        df=pd.read_csv(loc_full)

        df=df.loc[df['pA_q_value_BA']<=flr_filter]

        if loc==folder_list_temp[0]:
            df_temp=df
        else:
            df_temp=pd.concat([df_temp,df])


    df_temp=df_temp.sort_values(['Peptide_mod_pos','Binomial_final_score','pA_q_value_BA'],ascending=[True,True,False])
    df_temp=df_temp.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)

    df_temp=df_temp.sort_values(['Protein-pos','Binomial_final_score','pA_q_value_BA'],ascending=[True,True,False])

    if len(df_temp)>=1:
        df_temp['PTM_residue']=df_temp.apply(lambda x: x['Peptide'][x['PTM positions']-1],axis=1)
    else:
        df_temp['PTM_residue']=""
    df_temp['Protein_pos_res']=df_temp['Protein-pos']+"_"+df_temp['PTM_residue']
    df_temp=df_temp.drop_duplicates(subset=('Protein_pos_res'),keep="last",inplace=False)

    df_temp=df_temp[['Peptide_mod_pos', 'pA_q_value_BA','Binomial_final_score', 'Protein_pos_res']]
    df_temp.rename(columns = {'pA_q_value_BA':dataset+"_FLR",'Binomial_final_score':dataset+'_BinomialScore','Peptide_mod_pos':dataset+"_peptide_mod_pos"}, inplace = True)
    df_temp=df_temp.set_index('Protein_pos_res')

    if dataset==dataset_list[0]:
        df_final=df_temp
    else:
        df_final=pd.concat([df_final,df_temp],axis=1)


for i in df_final.index.values.tolist():
    df_final.loc[i,'Protein']=i.rsplit("-",1)[0]
    df_final.loc[i,'Protein_pos']=i.rsplit("-",1)[1].split("_")[0]
    df_final.loc[i,'PTM_residue']=i.rsplit("_",1)[-1]
    flr1_count=0
    for dataset in dataset_list:
        FLR_count_all=0
        FLR_cols = [x for x in dataset_list if dataset in x]
        for FLR_col in FLR_cols:
            if df_final.loc[i,FLR_col+"_FLR"]!="N/A":
                if float(df_final.loc[i,FLR_col+"_FLR"])<=0.01:
                    FLR_count_all+=1
        if FLR_count_all!=0:
            flr1_count+=1
    if flr1_count>=gold_count:
        df_final.loc[i,'PTM_FLR_category']="Gold"
    elif flr1_count>=silver_count:
        df_final.loc[i,'PTM_FLR_category']="Silver"
    else:
        df_final.loc[i,'PTM_FLR_category']="Bronze"

cols = list(df_final.columns.values)
cols.pop(cols.index('Protein'))
cols.pop(cols.index('Protein_pos'))
cols.pop(cols.index('PTM_residue'))
cols.pop(cols.index('PTM_FLR_category'))
df_final = df_final[['Protein','Protein_pos','PTM_residue','PTM_FLR_category']+cols]

df_final.to_csv("All_site_formats/GSB_"+str(flr_filter)+"_protein_pos.csv", index=False)
counts=pd.crosstab(df_final.PTM_residue,df_final.PTM_FLR_category).replace(0,np.nan).stack().reset_index().rename(columns={0:'Count of Sites'})
print(counts)
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
plt.savefig("All_site_formats/GSB_"+str(flr_filter)+"_protein_pos_categories.png",dpi=300)

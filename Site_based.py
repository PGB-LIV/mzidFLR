import pandas as pd
import numpy as np
# import requests
# import json

#read in csv from FLR pipeline
#DECOY and CONTAM removed, filtered for contains Phospho, exploded for site based, binomial adjustment and collapsed by peptidoform

## Usage - Site_based.py [optional: PMID] [optional: SDRF file]
## If no PMID or USI given, "NA" used
try:
    PMID = args[1]
    if ".sdrf.tsv" in SDRF:
        SDRF = args[1]
        try:
            PMID = args[2]
        except:
            PMID = "NA"
except:
    PMID = "NA"

try:
    SDRF = args[2]
    if ".sdrf.tsv" in SDRF:
        SDRF = args[2]
    else:
        SDRF = "NA"
        try:
            PMID = args[2]
        except:
            PMID = "NA"
except:
    SDRF = "NA"

r = lambda x, y : 0 if x[int(y)-1]!="A" else 1

### Peptide site format

file="FDR_0.01/binomial.csv"
df_temp = pd.read_csv(file)
search_mod="Phospho"
df = df_temp[['All_Proteins','All_PTM_protein_positions','Peptide','USI','pA_q_value','PTM Score','Score','PTM positions','PTM_info','Binomial_final_score','Peptide_mod','All_PTMs','All_PTM_positions',
              'All_PTM_scores','All_USI','Protein position','Source']].copy()
df = df.rename(columns={'All_Proteins':'Proteins','All_PTM_protein_positions':'Protein Modification Positions','pA_q_value':'Site Q-Value','Score':'PSM Probability','PTM Score':'PTM Probability',
                        'PTM positions':'Peptide Modification Position','Binomial_final_score':'Final Probability','All_PTMs':'Modifications','All_PTM_scores':'Modification probabilities',
                        'All_PTM_positions':'Modification positions','Peptide_mod':'Peptidoform','Peptide':'Unmodified Sequence','USI':'Universal Spectrum Identifier'}) #'All_USI':'Supporting PSM count',
df = df.reset_index(drop=True)
df['Proteins'] = df['Proteins'].str.replace(":",";")
df['Modification']=search_mod

######
if SDRF != "NA":
    meta=pd.read_csv(SDRF,sep="\t")
    meta['comment[data file]']=meta['comment[data file]'].str.split(".").str[0]
    meta['comment[data file]']=meta['comment[data file]'].str.replace(" ","_")
    meta.set_index('comment[data file]', inplace=True)
    dict=meta.to_dict()


PMID_list=[]
Sample=[]
dataset_ID=[]
organism=[]
organism_part=[]
cell_line=[]
disease=[]
for x in range(len(df)):
    protein_PTM=""
    for z,y in zip(df.loc[x,'Proteins'].split(";"),df.loc[x,'Protein Modification Positions'].split(":")):
        for a,b in zip(df.loc[x,'Modification positions'].split(";"),y.split(";")):
            if int(a)==int(df.loc[x,'Peptide Modification Position']):
                protein_PTM+=b+";"
    df.loc[x,'Protein Modification Positions']=protein_PTM[:-1]

    source=df.loc[x,"Source"]
    source=source.replace("_raw","")
    try:
        PXD=dict["comment[proteomexchange accession number]"][source]
    except:
        PXD = df.loc[x,'Universal Spectrum Identifier'].split(":")[1]
    dataset_ID.append(PXD)
    PMID_list.append(PMID)
    try:
        Sample.append(dict["source name"][source])
    except:
        Sample.append("NA")
    try:
        organism.append(dict["characteristics[organism]"][source])
    except:
        organism.append("NA")
    try:
        organism_part.append(dict["characteristics[organism part]"][source])
    except:
        organism_part.append("NA")
    try:
        cell_line.append(dict["characteristics[cell type]"][source])#
    except:
        cell_line.append("NA")
    try:
        disease.append(dict["characteristics[disease]"][source])
    except:
        disease.append("NA")
df['PubMedIDs']=PMID_list
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


df = df[['PSM Site ID', 'Proteins', 'Unmodified Sequence', 'Peptidoform','Modification','Peptide Modification Position', 'Protein Modification Positions','PSM Probability','PTM Probability',
         'Final Probability','Site Q-Value','Site Passes Threshold [0.05]','Site Passes Threshold [0.01]', 'Decoy Peptide','Decoy Modification Site',
         'PSM Count Passing Threshold 0.05','Source Dataset Identifier', 'Reanalysis Dataset Identifier', 'PubMedIDs','Sample ID', 'Organism', 'Organism Part',
         'Cell Line','Disease Information', 'Universal Spectrum Identifier']]
df = df.replace("not applicable","NA")

df.to_csv("FDR_0.01/Site_PSM_centric.tsv", sep="\t", index=False)

# #Check USI valid
# df=df.head(100)
# url = 'https://proteomecentral.proteomexchange.org/api/proxi/v0.1/usi_validator'
# x = requests.post(url, json=df['Universal Spectrum Identifier'].to_list())
# res=json.loads(x.text)
#
# for i in range(len(df)):
#     #df.loc[i,'Valid_USI']=res['validation_results'][df.loc[i,'USI']]['is_valid']
#     if res['validation_results'][df.loc[i,'Universal Spectrum Identifier']]['is_valid']==False:
#         print(df.loc[i,'Universal Spectrum Identifier'])



### Peptidoform site format

print("PSM format complete, creating Peptidoform format")

file="FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
df_temp = pd.read_csv(file)

search_mod="Phospho"
#Filter for decoy FLR? - not needed due to "Site Passes Threshold" column
df = df_temp[['All_Proteins','All_PTM_protein_positions','Peptide','USI','pA_q_value_BA','PTM Score','Score','PTM positions','PTM_info','Binomial_final_score','Peptide_mod','All_PTMs','All_PTM_positions',
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

df['Pep_pos'] = df['Unmodified Sequence'] + "-" + df['Peptide Modification Position'].astype(str)
df['Site Passes Threshold [0.05]']=np.where(df['Site Q-Value'] <= 0.05, 1, 0)
df['Site Passes Threshold [0.01]'] = np.where(df['Site Q-Value'] <= 0.01, 1, 0)
df['Decoy Peptide']=np.where(df['Proteins'].str.contains("DECOY")==True, 1, 0)
df['Decoy Modification Site']=df.apply(lambda x: r(x.Pep_pos.split("-")[0], x.Pep_pos.split("-")[1]), axis=1)
df['Peptidoform Site ID']=df.index

PMID=[]
Sample=[]
dataset_ID=[]
organism=[]
organism_part=[]
cell_line=[]
disease=[]


for x in range(len(df)):
    protein_PTM=""
    for z,y in zip(df.loc[x,'Proteins'].split(";"),df.loc[x,'Protein Modification Positions'].split(":")):
        for a,b in zip(df.loc[x,'Modification positions'].split(";"),y.split(";")):
            if int(a)==df.loc[x,'Peptide Modification Position']:
                protein_PTM+=b+";"
    df.loc[x,'Protein Modification Positions']=protein_PTM[:-1]

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
        try:
            Sample_temp.append(dict["source name"][source])
        except:
            Sample_temp.append("NA")
        try:
            PXD=dict["comment[proteomexchange accession number]"][source]
        except:
            PXD=df.loc[x,'Universal Spectrum Identifier'].split(":")[1]
        try:
            PMID_temp.append(dict_all[PXD][0])
        except:
            PMID_temp.append("NA")
        dataset_ID_temp.append(PXD)
        try:
            organism_temp.append(dict["characteristics[organism]"][source])
        except:
            organism_temp.append("NA")
        try:
            organism_part_temp.append(dict["characteristics[organism part]"][source])
        except:
            organism_part_temp.append("NA")
        try:
            cell_line_temp.append(dict["characteristics[cell type]"][source])
        except:
            cell_line_temp.append("NA")
        try:
            disease_temp.append(dict["characteristics[disease]"][source])
        except:
            disease_temp.append("NA")
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

df = df[['Peptidoform Site ID', 'Proteins', 'Unmodified Sequence', 'Peptidoform','Modification','Peptide Modification Position', 'Protein Modification Positions', 'PSM Probability', 'PTM Probability',
         'Final Probability','Site Q-Value','Site Passes Threshold [0.05]','Site Passes Threshold [0.01]', 'Decoy Peptide','Decoy Modification Site',
         'PSM Count Passing Threshold [0.05]','PSM Count Passing Threshold [0.01]', 'Source Dataset Identifier', 'Reanalysis Dataset Identifier','PubMedIDs',
         'Sample ID', 'Organism', 'Organism Part', 'Cell Line', 'Disease Information', 'Universal Spectrum Identifier','opt_PSM count 0.01<P<=0.05','opt_PSM count 0.05<P<=0.19',
         'opt_PSM count 0.19<P<=0.81','opt_PSM count 0.81<P<=0.95','opt_PSM count 0.95<P<0.99','opt_PSM count P>=0.99']]

df.to_csv("FDR_0.01/Site_Peptidoform_centric.tsv", sep="\t", index=False)

print("Peptidoform format complete")



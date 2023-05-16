import pandas as pd
#from scipy.stats import binom
import numpy as np
import csv
import os

#read in csv from FLR pipeline
#DECOY and CONTAM removed, filtered for contains Phospho, exploded for site based, binomial adjustment and collapsed by peptidoform

meta="C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/rice_ptm_metadata.csv"
with open(meta,'r') as f:
    reader = csv.reader(f)
    dict = {rows[0]:rows[1:] for rows in reader}
f.close()

wd="C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/"
folder_list=["PXD000923","PXD002222","PXD002756","PXD004705","PXD004939","PXD005241","PXD012764","PXD019291"]#,"PXD001168"]

for i in folder_list:
    file=wd+i+"/FDR_updated_0.01/binomial_peptidoform_collapsed_FLR.csv"
    df_temp = pd.read_csv(file)
    FLR_threshold=0.05
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
    df['Site Passes Threshold']=np.where(df['Binomial adjusted pA q value'] <= FLR_threshold, 1, 0)
    df['Dataset ID']=dataset_ID
    df['Site ID']=df.index
    df['Organism']=organism
    df['PubMedID']=PMID

    df=df.replace("pryo","pyro", regex=True)
    df=df.replace("Q\[Pyro_glu\]","Q[Gln->pyro-Glu]", regex=True)
    df=df.replace("E\[Pyro_glu\]","E[Glu->pyro-Glu]", regex=True)
    df=df.replace("C\[Pyro_glu\]","C[Ammonia-loss]", regex=True)
    df=df.replace("Carbamidomethylation","Carbamidomethyl", regex=True)
    df=df.replace("346.212775","iTRAQ8plex", regex=True)
    df=df.replace("346.212800","iTRAQ8plex", regex=True)
    df=df.replace("229.162975","TMT6plex", regex=True)
    df=df.replace("229.162932","TMT6plex", regex=True)
    df = df.replace("346.21", "iTRAQ8plex", regex=True)
    df = df.replace("229.16", "TMT6plex", regex=True)

    df = df[['Site ID', 'Proteins', 'Unmodified Sequence', 'Modified Peptide','Modification','Site position', 'Site probability', 'PSM Score',
            'Final site probability','Binomial adjusted pA q value','Site Passes Threshold', 'Supporting PSM count','Supporting PSM Count (Thresholded 0.05 Global FLR)', 'Dataset ID', 'PRIDE dataset URL', 'PubMedID',
             'Sample ID', 'Organism', 'USI','PRIDE dataset URL (including the USI)','PSM count 0.01<P<=0.05','PSM count 0.05<P<=0.19','PSM count 0.19<P<=0.81','PSM count 0.81<P<=0.95','PSM count 0.95<P<0.99','PSM count P>=0.99']]
    if not os.path.exists(wd+"Site_format_files"):
        os.mkdir(wd+"Site_format_files")
    df.to_csv(wd+"Site_format_files/"+dataset_ID[0]+"_Site_peptidoform_centric_binomial_FLR_no_pA_bins_"+str(FLR_threshold)+".tsv", sep="\t", index=False)

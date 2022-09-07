import pandas as pd
#from scipy.stats import binom
import numpy as np
import csv

#read in csv from FLR pipeline
#DECOY and CONTAM removed, filtered for contains Phospho, exploded for site based, binomial adjustment and collapsed by position
PXD000923 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD000923/FDR_0.01_PTM_score_0/binomial_collapsed_FLR.csv"
PXD002222 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD002222/FDR_0.01_PTM_score_0/binomial_collapsed_FLR.csv"
PXD002756 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD002756/FDR_0.01_PTM_score_0/binomial_collapsed_FLR.csv"
PXD004705 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD004705/FDR_0.01_PTM_score_0/binomial_collapsed_FLR.csv"
PXD004939 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD004939/FDR_0.01_PTM_score_0/binomial_collapsed_FLR.csv"
PXD005241 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD005241/FDR_0.01_PTM_score_0/binomial_collapsed_FLR.csv"
PXD012764 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD012764/FDR_0.01_PTM_score_0/binomial_collapsed_FLR.csv"
PXD019291 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD019291/FDR_0.01_PTM_score_0/binomial_collapsed_FLR.csv"

meta="C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/rice_ptm_metadata.csv"
with open(meta,'r') as f:
    reader = csv.reader(f)
    dict = {rows[0]:rows[1:] for rows in reader}
f.close()

for i in [PXD000923,PXD002222,PXD002756,PXD004705,PXD004939,PXD005241,PXD012764,PXD019291]:
    df_temp = pd.read_csv(i)
    FLR_threshold=0.05
    search_mod="Phosphorylation"
    #Filter for decoy FLR? - not needed due to "Site Passes Threshold" column
    df = df_temp[['All_Proteins','Peptide','USI','Binomial_final_prob_q_value','PTM Score','All_PTM_positions','Score','DecoyP','PTM positions','PTM_info','Binomial_final_score','Peptide_mod','All_PTMs','All_PTM_positions',
                  'All_PTM_scores','All_USI','Protein position','Peptide_start_Protein','Peptide_end_Protein','pA_prob_chance_binomial',
                  '0.05FLR_threshold_count','0.01<P<=0.05_count','0.05<P<=0.19_count','0.19<P<=0.81_count', '0.81<P<=0.95_count', '0.95<P<0.99_count', 'P>=0.99_count']].copy()
    df = df.rename(columns={'All_Proteins':'Proteins','Binomial_final_prob_q_value':'Binomial final q value','Score':'PSM Score','PTM Score':'Site probability','All_USI':'Supporting PSM count',
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

    #Repeat binomial adjustment for decoy replacements
    # df = df.sort_values(by=(['Proteins','Protein position']), ascending=[False,False])
    # df = df.reset_index(drop=True)
    # for i in range(len(df)):
    #     protein=df.loc[i,'Proteins'].split(";")[0]
    #     pos=df.loc[i,'Protein position']
    #     df.loc[i,'Protein_temp']=protein
    #     df.loc[i,'Protein-pos']=protein+"-"+str(pos)
    #     df.loc[i,'Modification']=search_mod
    # df2= df.groupby(['Protein-pos'])['Supporting PSM count'].sum()
    # dict=df2.to_dict()
    #
    # for i in range(len(df)):
    #     protein=df.loc[i,'Proteins'].split(";")[0]
    #     pos=df.loc[i,'Protein position']
    #     df_temp=df.loc[(df['Protein_temp']==protein) & (df['Peptide_start_Protein']<=pos) & (df['Peptide_end_Protein']>=pos)]
    #     df.loc[i,'Protein_occ']=df_temp['Supporting PSM count'].sum()
    #     df.loc[i,'Protein_loc']=dict[df.loc[i,'Protein-pos']]
    # #df.to_csv("temp.csv", sep=",", index=False)
    # df['Binomial_p'] = binom.pmf(df['Protein_loc'],df['Protein_occ'],df['pA_prob_chance_binomial'])
    # df['Binomial_i'] = df['Site probability'] * (1-df['Binomial_p'])
    # df['Final adjusted site probability'] = df['Binomial_i'] * df['PSM Score']

    #df.to_csv("Site_peptidoform_centric_final_updated_binomial_no_pA_temp.csv", sep="\t", index=False)

    #df.to_csv("Site_peptidoform_final_temp.csv", sep="\t", index=False)
    df['Pep_pos'] = df['Unmodified Sequence'] + "-" + df['Site position'].astype(str)
    df = df.sort_values(by=(['Pep_pos','Modified Peptide','Final site probability']), ascending=[True, True, True])
    #Repeat collapse for the updated decoy positions
    PSM_count=[]
    thresholded_PSM_count=[]
    for i in df.groupby(['Pep_pos','Modified Peptide'])['Supporting PSM count'].sum():
        PSM_count.append(i)
    for i in df.groupby(['Pep_pos','Modified Peptide'])['Supporting PSM Count (Thresholded 0.05 Global FLR)'].sum():
        thresholded_PSM_count.append(i)
    sig1=[]
    sig2=[]
    sig3=[]
    sig4=[]
    sig5=[]
    sig6=[]
    for i in df.groupby(['Pep_pos','Modified Peptide'])['PSM count 0.01<P<=0.05'].sum():
        sig1.append(i)
    for i in df.groupby(['Pep_pos','Modified Peptide'])['PSM count 0.05<P<=0.19'].sum():
        sig2.append(i)
    for i in df.groupby(['Pep_pos','Modified Peptide'])['PSM count 0.19<P<=0.81'].sum():
        sig3.append(i)
    for i in df.groupby(['Pep_pos','Modified Peptide'])['PSM count 0.81<P<=0.95'].sum():
        sig4.append(i)
    for i in df.groupby(['Pep_pos','Modified Peptide'])['PSM count 0.95<P<0.99'].sum():
        sig5.append(i)
    for i in df.groupby(['Pep_pos','Modified Peptide'])['PSM count P>=0.99'].sum():
        sig6.append(i)
    df=df.drop_duplicates(subset=(['Pep_pos','Modified Peptide']),keep='last',inplace=False)
    df['Supporting PSM count']=PSM_count
    df['Supporting PSM Count (Thresholded 0.05 Global FLR)']=thresholded_PSM_count
    df['PSM count 0.01<P<=0.05']=sig1
    df['PSM count 0.05<P<=0.19']=sig2
    df['PSM count 0.19<P<=0.81']=sig3
    df['PSM count 0.81<P<=0.95']=sig4
    df['PSM count 0.95<P<0.99']=sig5
    df['PSM count P>=0.99']=sig6
    df = df.reset_index(drop=True)


    #Recalculate FLR
    # df = df.sort_values(by=(['Final adjusted site probability']), ascending=[False])
    # df = df.reset_index(drop=True)
    # df['Count_binomial'] = (df.index) + 1
    # df['final_temp'] = 1 - df['Final adjusted site probability']
    # df['Binomial_final_prob_FLR'] = df['final_temp'].cumsum() / df['Count_binomial']
    # df['Binomial final adjusted q_value'] = df['Binomial_final_prob_FLR']
    # df['Binomial final adjusted q_value'] = df.iloc[::-1]['Binomial_final_prob_FLR'].cummin()

    organism=[]
    PMID=[]
    for i in range(len(df)):
        USI_temp=df.loc[i,'USI']
        PXD=USI_temp.split(":")[1]
        dataset_ID.append(PXD)
        organism.append(dict[PXD][1])
        PMID.append(dict[PXD][0])
    #print(df['Binomial_final_prob_q_value'])
    df['Site Passes Threshold']=np.where(df['Binomial final q value'] <= FLR_threshold, 1, 0)
    df['Dataset ID']=dataset_ID
    df['Site ID']=df.index
    df['Organism']=organism
    df['PubMedID']=PMID

    #Removed PTM_info, Decop_P,
    df = df[['Site ID', 'Proteins', 'Unmodified Sequence', 'Modified Peptide','Modification','Site position', 'Site probability', 'PSM Score',
            'Final site probability','Binomial final q value','Site Passes Threshold', 'Supporting PSM count','Supporting PSM Count (Thresholded 0.05 Global FLR)', 'Dataset ID', 'PRIDE dataset URL', 'PubMedID',
             'Sample ID', 'Organism', 'USI','PRIDE dataset URL (including the USI)','PSM count 0.01<P<=0.05','PSM count 0.05<P<=0.19','PSM count 0.19<P<=0.81','PSM count 0.81<P<=0.95','PSM count 0.95<P<0.99','PSM count P>=0.99']]
    df.to_csv("C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/Site_format_files/"+dataset_ID[0]+"_Site_peptidoform_centric_binomial_FLR_no_pA_bins_"+str(FLR_threshold)+"_updated2.tsv", sep="\t", index=False)

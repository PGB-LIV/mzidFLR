import pandas as pd
import os
import sys

#usage: python Uniprot_formats.py file_list.txt [optional: modification:target:decoy ie. phospho:STY:A]

folder_list_file = open(sys.argv[1], "r")
folder_list_file = folder_list_file.read()
folder_list = folder_list_file.replace('\n', ';').split(";")

dataset_list = []
for i in folder_list:
    index = [idx for idx, s in enumerate(i.split("/")) if 'PXD' in s][0]
    if i.split("/")[index] not in dataset_list:
        dataset_list.append(i.split("/")[index])

if not os.path.exists("Uniprot_files"):
    os.mkdir("Uniprot_files")

if len(sys.argv) > 1:
    decoy_aa = sys.argv[2].split(":")[2]
else:
    decoy_aa = "A"

for dataset in dataset_list:
    loc_full = "All_site_formats_Updated_peptidoform_decoy/" + dataset + "_merged_Site_Peptidoform_centric.tsv"
    df = pd.read_csv(loc_full, sep="\t")

    ### Peptidoform site format
    df = df.rename(columns={'All_Proteins': 'Proteins', 'All_PTM_protein_positions': 'Protein Modification Positions',
                            'Peptide_mod': 'Peptidoform', 'Peptide': 'Unmodified Sequence',
                            'p' + decoy_aa + '_q_value_BA': 'Site Q-Value', 'Score': 'PSM Probability',
                            'PTM Score': 'PTM Probability',
                            'PTM positions': 'Peptide Modification Position', 'Binomial_final_score': 'Final Probability',
                            'All_PTMs': 'Modifications', 'All_PTM_scores': 'Modification probabilities',
                            'All_PTM_positions': 'Modification positions',
                            '0.05FLR_threshold_count': 'PSM Count Passing Threshold [0.05]',
                            '0.01FLR_threshold_count': 'PSM Count Passing Threshold [0.01]',
                            '0.01<P<=0.05_count': 'opt_PSM count 0.01<P<=0.05',
                            '0.05<P<=0.19_count': 'opt_PSM count 0.05<P<=0.19',
                            '0.19<P<=0.81_count': 'opt_PSM count 0.19<P<=0.81',
                            '0.81<P<=0.95_count': 'opt_PSM count 0.81<P<=0.95',
                            '0.95<P<0.99_count': 'opt_PSM count 0.95<P<0.99', 'P>=0.99_count': 'opt_PSM count P>=0.99',
                            'USI': 'Universal Spectrum Identifier'})  # 'All_USI':'Supporting PSM count'


    df = df[['Peptidoform Site ID','Proteins','Unmodified Sequence','Peptidoform','Modification','Peptide Modification Position','Protein Modification Positions',
             'PSM Probability','PTM Probability','Final Probability','Site Q-Value','Site Passes Threshold [0.05]','Site Passes Threshold [0.01]','Decoy Peptide',
             'Decoy Modification Site','PSM Count Passing Threshold [0.05]','PSM Count Passing Threshold [0.01]','Source Dataset Identifier','Reanalysis Dataset Identifier',
             'PubMedIDs','Sample ID','Organism','Organism Part','Cell Line','Disease Information','Universal Spectrum Identifier','opt_PSM count 0.01<P<=0.05',
             'opt_PSM count 0.05<P<=0.19','opt_PSM count 0.19<P<=0.81','opt_PSM count 0.81<P<=0.95','opt_PSM count 0.95<P<0.99','opt_PSM count P>=0.99']]

    df.to_csv("Uniprot_files/" + dataset + "_merged_Site_Peptidoform_centric.tsv", sep="\t",
                            index=False)



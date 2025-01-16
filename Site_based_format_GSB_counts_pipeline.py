import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
# import requests
# import json
import sys
import time
import multiprocessing

"""
Requires csv file (meta.csv) with meta data, sdrf location (eg. "PRIDE/SDRFs/"), gold count and silver count - if meta or sdrf files not available, use "NA"
Alternatively, can use simple meta data files ("simple_meta/") when the meta data is the same for all files in a data set
Can also accept optional decoy and contam prefix as well as modification:target:decoy, if not specified "DECOY" and "CONTAM" prefixes will be used as default 
and phospho:STY:A used as search mod.
"""

# usage: python Site_based_format_GSB_counts_pipeline.py file_list.txt meta.csv PRIDE/SDRFs/ NA 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]
# OR
# python Site_based_format_GSB_counts_pipeline.py file_list.txt NA NA simple_meta/ 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]
# OR
# python Site_based_format_GSB_counts_pipeline.py file_list.txt NA NA NA 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]

# read in csv from FLR pipeline
# DECOY and CONTAM removed, filtered for contains Phospho, exploded for site based, binomial adjustment and collapsed by peptidoform

# TO DO: probably need some error handling here, to check using the right params
folder_list_file = open(sys.argv[1], "r")
folder_list_file = folder_list_file.read()
folder_list = folder_list_file.replace('\n', ';').split(";")

dataset_list = []
for i in folder_list:
    index = [idx for idx, s in enumerate(i.split("/")) if 'PXD' in s][0]
    if i.split("/")[index] not in dataset_list:
        dataset_list.append(i.split("/")[index])

meta_all = sys.argv[2]
if meta_all != "NA":
    with open(meta_all, 'r') as f:
        reader = csv.reader(f)
        dict_all = {rows[0]: rows[1:] for rows in reader}
    f.close()

gold_count = int(sys.argv[5])
silver_count = int(sys.argv[6])

if len(sys.argv) > 7:
    decoy_prefix = sys.argv[7]
    contam_prefix = sys.argv[8]
else:
    decoy_prefix = "DECOY"
    contam_prefix = "CONTAM"
if len(sys.argv) > 9:
    search_mod = sys.argv[9].split(":")[0]
    target_aa = sys.argv[9].split(":")[1]
    decoy_aa = sys.argv[9].split(":")[2]
else:
    decoy_aa = "A"
    search_mod = "Phospho"
    target_aa = "STY"

print(f"Using decoy prefix: {decoy_prefix}")
print(f"Using contam prefix: {contam_prefix}")
print(f"Using modification: {search_mod}")
print(f"Using target amino acids: {target_aa}")
print(f"Using decoy amino acid: {decoy_aa}")

num_workers=40

def ratio(df, targets, decoy):
    STY_count = 0
    for target in list(targets):
        T_count = df['Peptide'].str.count(target).sum()
        STY_count += T_count
    A_count = df['Peptide'].str.count(decoy).sum()
    STY_A_ratio = STY_count / A_count
    return STY_A_ratio

def per_folder (i,p, decoy_method):
    start_time = time.time()
    if os.path.isfile(output_location + "/" + i.replace("/", "_") + "_Site_" + p + "_centric.tsv"):
        print(f"{p}-Centric format for {i} already exists")
    else:
        print(f"Creating {p}-Centric format for {i}")
        if p == "PSM":
            file = i + "/FDR_0.01/binomial" + decoy_method + ".csv"
            df_temp = pd.read_csv(file)
            df = df_temp[
                ['All_Proteins', 'All_PTM_protein_positions', 'All_PTM_positions', 'All_PTMs', 'All_PTM_scores',
                 'All_USI', 'Peptide', 'Peptide_mod', 'Protein position', 'PTM Score', 'Score',
                 'PTM positions', 'PTM_info', 'Binomial_final_score', 'Source', 'USI',
                 'p' + decoy_aa + '_q_value']].copy()
        else:
            file = i + "/FDR_0.01/binomial_peptidoform_collapsed_FLR" + decoy_method + ".csv"
            df_temp = pd.read_csv(file)
            df = df_temp[
                ['All_Proteins', 'All_PTM_protein_positions', 'All_PTM_positions', 'All_PTMs', 'All_PTM_scores',
                 'All_USI', 'All_Source', 'Peptide', 'Protein', 'Protein position', 'PTM Score', 'Score', 'PTM positions',
                 'PTM_info', 'Binomial_final_score', 'Peptide_mod', 'Peptide_mod_pos', 'USI',
                 'p' + decoy_aa + '_q_value_BA','0.05FLR_threshold_count', '0.01FLR_threshold_count',
                 '0.01<P<=0.05_count', '0.05<P<=0.19_count', '0.19<P<=0.81_count', '0.81<P<=0.95_count',
                 '0.95<P<0.99_count', 'P>=0.99_count']].copy()
        df = df.reset_index(drop=True)

        df['Modification'] = search_mod
        PXD_index = [idx for idx, s in enumerate(i.split("/")) if 'PXD' in s][0]
        PXD = i.split("/")[PXD_index]

        # Simple meta data mode
        if sys.argv[4] != "NA":
            simple_meta_file = sys.argv[4] + "/" + PXD + ".tsv"
            simple_meta = pd.read_csv(simple_meta_file, sep="\t")
            for x in range(len(df)):
                protein_PTM = ""
                for z, y in zip(df.loc[x, 'All_Proteins'].split(":"),df.loc[x, 'All_PTM_protein_positions'].split(":")):
                    for a, b in zip(df.loc[x, 'All_PTM_positions'].split(";"), y.split(";")):
                        if int(a) == int(df.loc[x, 'PTM positions']):
                            protein_PTM += b + ";"
                    protein_PTM= protein_PTM[:-1] +":"
                df.loc[x, 'All_PTM_protein_positions'] = protein_PTM[:-1]
            df['PubMedIDs'] = simple_meta["PubMedIDs"].values[0]
            df['Sample ID'] = simple_meta["Sample ID"].values[0]
            df['Source Dataset Identifier'] = PXD
            df['Reanalysis Dataset Identifier'] = "NA"
            df['Organism'] = simple_meta["Organism"].values[0]
            df['Organism Part'] = simple_meta["Organism Part"].values[0]
            df['Cell Line'] = simple_meta["Cell Line"].values[0]
            df['Disease Information'] = simple_meta["Disease Information"].values[0]
        # SDRF mode
        elif sys.argv[3] != "NA":
            PMID = []
            Sample = []
            dataset_ID = []
            organism = []
            organism_part = []
            cell_line = []
            disease = []

            meta_file = sys.argv[3] + "/" + PXD + ".sdrf.tsv"
            meta = pd.read_csv(meta_file, sep="\t")
            meta['comment[data file]'] = meta['comment[data file]'].str.split(".").str[0]
            meta['comment[data file]'] = meta['comment[data file]'].str.replace(" ", "_")
            meta.set_index('comment[data file]', inplace=True)
            meta_dict = meta.to_dict()

            for x in range(len(df)):
                protein_PTM = ""
                for z, y in zip(df.loc[x, 'All_Proteins'].split(":"),
                                df.loc[x, 'All_PTM_protein_positions'].split(":")):
                    for a, b in zip(df.loc[x, 'All_PTM_positions'].split(";"), y.split(";")):
                        if int(a) == int(df.loc[x, 'PTM positions']):
                            protein_PTM += b + ";"
                    protein_PTM=protein_PTM[:-1]+":"
                df.loc[x, 'All_PTM_protein_positions'] = protein_PTM[:-1]

                if p == "PSM":
                    source = str(df.loc[x, "Source"])
                    source = source.replace("_raw", "")
                    source = source.replace("(", "")
                    source = source.replace(")", "")
                    df.loc[x,"Source_new"]=source

            if p=="PSM":
                df['Sample ID'] = df['Source_new'].map(meta_dict['source name'])
                df['Organism'] = df['Source_new'].map(meta_dict['characteristics[organism]'])

                if meta_all != "NA":
                    df['PubMedIDs'] = dict_all[PXD][0]
                else:
                    df['PubMedIDs'] = "NA"

                df['Source Dataset Identifier'] = PXD
                df['Reanalysis Dataset Identifier'] = "NA"
                try:
                    df['Organism Part'] = df['Source_new'].map(meta_dict["characteristics[organism part]"])
                except:
                    df['Organism Part'] = "NA"
                df['Cell Line'] = df['Source_new'].map(meta_dict["characteristics[disease]"])
                try:
                    df['Disease Information'] = df['Source_new'].map(meta_dict["characteristics[disease]"])
                except:
                    df['Disease Information'] = "NA"
                df=df.drop(['Source_new'], axis=1)
            else:
                for x in range(len(df)):
                    source_all = df.loc[x, "All_Source"]
                    for source in source_all.split(";"):
                        source = source.replace("_raw", "")
                        source = source.replace("(", "")
                        source = source.replace(")", "")
                        PMID_temp = []
                        Sample_temp = []
                        organism_temp = []
                        organism_part_temp = []
                        cell_line_temp = []
                        disease_temp = []

                        Sample_temp.append(meta_dict["source name"][source])
                        organism_temp.append(meta_dict["characteristics[organism]"][source])
                        try:
                            organism_part_temp.append(meta_dict["characteristics[organism part]"][source])
                        except:
                            organism_part_temp.append("NA")
                        cell_line_temp.append(meta_dict["characteristics[cell type]"][source])
                        try:
                            disease_temp.append(meta_dict["characteristics[disease]"][source])
                        except:
                            disease_temp.append("NA")

                        if meta_all != "NA":
                            PMID_temp.append(dict_all[PXD][0])
                        else:
                            PMID_temp.append("NA")

                    PMID.append(";".join(map(str, list(set(PMID_temp)))))
                    Sample.append(";".join(map(str, list(set(Sample_temp)))))
                    dataset_ID.append(PXD)
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

        # No meta data mode
        else:
            for x in range(len(df)):
                protein_PTM = ""
                for z, y in zip(df.loc[x, 'All_Proteins'].split(":"),
                                df.loc[x, 'All_PTM_protein_positions'].split(":")):
                    for a, b in zip(df.loc[x, 'All_PTM_positions'].split(";"), y.split(";")):
                        if int(a) == int(df.loc[x, 'PTM positions']):
                            protein_PTM += b + ";"
                    protein_PTM=protein_PTM[:-1]+":"
                df.loc[x, 'All_PTM_protein_positions'] = protein_PTM[:-1]
            df['PubMedIDs'] = "NA"
            df['Sample ID'] = "NA"
            df['Source Dataset Identifier'] = PXD
            df['Reanalysis Dataset Identifier'] = "NA"
            df['Organism'] = "NA"
            df['Organism Part'] = "NA"
            df['Cell Line'] = "NA"
            df['Disease Information'] = "NA"

        #####

        df['Pep_pos'] = df['Peptide'] + "-" + df['PTM positions'].astype(str)
        if p == "PSM":
            df['Site Passes Threshold [0.05]'] = np.where(df['p' + decoy_aa + '_q_value'] <= 0.05, 1, 0)
            df['Site Passes Threshold [0.01]'] = np.where(df['p' + decoy_aa + '_q_value'] <= 0.01, 1, 0)
        else:
            df['Site Passes Threshold [0.05]'] = np.where(df['p' + decoy_aa + '_q_value_BA'] <= 0.05, 1, 0)
            df['Site Passes Threshold [0.01]'] = np.where(df['p' + decoy_aa + '_q_value_BA'] <= 0.01, 1, 0)
        df['Protein_count'] = df['All_Proteins'].str.count(":") + 1
        df['Decoy_count'] = df['All_Proteins'].str.count(decoy_prefix)
        df['Contam_count'] = df['All_Proteins'].str.count(contam_prefix)
        # Decoy peptide = all proteins are decoys or contam
        df['Decoy Peptide'] = np.where(df['Protein_count'] == df['Decoy_count'] + df['Contam_count'], 1, 0)
        df['Decoy Modification Site'] = df.apply(lambda x: r(x.Pep_pos.split("-")[0], x.Pep_pos.split("-")[1]), axis=1)
        if p == "PSM":
            df['PSM Site ID'] = df.index
            df['PSM Count Passing Threshold [0.05]'] = "1"
        else:
            df['Peptidoform Site ID'] = df.index

        # peptidoform without reagent labels?
        df['Peptide'] = df['Peptide'].str.replace(r'\[TMT(.*?)\]', '', regex="True")
        df['Peptide'] = df['Peptide'].str.replace(r'\[iTRAQ.*?\]', '', regex="True")

        df = df.replace("not applicable", "NA")
        df.to_csv(output_location + "/" + i.replace("/", "_") + "_Site_" + p + "_centric.tsv", sep="\t", index=False)

        print(f"{p}-Centric format for {i} complete")
        print("--- %s seconds ---" % (time.time() - start_time))

        # #Test for valid USI
        # #Won't work unless can access url
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

def per_dataset (dataset, p):
    if os.path.isfile(output_location + "/" + dataset + "_merged_Site_" + p + "_centric.tsv"):
        print(f"Merged {p} files for {dataset} already exists")
    else:
        print(f"creating merged {p} files for {dataset}")
        folder_list_temp = [x for x in folder_list if dataset in x]
        for loc in folder_list_temp:
            loc_full = output_location + "/" + loc.replace("/", "_") + "_Site_" + p + "_centric.tsv"
            df = pd.read_csv(loc_full, sep="\t", dtype={'All_PTM_protein_positions': str})
            if loc == folder_list_temp[0]:
                df_temp = df
            else:
                df_temp = pd.concat([df_temp, df])
        if p == "PSM":
            df_temp = df_temp.reset_index(drop=True)
            df_temp.to_csv(output_location + "/" + dataset + "_merged_Site_PSM_centric.tsv", sep="\t", index=False)
        else:
            df_temp = df_temp.sort_values(
                ['Peptide_mod_pos', 'p' + decoy_aa + '_q_value_BA', 'Binomial_final_score'],
                ascending=[True, False, True])
            df_temp = df_temp.drop_duplicates(subset=('Peptide_mod_pos'), keep="last", inplace=False)
            df_temp = df_temp.reset_index()
            df_temp['Peptidoform Site ID'] = df_temp.index

            df_merge = df_temp.drop(['index'], axis=1)
            df_merge.to_csv(output_location + "/" + dataset + "_merged_Site_Peptidoform_centric.tsv", sep="\t",
                            index=False)

r = lambda x, y: 0 if x[int(y) - 1] != decoy_aa or int(y) == 0 else 1

for decoy_method in ["", "_peptidoform_decoy", "_site_decoy"]:
    file=folder_list[0]+"/FDR_0.01/"+"binomial"+decoy_method+".csv"
    if os.path.isfile(file):
        output_location = "All_site_formats_Updated" + decoy_method
        if not os.path.exists(output_location):
            os.mkdir(output_location)
        for p in ["PSM","Peptidoform"]:
            with multiprocessing.Pool(processes=num_workers) as pool:
                results = pool.starmap(per_folder, [(folder, p, decoy_method) for folder in folder_list])

            #Merged files
            print(f"Merging {p}-centric files per Dataset")
            with multiprocessing.Pool(processes=num_workers) as pool:
                results = pool.starmap(per_dataset, [(dataset, p) for dataset in dataset_list])

            if os.path.isfile(output_location + "/all_datasets_merged_Site_" + p + "_centric.tsv"):
                print(f"Merged {p}-centric files for all Datasets already exists")
            else:
                print(f"Merging {p}-centric files for all Datasets")
                for dataset in dataset_list:
                    loc = output_location + "/" + dataset + "_merged_Site_"+p+"_centric.tsv"
                    df_temp=pd.read_csv(loc, sep="\t")
                    if dataset == dataset_list[0]:
                        df_all = df_temp
                    else:
                        df_all = pd.concat([df_temp, df_all])
                df_all.to_csv(output_location + "/all_datasets_merged_Site_"+p+"_centric.tsv", sep="\t", index=False)
            if p=="Peptidoform":
                if os.path.isfile(output_location + "/all_datasets_merged_Site_" + p + "_centric_Uniprot.tsv"):
                    print(f"Uniprot format for merged {p}-centric files for all Datasets already exists")
                else:
                    print(f"Creating Uniprot {p}-centric merged format for all Datasets")
                    df_all = pd.read_csv(output_location + "/all_datasets_merged_Site_"+p+"_centric.tsv", sep="\t")
                    df_all = df_all.rename(columns={'All_Proteins': 'Proteins',
                                            'All_PTM_protein_positions': 'Protein Modification Positions',
                                            'Peptide_mod': 'Peptidoform', 'Peptide': 'Unmodified Sequence',
                                            'p' + decoy_aa + '_q_value_BA': 'Site Q-Value',
                                            'Score': 'PSM Probability',
                                            'PTM Score': 'PTM Probability',
                                            'PTM positions': 'Peptide Modification Position',
                                            'Binomial_final_score': 'Final Probability',
                                            'All_PTMs': 'Modifications',
                                            'All_PTM_scores': 'Modification probabilities',
                                            'All_PTM_positions': 'Modification positions',
                                            '0.05FLR_threshold_count': 'PSM Count Passing Threshold [0.05]',
                                            '0.01FLR_threshold_count': 'PSM Count Passing Threshold [0.01]',
                                            '0.01<P<=0.05_count': 'opt_PSM count 0.01<P<=0.05',
                                            '0.05<P<=0.19_count': 'opt_PSM count 0.05<P<=0.19',
                                            '0.19<P<=0.81_count': 'opt_PSM count 0.19<P<=0.81',
                                            '0.81<P<=0.95_count': 'opt_PSM count 0.81<P<=0.95',
                                            '0.95<P<0.99_count': 'opt_PSM count 0.95<P<0.99',
                                            'P>=0.99_count': 'opt_PSM count P>=0.99',
                                            'USI': 'Universal Spectrum Identifier'})

                    df_all = df_all[['Peptidoform Site ID', 'Proteins', 'Unmodified Sequence', 'Peptidoform', 'Modification','Peptide Modification Position',
                                     'Protein Modification Positions','PSM Probability', 'PTM Probability', 'Final Probability', 'Site Q-Value','Site Passes Threshold [0.05]',
                                     'Site Passes Threshold [0.01]', 'Decoy Peptide','Decoy Modification Site', 'PSM Count Passing Threshold [0.05]',
                                     'PSM Count Passing Threshold [0.01]', 'Source Dataset Identifier','Reanalysis Dataset Identifier','PubMedIDs', 'Sample ID', 'Organism',
                                     'Organism Part', 'Cell Line', 'Disease Information','Universal Spectrum Identifier', 'opt_PSM count 0.01<P<=0.05','opt_PSM count 0.05<P<=0.19',
                                     'opt_PSM count 0.19<P<=0.81', 'opt_PSM count 0.81<P<=0.95','opt_PSM count 0.95<P<0.99', 'opt_PSM count P>=0.99']]
                    df_all.to_csv(output_location + "/all_datasets_merged_Site_"+p+"_centric_Uniprot.tsv", sep="\t", index=False)

        print("PSM and Peptidoform formats done, creating GSB counts")
        print("Gold threshold: " + str(gold_count) + "\nSilver threshold: " + str(silver_count))
        flr_filter = 0.05

        df = pd.read_csv(output_location + "/all_datasets_merged_Site_Peptidoform_centric.tsv", sep="\t")
        # Calculate STY:A ratio from the PSM centric format - for calculating FLR at GSB levels
        STY_ratio = ratio(df, target_aa, decoy_aa)
        print(f"Using target:decoy ratio of: {STY_ratio}")

        ### GSB counts
        # filter for sites seen multiple times across dataset studies
        for m in ["single", "all"]:
            print(f"Creating GSB format mapped to {m} proteins")
            for dataset in dataset_list:
                loc_full = output_location + "/" + dataset + "_merged_Site_Peptidoform_centric.tsv"
                df = pd.read_csv(loc_full, sep="\t")
                df = df.loc[df['p' + decoy_aa + '_q_value_BA'] <= flr_filter]
                PSM_threshold = dict(df.groupby('Peptide_mod_pos')['0.05FLR_threshold_count'].sum())

                df = df.sort_values(['Peptide_mod_pos', 'p' + decoy_aa + '_q_value_BA', 'Binomial_final_score'],
                                         ascending=[True, False, True])
                df = df.drop_duplicates(subset=('Peptide_mod_pos'), keep="last", inplace=False)

                df['0.05FLR_threshold_count'] = df['Peptide_mod_pos'].map(PSM_threshold)

                if m == "all":
                    if len(df) >= 1:
                        df['PTM_residue'] = df.apply(lambda x: x['Peptide'][x['PTM positions'] - 1], axis=1)
                        df['All_PTM_protein_positions'] = df['All_PTM_protein_positions'].str.split(':')
                        df['All_Proteins'] = df['All_Proteins'].str.split(':')
                        df = df.explode(['All_PTM_protein_positions', 'All_Proteins'])
                        df["Protein-pos"] = df['All_Proteins'] + "-" + df['All_PTM_protein_positions']
                    else:
                        df['PTM_residue'] = ""
                else:
                    if len(df) >= 1:
                        df['PTM_residue'] = df.apply(lambda x: x['Peptide'][x['PTM positions'] - 1], axis=1)
                        df['Protein-pos'] = df['Protein'] + "-" + df['Protein position'].astype(str)
                    else:
                        df['PTM_residue'] = ""
                df['PTM_residue'] = np.where(df['PTM positions'] == 0, "N-term", df['PTM_residue'])
                df = df.sort_values(['Protein-pos', 'p' + decoy_aa + '_q_value_BA', 'Binomial_final_score'],
                                              ascending=[True, False, True])

                df['Protein_pos_res'] = df['Protein-pos'] + "_" + df['PTM_residue']

                PSM_threshold_2 = dict(df.groupby('Protein_pos_res')['0.05FLR_threshold_count'].sum())

                df = df.drop_duplicates(subset=('Protein_pos_res'), keep="last", inplace=False)

                # column for all PSM counts at 5%FLR
                df['Sum_of_PSM_counts(5%FLR)'] = df['Protein_pos_res'].map(PSM_threshold_2)

                df['PXD'] = dataset
                if dataset == dataset_list[0]:
                    df_counts = df
                else:
                    df_counts = pd.concat([df_counts, df])

                df_temp = df[
                    ['Peptide_mod_pos', 'p' + decoy_aa + '_q_value_BA', 'Binomial_final_score', 'Protein_pos_res',
                     'Sum_of_PSM_counts(5%FLR)', '0.05FLR_threshold_count']]
                df_temp.rename(columns={'p' + decoy_aa + '_q_value_BA': dataset + "_FLR",
                                        'Binomial_final_score': dataset + '_BinomialScore',
                                        'Peptide_mod_pos': dataset + "_peptide_mod_pos",
                                        'Sum_of_PSM_counts(5%FLR)': dataset + "_Sum_of_PSM_counts(5%FLR)",
                                        '0.05FLR_threshold_count': dataset + "_peptidoform_PSMcount(5%FLR)"}, inplace=True)
                df_temp = df_temp.set_index('Protein_pos_res')

                if dataset == dataset_list[0]:
                    df_final = df_temp
                else:
                    df_final = pd.concat([df_final, df_temp], axis=1)

            column_list = [PXD + "_Sum_of_PSM_counts(5%FLR)" for PXD in dataset_list]
            df_final["Sum_of_PSM_counts(5%FLR)"] = df_final[column_list].sum(axis=1)

            cols = df_final.columns.tolist()
            cols = [cols[-1]] + cols[:-1]
            df_final = df_final[cols]
            df_final = df_final.drop(column_list, axis=1)

            df_final = df_final[~df_final.index.str.contains(decoy_prefix)]
            df_final = df_final[~df_final.index.str.contains(contam_prefix)]

            df_counts['PTM_residue'] = df_counts['Protein_pos_res'].str.rsplit("_", n=1).str[-1]
            counts_res = pd.crosstab(df_counts['PTM_residue'], df_counts['PXD']).replace(0,
                                                                                         np.nan).stack().reset_index().rename(
                columns={0: 'Count'})
            print(counts_res)
            counts_res.to_csv(output_location + "/G" + str(gold_count) + "S" + str(silver_count) + "B_" + str(
                flr_filter) + "_Residue_counts_" + m+".csv", index=False)

            for i in df_final.index.values.tolist():
                df_final.loc[i, 'Protein'] = i.rsplit("-", 1)[0]
                df_final.loc[i, 'Protein_pos'] = i.rsplit("-", 1)[1].split("_")[0]
                df_final.loc[i, 'PTM_residue'] = i.rsplit("_", 1)[-1]
                flr1_count = 0
                for dataset in dataset_list:
                    FLR_count_all = 0
                    FLR_cols = [x for x in dataset_list if dataset in x]
                    for FLR_col in FLR_cols:
                        if df_final.loc[i, FLR_col + "_FLR"] != "N/A":
                            if float(df_final.loc[i, FLR_col + "_FLR"]) <= 0.01:
                                FLR_count_all += 1
                    if FLR_count_all != 0:
                        flr1_count += 1
                if flr1_count >= gold_count:
                    df_final.loc[i, 'PTM_FLR_category'] = "Gold"
                elif flr1_count >= silver_count:
                    df_final.loc[i, 'PTM_FLR_category'] = "Silver"
                else:
                    df_final.loc[i, 'PTM_FLR_category'] = "Bronze"

            cols = list(df_final.columns.values)
            cols.pop(cols.index('Protein'))
            cols.pop(cols.index('Protein_pos'))
            cols.pop(cols.index('PTM_residue'))
            cols.pop(cols.index('PTM_FLR_category'))

            df_final['Decoy_mod'] = np.where(df_final['PTM_residue'] == decoy_aa, 1, 0)

            df_final = df_final[['Protein', 'Protein_pos', 'PTM_residue', 'Decoy_mod', 'PTM_FLR_category'] + cols]

            # replace 0 PSM counts -> 1
            column_list = [PXD + "_peptidoform_PSMcount(5%FLR)" for PXD in dataset_list]
            column_list.append("Sum_of_PSM_counts(5%FLR)")
            df_final[column_list] = df_final[column_list].replace(0, 1)

            df_final.to_csv(output_location + "/G" + str(gold_count) + "S" + str(silver_count) + "B_" + str(
                flr_filter) + "_protein_pos" + decoy_method +"_"+ m +"_prot_mapping.csv", index=False)

            counts = pd.crosstab(df_final.PTM_residue, df_final.PTM_FLR_category).replace(0,
                                                                                          np.nan).stack().reset_index().rename(
                columns={0: 'Count'})

            # calculate FLR estimate using STY:A ratio from PSM centric format
            f = lambda x, y, z: np.nan if x != decoy_aa else str(
                round(z * STY_ratio / counts.loc[counts['PTM_FLR_category'] == y, 'Count'].sum() * 100, 2)) + "%"

            counts['FLR'] = counts.apply(lambda x: f(x.PTM_residue, x.PTM_FLR_category, x.Count), axis=1)
            counts['FLR'] = counts.groupby('PTM_FLR_category')['FLR'].ffill()
            counts['PTM_FLR_category_FLR'] = counts['PTM_FLR_category'] + " - " + counts['FLR']
            counts['PTM_FLR_category'] = pd.Categorical(counts['PTM_FLR_category'], ["Bronze", "Silver", "Gold"])
            counts = counts.sort_values("PTM_FLR_category")
            print(counts)

            fig, axes = plt.subplots(ncols=3)
            for i, (name, group) in enumerate(counts.groupby("PTM_FLR_category_FLR", sort=False)):
                axes[i].set_title(name)
                group.plot(kind="bar", x="PTM_residue", y="Count", ax=axes[i], legend=False)
                axes[i].set_ylabel("count")
                axes[i].set_xlabel("")
            plt.tight_layout()
            # plt.show()
            plt.savefig(output_location + "/G" + str(gold_count) + "S" + str(silver_count) + "B_" + str(
                flr_filter) + "_protein_pos_categories" + decoy_method +"_"+ m +"_prot_mapping.png", dpi=300)

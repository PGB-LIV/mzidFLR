import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
#import requests
#import json
import sys

"""
Requires csv file (meta.csv) with meta data, sdrf location (eg. "PRIDE/SDRFs/"), gold count and silver count - if meta or sdrf files not available, use "NA"
Can also accept optional decoy and contam prefix, if not specified "DECOY" and "CONTAM" will be used as default
"""

# usage: python Site_based_format_GSB_counts_pipeline.py file_list.txt meta.csv PRIDE/SDRFs/ 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM]
# OR
# python Site_based_format_GSB_counts_pipeline.py file_list.txt NA NA 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM]

#read in csv from FLR pipeline
#DECOY and CONTAM removed, filtered for contains Phospho, exploded for site based, binomial adjustment and collapsed by peptidoform

folder_list_file = open(sys.argv[1],"r")
folder_list_file = folder_list_file.read()
folder_list = folder_list_file.replace('\n', ';').split(";")

dataset_list=[]
for i in folder_list:
    index = [idx for idx, s in enumerate(i.split("/")) if 'PXD' in s][0]
    if i.split("/")[index] not in dataset_list:
        dataset_list.append(i.split("/")[index])

meta_all = sys.argv[2]
#meta_all=wd+"/meta.csv"
if meta_all!="NA":
    with open(meta_all,'r') as f:
        reader = csv.reader(f)
        dict_all = {rows[0]:rows[1:] for rows in reader}
    f.close()

gold_count = int(sys.argv[4])
silver_count = int(sys.argv[5])

if len(sys.argv)>6:
    decoy_prefix = sys.argv[6]
    contam_prefix = sys.argv[7]
else:
    decoy_prefix = "DECOY"
    contam_prefix = "CONTAM"

print(f"Using decoy prefix: {decoy_prefix}")
print(f"Using contam prefix: {contam_prefix}")

r = lambda x, y: 0 if x[int(y) - 1] != "A" or int(y)==0 else 1

for decoy_method in ["","_peptidoform_decoy","_site_decoy"]:
    file=folder_list[0]+"/FDR_0.01/"+"binomial"+decoy_method+".csv"
    if os.path.isfile(file):
        output_location="All_site_formats"+decoy_method
        if not os.path.exists(output_location):
            os.mkdir(output_location)
        ## Peptide site format
        for i in folder_list:
            if os.path.isfile(output_location + "/"+i.replace("/","_")+"_Site_PSM_centric.tsv"):
                continue
            print(i)
            file=i+"/FDR_0.01/binomial"+decoy_method+".csv"
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

                #try:
                #    PXD = dict["comment[proteomexchange accession number]"][source]
                #except:
                #    PXD="NA"

                PXD_index = [idx for idx, s in enumerate(i.split("/")) if 'PXD' in s][0]
                PXD=i.split("/")[PXD_index]

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

            df.to_csv(output_location + "/"+i.replace("/","_")+"_Site_PSM_centric.tsv", sep="\t", index=False)

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
            if os.path.isfile(output_location + "/"+dataset+"_merged_Site_PSM_centric.tsv"):
                continue
            folder_list_temp=[x for x in folder_list if dataset in x]
            for loc in folder_list_temp:
                loc_full=output_location + "/"+loc.replace("/","_")+"_Site_PSM_centric.tsv"
                df=pd.read_csv(loc_full, sep="\t", dtype={'Protein Modification Positions': str})
                if loc==folder_list_temp[0]:
                    df_temp=df
                else:
                    df_temp=pd.concat([df_temp,df])
            df_temp=df_temp.reset_index(drop=True)
            df_temp['PSM Site ID']=df_temp.index

            df_temp.to_csv(output_location + "/"+dataset+"_merged_Site_PSM_centric.tsv", sep="\t", index=False)

            if dataset==dataset_list[0]:
                df_all=df_temp
            else:
                df_all=pd.concat([df_temp,df_all])
        if not os.path.isfile(output_location + "/all_datasets_merged_Site_PSM_centric.tsv"):
            df_all.to_csv(output_location + "/all_datasets_merged_Site_PSM_centric.tsv", sep="\t", index=False)


        ### Peptidoform site format
        print("PSM format complete, creating Peptidoform format")
        for i in folder_list:
            if os.path.isfile(output_location + "/"+i.replace("/","_")+"_Site_Peptidoform_centric.tsv"):
                continue
            print(i)
            file=i+"/FDR_0.01/binomial_peptidoform_collapsed_FLR"+decoy_method+".csv"
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

                PXD_index = [idx for idx, s in enumerate(i.split("/")) if 'PXD' in s][0]
                PXD=i.split("/")[PXD_index]

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
                        #PXD = dict["comment[proteomexchange accession number]"][source]
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
                #dataset_ID.append(";".join(map(str, list(set(dataset_ID_temp)))))
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
            #####

            df = df[['Peptidoform Site ID', 'Proteins', 'Unmodified Sequence', 'Peptidoform','Modification','Peptide Modification Position','Protein Modification Positions', 'PSM Probability', 'PTM Probability',
                     'Final Probability','Site Q-Value','Site Passes Threshold [0.05]','Site Passes Threshold [0.01]', 'Decoy Peptide','Decoy Modification Site',
                     'PSM Count Passing Threshold [0.05]','PSM Count Passing Threshold [0.01]', 'Source Dataset Identifier', 'Reanalysis Dataset Identifier','PubMedIDs',
                     'Sample ID', 'Organism', 'Organism Part', 'Cell Line', 'Disease Information', 'Universal Spectrum Identifier','opt_PSM count 0.01<P<=0.05','opt_PSM count 0.05<P<=0.19',
                     'opt_PSM count 0.19<P<=0.81','opt_PSM count 0.81<P<=0.95','opt_PSM count 0.95<P<0.99','opt_PSM count P>=0.99']]
            df = df.replace("not applicable","NA")
            df.to_csv(output_location + "/"+i.replace("/","_")+"_Site_Peptidoform_centric.tsv", sep="\t", index=False)

        #Merged files
        for dataset in dataset_list:
            if os.path.isfile(output_location + "/"+dataset+"_merged_Site_Peptidoform_centric.tsv"):
                continue
            folder_list_temp=[x for x in folder_list if dataset in x]
            for loc in folder_list_temp:
                loc_full=output_location + "/"+loc.replace("/","_")+"_Site_Peptidoform_centric.tsv"
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
            df_merge.to_csv(output_location + "/"+dataset+"_merged_Site_Peptidoform_centric.tsv", sep="\t", index=False)

            if dataset==dataset_list[0]:
                df_all=df_merge
            else:
                df_all=pd.concat([df_merge,df_all])

        if not os.path.isfile(output_location + "/all_datasets_merged_Site_Peptidoform_centric.tsv"):
            df_all.to_csv(output_location + "/all_datasets_merged_Site_Peptidoform_centric.tsv", sep="\t", index=False)


        print("Peptidoform formats done, creating GSB counts")
        print("Gold threshold: "+ str(gold_count) + "\nSilver threshold: "+ str(silver_count))
        flr_filter=0.05
        ### GSB counts
        #filter for sites seen multiple times across dataset studies
        for choice in ["","_no_choice"]:
            if not os.path.isfile(folder_list[0]+"/FDR_0.01/binomial_peptidoform_collapsed_FLR"+decoy_method+choice+".csv"):
                continue
            else:
                for dataset in dataset_list:
                    folder_list_temp=[x for x in folder_list if dataset in x]
                    for loc in folder_list_temp:
                        loc_full=loc+"/FDR_0.01/binomial_peptidoform_collapsed_FLR"+decoy_method+choice+".csv"
                        df=pd.read_csv(loc_full)

                        df=df.loc[df['pA_q_value_BA']<=flr_filter]

                        df['Experiment']=loc

                        if loc==folder_list_temp[0]:
                            df_temp=df
                        else:
                            df_temp=pd.concat([df_temp,df])

                    PSM_threshold=dict(df_temp.groupby('Peptide_mod_pos')['0.05FLR_threshold_count'].sum())

                    df_temp=df_temp.sort_values(['Peptide_mod_pos','Binomial_final_score','pA_q_value_BA'],ascending=[True,True,False])
                    df_temp=df_temp.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)

                    df_temp['0.05FLR_threshold_count']=df_temp['Peptide_mod_pos'].map(PSM_threshold)

                    df_temp=df_temp.sort_values(['Protein-pos','Binomial_final_score','pA_q_value_BA'],ascending=[True,True,False])

                    if len(df_temp)>=1:
                        df_temp['PTM_residue']=df_temp.apply(lambda x: x['Peptide'][x['PTM positions']-1],axis=1)
                    else:
                        df_temp['PTM_residue']=""
                    df_temp['Protein_pos_res']=df_temp['Protein-pos']+"_"+df_temp['PTM_residue']

                    PSM_threshold_2=dict(df_temp.groupby('Protein_pos_res')['0.05FLR_threshold_count'].sum())

                    df_temp=df_temp.drop_duplicates(subset=('Protein_pos_res'),keep="last",inplace=False)

                    #column for all PSM counts at 5%FLR
                    df_temp['Sum_of_PSM_counts(5%FLR)']=df_temp['Protein_pos_res'].map(PSM_threshold_2)

                    df_temp['PXD']=dataset
                    if dataset==dataset_list[0]:
                        df_counts=df_temp
                    else:
                        df_counts=pd.concat([df_counts,df_temp])

                    df_temp=df_temp[['Peptide_mod_pos', 'pA_q_value_BA','Binomial_final_score', 'Protein_pos_res','Sum_of_PSM_counts(5%FLR)','0.05FLR_threshold_count']]
                    df_temp.rename(columns = {'pA_q_value_BA':dataset+"_FLR",'Binomial_final_score':dataset+'_BinomialScore','Peptide_mod_pos':dataset+"_peptide_mod_pos",'Sum_of_PSM_counts(5%FLR)':dataset+"_Sum_of_PSM_counts(5%FLR)",'0.05FLR_threshold_count':dataset+"_peptidoform_PSMcount(5%FLR)"}, inplace = True)
                    df_temp=df_temp.set_index('Protein_pos_res')

                    if dataset==dataset_list[0]:
                        df_final=df_temp
                    else:
                        df_final=pd.concat([df_final,df_temp],axis=1)

                column_list=[PXD+"_Sum_of_PSM_counts(5%FLR)" for PXD in dataset_list]
                df_final["Sum_of_PSM_counts(5%FLR)"]=df_final[column_list].sum(axis=1)

                cols=df_final.columns.tolist()
                cols=[cols[-1]]+cols[:-1]
                df_final=df_final[cols]
                df_final=df_final.drop(column_list, axis=1)

                df_counts['PTM_residue']=df_counts['Protein_pos_res'].str.rsplit("_",n=1).str[-1]
                counts_res=pd.crosstab(df_counts['PTM_residue'],df_counts['Experiment']).replace(0,np.nan).stack().reset_index().rename(columns={0:'Count'})
                print(counts_res)
                counts_res.to_csv(output_location + "/G"+str(gold_count)+"S"+str(silver_count)+"B_"+str(flr_filter)+"_Residue_counts_"+decoy_method+choice+".csv", index=False)

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

                df_final.to_csv(output_location + "/G"+str(gold_count)+"S"+str(silver_count)+"B_"+str(flr_filter)+"_protein_pos"+decoy_method+choice+".csv", index=False)

                counts=pd.crosstab(df_final.PTM_residue,df_final.PTM_FLR_category).replace(0,np.nan).stack().reset_index().rename(columns={0:'Count'})

                f = lambda x, y, z : np.nan if x!="A" else str(round(z*3/counts.loc[counts['PTM_FLR_category']==y,'Count'].sum()*100,2))+"%"

                counts['FLR']=counts.apply(lambda x: f(x.PTM_residue, x.PTM_FLR_category, x.Count), axis=1)
                counts['FLR']=counts.groupby('PTM_FLR_category')['FLR'].ffill()
                counts['PTM_FLR_category_FLR']=counts['PTM_FLR_category'] + " - " + counts['FLR']
                counts['PTM_FLR_category'] = pd.Categorical(counts['PTM_FLR_category'], ["Bronze", "Silver", "Gold"])
                counts=counts.sort_values("PTM_FLR_category")
                print(choice)
                print(counts)

                fig, axes = plt.subplots(ncols=3)
                for i, (name, group) in enumerate(counts.groupby("PTM_FLR_category_FLR", sort=False)):
                    axes[i].set_title(name)
                    group.plot(kind="bar", x = "PTM_residue", y="Count", ax=axes[i], legend=False)
                    axes[i].set_ylabel("count")
                    axes[i].set_xlabel("")
                plt.tight_layout()
                #plt.show()
                plt.savefig(output_location + "/G"+str(gold_count)+"S"+str(silver_count)+"B_"+str(flr_filter)+"_protein_pos_categories"+decoy_method+choice+".png",dpi=300)


        print("Creating GSB counts, mapped to all proteins")
        print("Gold threshold: "+ str(gold_count) + "\nSilver threshold: "+ str(silver_count))
        flr_filter=0.05
        ### GSB counts
        #filter for sites seen multiple times across dataset studies
        for choice in ["","_no_choice"]:
            if not os.path.isfile(folder_list[0]+"/FDR_0.01/binomial_peptidoform_collapsed_FLR"+decoy_method+choice+".csv"):
                continue
            else:
                for dataset in dataset_list:
                    folder_list_temp=[x for x in folder_list if dataset in x]
                    for loc in folder_list_temp:
                        loc_full=loc+"/FDR_0.01/binomial_peptidoform_collapsed_FLR"+decoy_method+choice+".csv"
                        df=pd.read_csv(loc_full)

                        df=df.loc[df['pA_q_value_BA']<=flr_filter]

                        if loc==folder_list_temp[0]:
                            df_temp=df
                        else:
                            df_temp=pd.concat([df_temp,df])


                    PSM_threshold=dict(df_temp.groupby('Peptide_mod_pos')['0.05FLR_threshold_count'].sum())

                    df_temp=df_temp.sort_values(['Peptide_mod_pos','Binomial_final_score','pA_q_value_BA'],ascending=[True,True,False])
                    df_temp=df_temp.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)

                    df_temp['0.05FLR_threshold_count']=df_temp['Peptide_mod_pos'].map(PSM_threshold)

                    if len(df_temp)>=1:
                        df_temp['PTM_residue']=df_temp.apply(lambda x: x['Peptide'][x['PTM positions']-1],axis=1)
                        df_temp['All_PTM_protein_positions']=df_temp['All_PTM_protein_positions'].str.split(':')
                        df_temp['All_Proteins']=df_temp['All_Proteins'].str.split(':')
                        df_temp = df_temp.explode(['All_PTM_protein_positions','All_Proteins'])

                        p = lambda x, y ,z: z[y.index(str(x))] #x=mod_pos, y=All_PTM_positions, z=All_PTM_protein_positions

                        df_temp["Protein_PTM_pos"]=df_temp.apply(lambda x: p(x['PTM positions'],x['All_PTM_positions'].split(";"), x['All_PTM_protein_positions'].split(";")), axis=1)
                        df_temp["Protein-pos"]=df_temp['All_Proteins']+"-"+df_temp['Protein_PTM_pos']

                    else:
                        df_temp['PTM_residue']=""
						
					
                    df_temp=df_temp.sort_values(['Protein-pos','Binomial_final_score','pA_q_value_BA'],ascending=[True,True,False])
                    df_temp['Protein_pos_res']=df_temp['Protein-pos']+"_"+df_temp['PTM_residue']

                    PSM_threshold_2=dict(df_temp.groupby('Protein_pos_res')['0.05FLR_threshold_count'].sum())

                    df_temp=df_temp.drop_duplicates(subset=('Protein_pos_res'),keep="last",inplace=False)

                    #column for all PSM counts at 5%FLR
                    df_temp['Sum_of_PSM_counts(5%FLR)']=df_temp['Protein_pos_res'].map(PSM_threshold_2)

                    df_temp=df_temp[['Peptide_mod_pos', 'pA_q_value_BA','Binomial_final_score', 'Protein_pos_res','Sum_of_PSM_counts(5%FLR)','0.05FLR_threshold_count']]
                    df_temp.rename(columns = {'pA_q_value_BA':dataset+"_FLR",'Binomial_final_score':dataset+'_BinomialScore','Peptide_mod_pos':dataset+"_peptide_mod_pos",'Sum_of_PSM_counts(5%FLR)':dataset+"_Sum_of_PSM_counts(5%FLR)",'0.05FLR_threshold_count':dataset+"_peptidoform_PSMcount(5%FLR)"}, inplace = True)
                    df_temp=df_temp.set_index('Protein_pos_res')


                    if dataset==dataset_list[0]:
                        df_final=df_temp
                    else:
                        df_final=pd.concat([df_final,df_temp],axis=1)

                column_list=[PXD+"_Sum_of_PSM_counts(5%FLR)" for PXD in dataset_list]
                df_final["Sum_of_PSM_counts(5%FLR)"]=df_final[column_list].sum(axis=1)

                cols=df_final.columns.tolist()
                cols=[cols[-1]]+cols[:-1]
                df_final=df_final[cols]
                df_final=df_final.drop(column_list, axis=1)

                df_final = df_final[~df_final.index.str.contains(decoy_prefix)]
                df_final = df_final[~df_final.index.str.contains(contam_prefix)]

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

                df_final.to_csv(output_location + "/G"+str(gold_count)+"S"+str(silver_count)+"B_"+str(flr_filter)+"_protein_pos"+decoy_method+choice+"_allProts.csv", index=False)

                counts=pd.crosstab(df_final.PTM_residue,df_final.PTM_FLR_category).replace(0,np.nan).stack().reset_index().rename(columns={0:'Count'})

                f = lambda x, y, z : np.nan if x!="A" else str(round(z*3/counts.loc[counts['PTM_FLR_category']==y,'Count'].sum()*100,2))+"%"

                counts['FLR']=counts.apply(lambda x: f(x.PTM_residue, x.PTM_FLR_category, x.Count), axis=1)
                counts['FLR']=counts.groupby('PTM_FLR_category')['FLR'].ffill()
                counts['PTM_FLR_category_FLR']=counts['PTM_FLR_category'] + " - " + counts['FLR']
                counts['PTM_FLR_category'] = pd.Categorical(counts['PTM_FLR_category'], ["Bronze", "Silver", "Gold"])
                counts=counts.sort_values("PTM_FLR_category")
                print(choice)
                print(counts)

                fig, axes = plt.subplots(ncols=3)
                for i, (name, group) in enumerate(counts.groupby("PTM_FLR_category_FLR", sort=False)):
                    axes[i].set_title(name)
                    group.plot(kind="bar", x = "PTM_residue", y="Count", ax=axes[i], legend=False)
                    axes[i].set_ylabel("count")
                    axes[i].set_xlabel("")
                plt.tight_layout()
                #plt.show()
                plt.savefig(output_location + "/G"+str(gold_count)+"S"+str(silver_count)+"B_"+str(flr_filter)+"_protein_pos_categories"+decoy_method+choice+"_allProts.png",dpi=300)
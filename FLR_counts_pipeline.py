import pandas as pd
import xml.sax
import numpy as np
import sys
import os

#usage: FLR_counts_pipeline.py filenames.txt
# where filenames.txt is a .txt file listing paths to analysis files 

class TPPHandler(xml.sax.ContentHandler):
    def __init__(self):
        self.Counter=0
    # Call when an element starts
    def startElement(self, tag, attributes):
        if tag == "SpectrumIdentificationResult":
            self.Counter +=1
    def endElement(self, name):
        if name == "MzIdentML":
            #print("END"+ str(self.Counter))
            output.write(str(self.Counter)+",")
            output2.write(str(self.Counter)+",")
            self.Counter=0

parser = xml.sax.make_parser()
Handler = TPPHandler()
parser.setContentHandler(Handler)

print(sys.argv[1])
decoy_list_file = open(sys.argv[1],"r")
decoy_list_file = decoy_list_file.read()
decoy_list = decoy_list_file.replace('\n', '.').split(".") 

output_f="FLRcounts_no_choice_noA_decoy_methods.csv"
output=open(output_f,"w")
output.write("Dataset,No. RAW files,Spectral count,0.01_FDR_PSM_Count,0.01_FDR_PSM_Count_Unique_Mod_Peptide,0.01_Phosphopeptide_Count,"
             "0.01_FDR_Phosphopeptide_Unique_Mod_Peptide_Count,PSM-site count,Peptidoform-site count")
             #"Binomial0.01_FLR_Peptidoform_Site_Count,Binomial0.05_FLR_Peptidoform_Site_Count,Binomial0.1_FLR_Peptidoform_Site_Count,"
             #"Binomial_pA0.01_FLR_Peptidoform_Site_Count,Binomial_pA0.05_FLR_Peptidoform_Site_Count,Binomial_pA0.1_FLR_Peptidoform_Site_Count,"
             #"Peptidoform-site count_without_noChoice,Binomial0.01_FLR_Peptidoform_Site_Count_without_NoChoice,Binomial0.05_FLR_Peptidoform_Site_Count_without_NoChoice,Binomial0.1_FLR_Peptidoform_Site_Count_without_NoChoice,"
             #"Binomial_pA0.01_FLR_Peptidoform_Site_Count_without_NoChoice,Binomial_pA0.05_FLR_Peptidoform_Site_Count_without_NoChoice,Binomial_pA0.1_FLR_Peptidoform_Site_Count_without_NoChoice")
for decoy_method in ["","_peptidoform_decoy","_site_decoy"]:
    file=decoy_list[0]+"/FDR_0.01/"+"binomial_peptidoform_collapsed_FLR"+decoy_method+".csv"
    if os.path.isfile(file):
        output.write(",Binomial0.01_FLR_Peptidoform_Site_Count"+decoy_method+ ",Binomial0.05_FLR_Peptidoform_Site_Count"+decoy_method+ ",Binomial0.1_FLR_Peptidoform_Site_Count"+decoy_method+
                     ",Binomial_pA0.01_FLR_Peptidoform_Site_Count"+decoy_method+ ",Binomial_pA0.05_FLR_Peptidoform_Site_Count"+decoy_method+ ",Binomial_pA0.1_FLR_Peptidoform_Site_Count"+decoy_method)
        if os.path.isfile(file.replace(".csv","_no_choice.csv")):
            output.write(",Binomial_pA0.01_FLR_Peptidoform_Site_Count_without_NoChoice"+decoy_method+ ",FLR0.01 combined probability cutoff_without_NoChoice"+decoy_method+
                     ",Binomial_pA0.05_FLR_Peptidoform_Site_Count_without_NoChoice"+decoy_method+  ",FLR0.05 combined probability cutoff_without_NoChoice"+decoy_method+
                     ",Binomial_pA0.1_FLR_Peptidoform_Site_Count_without_NoChoice"+decoy_method+ ",FLR0.1 combined probability cutoff_without_NoChoice"+decoy_method)
output.write("\n")

output_f2="FLRcounts_noA_decoy_methods.csv"
output2=open(output_f2,"w")
output2.write("Dataset,No. RAW files,Spectral Count,FDR 0.01 PSM Count,FDR 0.01 Phosphopeptide Count,Peptidoform-site count")
             #"FLR pA0.01 Peptidoform-site Count, FLR0.01 combined probability cutoff,"
              #"FLR pA0.05 Peptidoform-site Count, FLR0.05 combined probability cutoff,"
              #"FLR pA0.1 Peptidoform-site Count, FLR0.1 combined probability cutoff")
for decoy_method in ["","_peptidoform_decoy","_site_decoy"]:
    file=decoy_list[0]+"/FDR_0.01/"+"binomial_peptidoform_collapsed_FLR"+decoy_method+".csv"
    if os.path.isfile(file):
        output2.write(",FLR pA0.01 Peptidoform-site Count"+decoy_method+ ",FLR0.01 combined probability cutoff"+decoy_method+
                      ",FLR pA0.05 Peptidoform-site Count"+decoy_method+ ",FLR0.05 combined probability cutoff"+decoy_method+
                      ",FLR pA0.1 Peptidoform-site Count"+decoy_method+ ",FLR0.1 combined probability cutoff"+decoy_method)
output2.write("\n")

for decoy in decoy_list:
    #Dataset
    output.write(decoy+",")
    output2.write(decoy+",")

    file=decoy+"/FDR_0.01/FDR_output.csv"
    df=pd.read_csv(file)

    #no RAW files
    file_count=df['Source'].unique()
    output.write(str(len(file_count))+",")
    output2.write(str(len(file_count))+",")

    xml_file=decoy+"/interact_ipro_ptm.pep.mzid"
    print(xml_file)
    spectral_count=0
    try:
        parser.parse(xml_file)
    except:
        xml_file=xml_file.replace("interact_ipro_ptm.pep.mzid","interact-ipro-ptm.pep.mzid")
        print(xml_file)
        parser.parse(xml_file)

    #Filter for 1%FDR
    df=df.loc[df["q_value"]<=0.01]
    print(decoy + " 1%FDR PSM: "+str(len(df)))
    #PSM.append(len(df))
    #0.01_FDR_PSM_Count
    output.write(str(len(df))+",")
    output2.write(str(len(df))+",")
    df_temp=df.drop_duplicates(subset=["Peptide_mod"])
    #PSM_unique.append(len(df_temp))
    #0.01_FDR_PSM_Count_Unique_Mod_Peptide
    output.write(str(len(df_temp))+",")
    df=df[df['PTM'].notna()]
    df=df.loc[df['PTM'].str.contains("Phospho")]
    print(decoy + " 1%FDR phosphopep: "+str(len(df)))
    #phosphopep.append(len(df))
    #0.01_Phosphopeptide_Count
    output.write(str(len(df))+",")
    output2.write(str(len(df))+",")
    df_temp=df.drop_duplicates(subset=["Peptide_mod"])
    #phosphopep_unique.append(len(df_temp))
    #0.01_FDR_Phosphopeptide_Unique_Mod_Peptide_Count
    output.write(str(len(df_temp))+",")

    file=decoy+"/FDR_0.01/binomial.csv"
    df=pd.read_csv(file)
    print("Binomial PSM: "+str(len(df)))
    #PSM-site count
    output.write(str(len(df))+",")
    counter=0
    for decoy_method in ["","_peptidoform_decoy","_site_decoy"]:
        file=decoy+"/FDR_0.01/"+"binomial_peptidoform_collapsed_FLR"+decoy_method+".csv"
        if os.path.isfile(file):
            counter+=1
            df=pd.read_csv(file)
            print("Peptidoform PSM: "+str(len(df)))
            #Peptidoform-site count
            if counter==1:
                output.write(str(len(df))+",")
                output2.write(str(len(df))+",")

            # remove pA from counts
            r = lambda x, y: 0 if x[int(y) - 1] != "A" else 1
            df['Decoy Modification Site'] = df.apply(lambda x: r(x.Peptide_pos.split("-")[0], x.Peptide_pos.split("-")[1]), axis=1)
            df=df.loc[df["Decoy Modification Site"]!=1]

            print(decoy+" 0.01: "+ str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.01])))
            #Binomial0.01_FLR_Peptidoform_Site_Count
            output.write(str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.01]))+",")
            print(decoy+" 0.05: "+ str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.05])))
            #Binomial0.05_FLR_Peptidoform_Site_Count
            output.write(str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.05]))+",")
            print(decoy+" 0.10: "+ str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.10])))
            #Binomial0.1_FLR_Peptidoform_Site_Count
            output.write(str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.10]))+",")

            print(decoy+" 0.01: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.01])))
            #Binomial_pA0.01_FLR_Peptidoform_Site_Count
            output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.01]))+",")
            output2.write(str(len(df.loc[df["pA_q_value_BA"]<=0.01]))+",")
            #FLR0.01 combined probability cutoff
            output2.write(str(df.loc[df["pA_q_value_BA"]<=0.01]['PTM_final_prob'].min())+",")
            print(decoy+" 0.05: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.05])))
            #Binomial_pA0.05_FLR_Peptidoform_Site_Count
            output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.05]))+",")
            output2.write(str(len(df.loc[df["pA_q_value_BA"]<=0.05]))+",")
            #FLR0.05 combined probability cutoff
            output2.write(str(df.loc[df["pA_q_value_BA"]<=0.05]['PTM_final_prob'].min())+",")
            print(decoy+" 0.10: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.10])))
            #Binomial_pA0.1_FLR_Peptidoform_Site_Count
            output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.10]))+",")
            output2.write(str(len(df.loc[df["pA_q_value_BA"]<=0.10]))+",")
            #FLR0.1 combined probability cutoff
            output2.write(str(df.loc[df["pA_q_value_BA"]<=0.10]['PTM_final_prob'].min())+",")

            #no-choice
            if os.path.isfile(file.replace(".csv","_no_choice.csv")):
                df=pd.read_csv(file.replace(".csv","_no_choice.csv"))
                print("Without no choice:")
                print(decoy+" 0.01: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.01])))
                #Binomial_pA0.01_FLR_Peptidoform_Site_Count_without_NoChoice
                output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.01]))+",")
                #FLR0.01 combined probability cutoff
                output.write(str(df.loc[df["pA_q_value_BA"]<=0.01]['PTM_final_prob'].min())+",")
                print(decoy+" 0.05: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.05])))
                #Binomial_pA0.05_FLR_Peptidoform_Site_Count_without_NoChoice
                output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.05]))+",")
                #FLR0.05 combined probability cutoff
                output.write(str(df.loc[df["pA_q_value_BA"]<=0.05]['PTM_final_prob'].min())+",")
                print(decoy+" 0.10: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.10])))
                #Binomial_pA0.1_FLR_Peptidoform_Site_Count_without_NoChoice
                output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.10]))+",")
                #FLR0.1 combined probability cutoff
                output.write(str(df.loc[df["pA_q_value_BA"]<=0.1]['PTM_final_prob'].min()))

    output2.write("\n")
    output.write("\n")

output.close()
output2.close()

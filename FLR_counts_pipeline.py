import pandas as pd
import xml.sax
import numpy as np
import sys

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

output_f="FLRcounts_no_choice_noA.csv"
output=open(output_f,"w")
output.write("Dataset,No. RAW files,Spectral count,0.01_FDR_PSM_Count,0.01_FDR_PSM_Count_Unique_Mod_Peptide,0.01_Phosphopeptide_Count,"
             "0.01_FDR_Phosphopeptide_Unique_Mod_Peptide_Count,PSM-site count,Peptidoform-site count,"
             "Binomial0.01_FLR_Peptidoform_Site_Count,Binomial0.05_FLR_Peptidoform_Site_Count,Binomial0.1_FLR_Peptidoform_Site_Count,"
             "Binomial_pA0.01_FLR_Peptidoform_Site_Count,Binomial_pA0.05_FLR_Peptidoform_Site_Count,Binomial_pA0.1_FLR_Peptidoform_Site_Count,"
             "Peptidoform-site count_without_noChoice,Binomial0.01_FLR_Peptidoform_Site_Count_without_NoChoice,Binomial0.05_FLR_Peptidoform_Site_Count_without_NoChoice,Binomial0.1_FLR_Peptidoform_Site_Count_without_NoChoice,"
             "Binomial_pA0.01_FLR_Peptidoform_Site_Count_without_NoChoice,Binomial_pA0.05_FLR_Peptidoform_Site_Count_without_NoChoice,Binomial_pA0.1_FLR_Peptidoform_Site_Count_without_NoChoice"
             "\n")
output_f2="FLRcounts_noA.csv"
output2=open(output_f2,"w")
output2.write("Dataset,No. RAW files,Spectral Count,FDR 0.01 PSM Count,FDR 0.01 Phosphopeptide Count,Peptidoform-site count,"
             "FLR pA0.01 Peptidoform-site Count,FLR pA0.05 Peptidoform-site Count,FLR pA0.1 Peptidoform-site Count"
             "\n")

for decoy in decoy_list:

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
    output.write(str(len(df))+",")
    output2.write(str(len(df))+",")
    df_temp=df.drop_duplicates(subset=["Peptide_mod"])
    #PSM_unique.append(len(df_temp))
    output.write(str(len(df_temp))+",")
    df=df[df['PTM'].notna()]
    df=df.loc[df['PTM'].str.contains("Phospho")]
    print(decoy + " 1%FDR phosphopep: "+str(len(df)))
    #phosphopep.append(len(df))
    output.write(str(len(df))+",")
    output2.write(str(len(df))+",")
    df_temp=df.drop_duplicates(subset=["Peptide_mod"])
    #phosphopep_unique.append(len(df_temp))
    output.write(str(len(df_temp))+",")

    file=decoy+"/FDR_0.01/binomial.csv"
    df=pd.read_csv(file)
    print("Binomial PSM: "+str(len(df)))
    output.write(str(len(df))+",")

    file=decoy+"/FDR_0.01/"+"binomial_peptidoform_collapsed_FLR.csv"
    df=pd.read_csv(file)
    print("Peptidoform PSM: "+str(len(df)))
    output.write(str(len(df))+",")
    output2.write(str(len(df))+",")

    # remove pA from counts
    r = lambda x, y: 0 if x[int(y) - 1] != "A" else 1
    df['Decoy Modification Site'] = df.apply(lambda x: r(x.Peptide_pos.split("-")[0], x.Peptide_pos.split("-")[1]), axis=1)
    df=df.loc[df["Decoy Modification Site"]!=1]

    print(decoy+" 0.01: "+ str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.01])))
    output.write(str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.01]))+",")
    print(decoy+" 0.05: "+ str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.05])))
    output.write(str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.05]))+",")
    print(decoy+" 0.10: "+ str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.10])))
    output.write(str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.10]))+",")

    print(decoy+" 0.01: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.01])))
    output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.01]))+",")
    output2.write(str(len(df.loc[df["pA_q_value_BA"]<=0.01]))+",")
    print(decoy+" 0.05: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.05])))
    output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.05]))+",")
    output2.write(str(len(df.loc[df["pA_q_value_BA"]<=0.05]))+",")
    print(decoy+" 0.10: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.10])))
    output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.10]))+",")
    output2.write(str(len(df.loc[df["pA_q_value_BA"]<=0.10]))+"\n")

    #Remove no choice
    df = df[df.Peptide.str.count('S|T|Y|A')>df.Peptide_mod.str.count("Phospho")]
    print("Without no choice:")

    print("Peptidoform PSM: "+str(len(df)))
    output.write(str(len(df))+",")

    print(decoy+" 0.01: "+ str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.01])))
    output.write(str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.01]))+",")
    print(decoy+" 0.05: "+ str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.05])))
    output.write(str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.05]))+",")
    print(decoy+" 0.10: "+ str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.10])))
    output.write(str(len(df.loc[df["Binomial_final_prob_q_value"]<=0.10]))+",")

    print(decoy+" 0.01: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.01])))
    output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.01]))+",")
    print(decoy+" 0.05: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.05])))
    output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.05]))+",")
    print(decoy+" 0.10: "+ str(len(df.loc[df["pA_q_value_BA"]<=0.10])))
    output.write(str(len(df.loc[df["pA_q_value_BA"]<=0.10]))+"\n")
output.close()
output2.close()

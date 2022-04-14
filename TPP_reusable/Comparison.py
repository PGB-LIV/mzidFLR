#Gives a spectrum by spectrum comparison of all search files specified in "params_TPP"
#Creates a subfolder within "Spectrum_Comparisons" with the search names (specified within "params_TPP") separated by "_" (eg. pA_pG_pL_pD_pE_pP)
#output="FDR_[cutoff]_PTM_score_0_spectrum_comparison_[searchnames].csv" (where each row represents each spectrum and gives the peptide, protein, score, mods et. found in each search for that spectrum),"All_confident_PTM_no_collapse_Site-based_spectrum_match.csv" (binary spectrum match column, do search comparisons provide same spectra results)


import os
import pandas as pd
import csv
import re

#Spectrum comparison between search files
def per_spectrum(wd,file_list,file_name_list,output):
    counter=1
    os.chdir(wd)
    search_spectrum = []
    dict = {}
    header = "Spectrum,"
    for f, n in zip(file_list, file_name_list):
        dict[f]={}
        reader = csv.reader(open(f, 'r'))
        counter=1
        for row in reader:
            #Spectrum,Peptide_mod_TPP,Protein_TPP,USI_TPP,Score_TPP,PTMs_TPP,PTM_scores_all_TPP,PTM_positions_TPP,PTM_Protein_Positions_TPP,PTM_info_TPP,
            counter+=1
            if row[0] == "Peptide_mod":
                list = [0, 2, 10, 4, 5, 7, 6, 13, 8]
                for i in list:
                    header += row[i] + "_" + n + ","
            else:
                spectrum = row[9].replace(".mzML", "").replace(".0", ".")
                dict[f][spectrum] = row
                if spectrum not in search_spectrum:
                    search_spectrum.append(spectrum)
    output_file = open(output, 'w')
    output_file.write(header + "\n")
    output_file.close()

    output_file = open(output, 'a')
    for s in search_spectrum:
        line = s + ","
        for f in file_list:
            if s in dict[f]:
                temp = dict[f][s]
                temp_all = temp[0] + "," + temp[2] + "," + temp[4] + "," + temp[6] + "," + temp[7] + "," + temp[8] + "," + temp[
                    9] + "," + temp[10] + "," + temp[11] + ","
                line += temp_all
            else:
                line += ',' + ',' + ',' + ',' + ',' + ',' + ',' + ',' + ','
        output_file.write(line + "\n")
        counter += 1
    output_file.close()

#Create csv of matching spectrum 
def matches(input, file_list, final_output):
    df = pd.read_csv(input, sep=",")
    for i in file_list:
        df['PTM_'+i] = df['PTM_'+i].astype(str)
        df['Peptide_mod_' + i] = df['Peptide_mod_' + i].astype(str)
    matches = []
    pep_matches = []
    mass_matches = []
    for a in range(len(df)):
        peptides=[]
        pep_unmod=[]
        pep_mass=[]
        mass=0
        for i in file_list:
            peptides.append(df.loc[a,'Peptide_mod_'+i])
            pep_temp= re.sub("[\[].*?[\]]", "",df.loc[a,'Peptide_mod_'+i])
            pep_unmod.append(pep_temp)
            if df.loc[a,'PTM_'+i]=="nan":
                PTM_list=""
            else:
                PTM_list=df.loc[a,'PTM_'+i]
            for p in PTM_list.split(";"):
                if "Phosphorylation" in p:
                    mass += 80
                if "Pyrophosphorylation" in p:
                    mass += 160
                if "Oxidation" in p:
                    mass += 16
            pep_mass.append(mass)
        if len(set(peptides)) == 1:
            matches.append("TRUE")
        else:
            matches.append("FALSE")
        if len(set(pep_unmod)) == 1:
            pep_matches.append("TRUE")
        else:
            pep_matches.append("FALSE")
        if len(set(pep_mass)) == 1:
            mass_matches.append("TRUE")
        else:
            mass_matches.append("FALSE")
    df['Match'] = matches
    df['Pep_Match'] = pep_matches
    df['Mass_shift_Match'] = mass_matches
    df.to_csv(final_output, index=False)


def comparison(working,sub,file_list,file_name_list):
    os.chdir(working)    
    if not os.path.exists("Spectrum_Comparisons"):
        os.mkdir("Spectrum_Comparisons")
    file_list_output=""
    for i in file_name_list:
        file_list_output+=i+"_"
    spectrum_file="Spectrum_Comparisons/"+file_list_output[:-1]
    if not os.path.exists(spectrum_file):
        os.mkdir(spectrum_file)
    output=working+"/"+spectrum_file+"/"+sub+"_spectrum_comparison_"+file_list_output[:-1]+".txt"
    per_spectrum(working,file_list,file_name_list,output)
    #Seperate matches
    final=working+"/"+spectrum_file+"/"+sub+"_spectrum_comparison_"+file_list_output[:-1]+".csv"
    matches(output,file_name_list,final)
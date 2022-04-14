import pandas as pd
import os
import time
import glob

import TPP_reusable.FDR as FDR
import TPP_reusable.Unique_PTMs as Unique_PTMs
import TPP_reusable.convert_mzIdentML_sax as convert_mzIdentML_sax
import TPP_reusable.Comparison as Comparison
import TPP_reusable.Post_analysis as Post_analysis
import TPP_reusable.pAla as pAla
import TPP_reusable.Binomial_adjustment as Binomial_adjustment														  
from params_TPP import *

start_time = time.time()

PTM_score_cutoff = 0
sub = "FDR_" + str(FDR_cutoff) + "_PTM_score_" + str(PTM_score_cutoff)
decoy_list = []
for wd, software in zip(w, s):
    try:
        os.chdir(wd)
    except OSError:
        print("Can't find results directory, please check location: ")
    print(wd)
    # Run pepXML converter - if csv not in dir
    results_file = "interact-ipro-ptm.pep.mzid.csv"
    if os.path.isfile(results_file):
        print("XML converted file exists")
    else:
        xml = glob.glob("*.mzid")
        try:
            print("Converting XML file")
            convert_mzIdentML_sax.convert(xml[0])
            print("XML converted")
        except ValueError:
            raise ValueError("Please provide TPP '.mzid' results file")
    # calculate FDR
    if not os.path.exists(sub):
        os.mkdir(sub)
    FDR_output = sub + "/FDR_output.csv"
    if "pL" in software or "pLeu" in software or "pL" in wd or "pLeu" in wd:
        search = "STYL"
        decoy = "pLeu"
    elif "pG" in software or "pGly" in software or "pG" in wd or "pGly" in wd:
        search = "STYG"
        decoy = "pGly"
    elif "pD" in software or "pAsp" in software or "pD" in wd or "pAsp" in wd:
        search = "STYD"
        decoy = "pAsp"
    elif "pE" in software or "pGlu" in software or "pE" in wd or "pGlu" in wd:
        search = "STYE"
        decoy = "pGlu"
    elif "pP" in software or "pPro" in software or "pP" in wd or "pPro" in wd:
        search = "STYP"
        decoy = "pPro"
    else:
        search = "STYA"
        decoy = "pAla"
    decoy_list.append(decoy)
    search_software = "TPP"
    FDR.calculateFDR(results_file, FDR_output, PXD)
    FDR.ppm_error(FDR_output)
    print("FDR calculation done")
    print("--- %s seconds ---" % (time.time() - start_time))
    #Filter for FDR
    # Collapse for best scoring for each peptide, protein site, mass shift on protein or no collapse
    unique_peptide = sub + "/Peptide_confident_PTM_unique.csv"
    unique_site = sub + "/Site_confident_PTM_unique.csv"
    unique_mass = sub + "/Peptide_mass_confident_PTM_unique.csv"
    non_collapse = sub + "/All_confident_PTM_no_collapse.csv"
    Unique_PTMs.unique(FDR_output,FDR_cutoff, unique_peptide, unique_site, unique_mass, non_collapse)
    print("Unique PTMs done")
    print(FDR_cutoff, software, "--- %s seconds ---" % (time.time() - start_time))

# Compare spectrum between different searches
print("Spectrum Comparison --- %s seconds ---" % (time.time() - start_time))
files = []
for i in w:
    files.append(i + "/" + sub + "/All_confident_PTM_no_collapse.csv")
Comparison.comparison(p, sub, files, s)
print("Spectrum Comparison Done --- %s seconds ---" % (time.time() - start_time))

file_input_list = ["All_confident_PTM_no_collapse.csv"]
files = []
for i in w:
    files.append(i + "/" + sub + "/")
software_list_edit = ""
for i in s:
    software_list_edit += i + "_"
os.chdir(p)
if not os.path.exists("Comparisons"):
    os.mkdir("Comparisons")
# load in list of hits where spectrum comparisons gave match between all searches
spectrum_match_list = []
spectrum = p + "/Spectrum_Comparisons/" + software_list_edit[
                                          :-1] + "/" + sub + "_spectrum_comparison_" + software_list_edit[
                                                                                       :-1] + ".csv"
df = pd.read_csv(spectrum, dtype=str)
for i in range(len(df)):
    if df.loc[i, 'Match'] == "TRUE":
        spectrum_match_list.append(df.loc[i, 'Spectrum'])

print("Starting FLR calculations --- %s seconds ---" % (time.time() - start_time))
# Post analysis - FLR calulations and plots
files = []
for i in w:
    for f in file_input_list:
        files.append(i + "/" + sub + "/" + f)

        # print(i + "/" + sub)
for file in files:
    Post_analysis.site_input = Post_analysis.site_based(file)
    Post_analysis.site_all_input = Post_analysis.site_based_all(file)
    output = file.replace(".csv", "_Site-based.csv")
    Post_analysis.spectrum_comparison(output, spectrum_match_list)
final_file_list = []
final_file_list_collapsed = []
final_file_list_filtered = []
final_file_list_PTM = []
for i, decoy in zip(w, decoy_list):
    Post_analysis.model_FLR(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match.csv")
    Post_analysis.model_FLR_new(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match.csv")
    Post_analysis.model_FLR_filter(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match.csv")
    Post_analysis.PTM_sort_FLR(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match.csv")
    pAla.calulate_decoy_FLR(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR.csv",decoy)
    pAla.calulate_decoy_FLR(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match_new_FLR_collapse.csv", decoy)
    pAla.calulate_decoy_FLR(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_filtered.csv", decoy)
    pAla.calulate_decoy_FLR(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_PTM_sort.csv", decoy)
    final_file_list.append(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_" + decoy + ".csv")
    final_file_list_collapsed.append(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match_new_FLR_collapse_" + decoy + ".csv")
    final_file_list_filtered.append(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_filtered_" + decoy + ".csv")
    final_file_list_PTM.append(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_PTM_sort_" + decoy + ".csv")

Post_analysis.plot_FLR_comparisons(sub, final_file_list, s, p)
Post_analysis.plot_FLR_comparisons(sub, final_file_list_collapsed, s, p)
Post_analysis.plot_FLR_comparisons_PTM_prob(sub, final_file_list_PTM, s, p)
Post_analysis.plot_FLR_comparisons(sub, final_file_list_filtered, s, p)

final_file_list = []
final_file_list_collapse = []
for i, decoy in zip(w, decoy_list):
    Post_analysis.model_FLR(i + "/" + sub + "/" + "Site_confident_PTM_unique.csv")
    Post_analysis.model_FLR_new(i + "/" + sub + "/" + "Site_confident_PTM_unique.csv")
    pAla.calulate_decoy_FLR(i + "/" + sub + "/" + "Site_confident_PTM_unique_FLR.csv", decoy)
    pAla.calulate_decoy_FLR(i + "/" + sub + "/" + "Site_confident_PTM_unique_new_FLR_collapse.csv", decoy)
    final_file_list.append(i + "/" + sub + "/" + "Site_confident_PTM_unique_FLR_" + decoy + ".csv")
    final_file_list_collapse.append(i + "/" + sub + "/" + "Site_confident_PTM_unique_new_FLR_collapse_" + decoy + ".csv")
Post_analysis.plot_FLR_comparisons(sub, final_file_list, s, p)
Post_analysis.plot_FLR_comparisons(sub, final_file_list_collapse, s, p)

for i, decoy in zip(w, decoy_list):
   Binomial_adjustment.Binomial(i + "/" + sub + "/" + "All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_" + decoy + ".csv",decoy)
binomial_list=[]
binomial_name_list=[]
for a,b in zip(w,s):
    collapsed_output = a + "/" + sub + "/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla_FLR_collapse_pep_pos_final_score.csv"
    binomial_list.append(collapsed_output)
    binomial_name_list.append(b)
    FLR_output = a+ "/" + sub  + "/binomial_collapsed_FLR.csv"
    binomial_list.append(FLR_output)
    binomial_name_list.append(b)

Binomial_adjustment.plot_FLR_comparisons(binomial_list,binomial_name_list,p)
print("FLR calculations done --- %s seconds ---" % (time.time() - start_time))

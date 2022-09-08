import os
import time
import sys

import TPP_reusable.FDR as FDR
import TPP_reusable.convert_mzIdentML_sax as convert_mzIdentML_sax
import TPP_reusable.Post_analysis as Post_analysis
import TPP_reusable.pAla as pAla
import TPP_reusable.Binomial_adjustment as Binomial_adjustment

start_time = time.time()

#TPP_comparison.py [mzid_file] [PXD] [optional: FDR_cutoff] [optional: PTM_score_cutoff]
if len(sys.argv)<2:
    sys.exit("Provide mzid file and PXD identifier. TPP_comparison.py [mzid_file] [PXD] [optional: FDR_cutoff] [optional: PTM_score_cutoff]")
mzid_file=sys.argv[1]
if mzid_file[-5:]!=".mzid":
    print(mzid_file[-5:])
    sys.exit("Provide mzid file. TPP_comparison.py [mzid_file] [PXD] [optional: FDR_cutoff] [optional: PTM_score_cutoff]")
PXD = sys.argv[2]
if "PXD" not in PXD:
    sys.exit("Provide PXD identifier. TPP_comparison.py [mzid_file] [PXD] [optional: FDR_cutoff] [optional: PTM_score_cutoff]")
try:
    FDR_cutoff = sys.argv[3]
except:
    FDR_cutoff=0.01
    print("FDR cutoff not specified, 0.01 default FDR cutoff used")
try:
    PTM_score_cutoff = sys.arg[4]
except:
    PTM_score_cutoff = 0
    print("PTM score cutoff not specified, 0 default PTM score cutoff used")

sub = "FDR_" + str(FDR_cutoff) + "_PTM_score_" + str(PTM_score_cutoff)
decoy_list = []

results_file = "interact-ipro-ptm.pep.mzid.csv"

if os.path.isfile(results_file):
    print("mzid converted file exists")
else:
    try:
        print("Converting mzid file")
        convert_mzIdentML_sax.convert(mzid_file)
        print("mzid converted")
    except ValueError:
        raise ValueError("Please provide TPP '.mzid' results file")
# calculate FDR
if not os.path.exists(sub):
    os.mkdir(sub)
FDR_output = sub + "/FDR_output.csv"
search = "STYA"
decoy = "pAla"
decoy_list.append(decoy)
search_software = "TPP"
print("Starting FDR calculations")
print("--- %s seconds ---" % (time.time() - start_time))
FDR.calculateFDR(results_file, FDR_output, PXD)
FDR.ppm_error(FDR_output)
print("FDR calculation done")
print("Starting FLR calculations")
print("--- %s seconds ---" % (time.time() - start_time))
# Post analysis - FLR calulations and plots
Post_analysis.site_input = Post_analysis.site_based(FDR_output,FDR_cutoff)
Post_analysis.model_FLR(sub + "/" + "Site-based.csv")
pAla.calulate_decoy_FLR(sub + "/" + "Site-based_FLR.csv",decoy)
Binomial_adjustment.Binomial(sub + "/" + "Site-based_FLR_" + decoy + ".csv",decoy)
print("FLR calculations done --- %s seconds ---" % (time.time() - start_time))

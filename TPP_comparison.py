import os
import time
import sys
import glob

import TPP_reusable.FDR as FDR
import TPP_reusable.convert_mzIdentML_sax as convert_mzIdentML_sax
import TPP_reusable.convert_mzIdentML_sax_MSFrag as convert_mzIdentML_sax_MSFrag
import TPP_reusable.Post_analysis as Post_analysis
import TPP_reusable.pAla as pAla
import TPP_reusable.Binomial_adjustment as Binomial_adjustment

start_time = time.time()

#TPP_comparison.py [mzid_file] [PXD] [mod] [optional: PSM FDR_cutoff]
if len(sys.argv)<3 or len(sys.argv)>5:
    sys.exit("Provide mzid file, PXD identifier and modification of interest. TPP_comparison.py [mzid_file] [PXD] [modification:target:decoy] [optional: PSM_FDR_cutoff]")
mzid_file=sys.argv[1]
if mzid_file[-5:]!=".mzid":
    print(mzid_file[-5:])
    sys.exit("Provide mzid file. TPP_comparison.py [mzid_file] [PXD] [modification:target:decoy]] [optional: PSM_FDR_cutoff]")
PXD = sys.argv[2]
if "PXD" not in PXD:
    sys.exit("Provide PXD identifier. TPP_comparison.py [mzid_file] [modification:target:decoy]] [PXD] [optional: PSM_FDR_cutoff]")
mod_info = sys.argv[3]
if len(mod_info.split(":"))!=3:
	sys.exit("Provide modification in format modification:target:decoy. E.g. Phospho:STY:A. TPP_comparison.py [mzid_file] [modification:target:decoy]] [PXD] [optional: FPSM_DR_cutoff]")
try:
    FDR_cutoff = sys.argv[4]
except:
    FDR_cutoff=0.01
    print("PSM FDR cutoff not specified, 0.01 default PSM FDR cutoff used")

sub = "FDR_" + str(FDR_cutoff)
decoy_list = []

results_file = mzid_file+".csv"

if os.path.isfile(results_file):
        print("XML converted file exists")
else:
	xml = glob.glob("*.mzid")
	try:
		print("Converting XML file")
		convert_mzIdentML_sax.convert(xml[0])
	except ValueError:
		raise ValueError("Please provide TPP '.mzid' results file")
	if os.path.isfile(results_file):
		print("XML converted")
	else:
		try:
			convert_mzIdentML_sax_MSFrag.convert(xml[0])
			print("XML converted")
		except ValueError:
			raise ValueError("Please provide TPP '.mzid' results file")
# calculate FDR
if not os.path.exists(sub):
	os.mkdir(sub)
FDR_output = sub + "/FDR_output.csv"
mod=mod_info.split(":")[0]
targets=mod_info.split(":")[1]
decoy=mod_info.split(":")[2]
print("Starting FDR calculations")
print("--- %s seconds ---" % (time.time() - start_time))
FDR.calculateFDR(results_file, FDR_output, PXD, mod)
FDR.ppm_error(FDR_output)
print("FDR calculation done")
print("Starting FLR calculations")
print("--- %s seconds ---" % (time.time() - start_time))
# Post analysis - FLR calulations and plots
Post_analysis.site_input = Post_analysis.site_based(FDR_output,FDR_cutoff,mod)
Post_analysis.model_FLR(sub + "/" + "Site-based.csv",mod)
pAla.calulate_decoy_FLR(sub + "/" + "Site-based_FLR.csv",decoy,targets)
Binomial_adjustment.Binomial(sub + "/" + "Site-based_FLR_p" + decoy + ".csv",decoy, targets)
print("FLR calculations done --- %s seconds ---" % (time.time() - start_time))

os.remove(sub+"/Site-based.csv")
os.remove(sub+"/Site-based_FLR.csv")
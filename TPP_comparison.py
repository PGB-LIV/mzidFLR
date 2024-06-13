import os
import time
import sys
import pandas as pd

import TPP_reusable.FDR as FDR
import TPP_reusable.convert_mzIdentML_sax as convert_mzIdentML_sax
import TPP_reusable.convert_mzIdentML_scores as convert_mzIdentML_sax_scores
import TPP_reusable.convert_mzIdentML_sax_MSFrag as convert_mzIdentML_sax_MSFrag
import TPP_reusable.Post_analysis as Post_analysis
import TPP_reusable.Binomial_adjustment as Binomial_adjustment

start_time = time.time()

#--verbose function
opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
if len(opts)==1 and opts[0]=="--verbose":
	print("Using verbose mode")
	verbose=True
else:
	verbose=False
#everything else args, handled below
args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

mod_id="NA"
mod_mass="NA"
#TPP_comparison.py [mzid_file] [PXD] [modification:target:decoy] [optional: decoy prefix] [optional: modification:mass(2dp)] [optional: PSM FDR_cutoff] [optional: Decoy method]
if len(args)<3 or len(args)>7:
	sys.exit("Provide mzid file, PXD identifier and modification of interest. TPP_comparison.py [mzid_file] [PXD] [modification:target:decoy] [optional: decoy prefix]  [optional: modification:mass(2dp)] [optional: PSM FDR_cutoff] [optional: Decoy method 'peptidoform' or 'site']")
mzid_file=args[0]
if mzid_file[-5:]!=".mzid":
    print(mzid_file[-5:])
    sys.exit("Provide mzid file. TPP_comparison.py [mzid_file] [PXD] [modification:target:decoy] [optional: decoy prefix]  [optional: modification:mass(2dp)] [optional: PSM FDR_cutoff] [optional: Decoy method 'peptidoform' or 'site']")
PXD = args[1]
if "PXD" not in PXD:
    sys.exit("Provide PXD identifier. TPP_comparison.py [mzid_file] [PXD] [modification:target:decoy] [optional: decoy prefix]  [optional: modification:mass(2dp)] [optional: PSM FDR_cutoff] [optional: Decoy method 'peptidoform' or 'site']")
mod_info = args[2]
if len(mod_info.split(":"))!=3:
	sys.exit("Provide modification in format modification:target:decoy. E.g. Phospho:STY:A. TPP_comparison.py [mzid_file] [PXD] [modification:target:decoy] [optional: decoy prefix]  [optional: modification:mass(2dp)] [optional: PSM FDR_cutoff] [optional: Decoy method 'peptidoform' or 'site']")

if len(args)>3:
	for i in range(3,len(args)):
		temp = args[i]
		try:
			float(temp)
			FDR_cutoff=temp
		except:
			if (":") in temp:
				mod_id=temp.split(":")[0]
				mod_mass=temp.split(":")[1]
			else:
				if temp.lower()=="peptidoform" or temp.lower()=="site":
					decoy_method = temp.lower()
				else:
					decoy_prefix=temp
try:
	print("Using decoy prefix: " + decoy_prefix)
except:
	print("Decoy prefix not specified, \"DECOY\" default prefix used")
	decoy_prefix="DECOY"
try:
	print("Using FDR cutoff: " + FDR_cutoff)
except:
	print("PSM FDR cutoff not specified, 0.01 default PSM FDR cutoff used")
	FDR_cutoff=0.01
try:
	print("Using decoy mode: "+decoy_method)
except:
	print("Decoy method not specified, using \"Peptidoform\" decoy method")
	decoy_method="peptidoform"

if mod_id!="NA" and mod_mass!="NA":
	print("Modification of interest  "+mod_id+":"+mod_mass)
sub = "FDR_" + str(FDR_cutoff)
decoy_list = []

results_file = mzid_file+".csv"

if os.path.isfile(results_file):
        print("mzid converted file exists")
else:
	try:
		print("Converting mzid file")
		if verbose==True:
			convert_mzIdentML_sax_scores.convert(mzid_file)
		else:
			convert_mzIdentML_sax.convert(mzid_file)
	except ValueError:
		raise ValueError("Please provide TPP '.mzid' results file")
	if os.path.isfile(results_file):
		print("mzid converted")
	else:
		try:
			convert_mzIdentML_sax_MSFrag.convert(mzid_file)
			print("mzid converted")
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
FDR.calculateFDR(results_file, FDR_output, PXD, decoy_prefix, mod, mod_id, mod_mass, verbose)
FDR.ppm_error(FDR_output)
print("FDR calculation done")
print("Starting FLR calculations")
print("--- %s seconds ---" % (time.time() - start_time))
# Post analysis - FLR calulations and plots
Post_analysis.site_input = Post_analysis.site_based(FDR_output,FDR_cutoff,mod,verbose,decoy_prefix)
Post_analysis.model_FLR(sub + "/Site-based.csv",mod,verbose, decoy_prefix)
Post_analysis.calculate_decoy_FLR(sub + "/Site-based_FLR.csv",decoy,targets, verbose,decoy_method)
Binomial_adjustment.Binomial(sub + "/Site-based_FLR.csv",decoy, targets, verbose, decoy_method)

FLR_output = sub + "/binomial_peptidoform_collapsed_FLR.csv"
if verbose:
	FLR_output=FLR_output.replace(".csv","_verbose.csv")
elif decoy_method=="site":
	FLR_output=FLR_output.replace(".csv","_site_decoy.csv")

#Peptidoform to peptide
Post_analysis.peptidoform_to_peptide(FLR_output,mod, verbose)
Binomial_adjustment.calulate_decoy_FLR(FLR_output,decoy,targets,verbose, decoy_method)

df = pd.read_csv(FLR_output)
#Remove no choice - only BA pA FLR
df = df[df.Peptide.str.count('S|T|Y|A')>df.Peptide_mod.str.count("Phospho")]
df.to_csv(FLR_output.replace(".csv","_no_choice.csv"),index=False)
Binomial_adjustment.calulate_decoy_FLR(FLR_output.replace(".csv","_no_choice.csv"),decoy,targets,verbose, decoy_method)

print("Workflow complete --- %s seconds ---" % (time.time() - start_time))

try:
	os.remove(sub+"/Site-based.csv")
except:
	os.remove(sub+"/Site-based_verbose.csv")

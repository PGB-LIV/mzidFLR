# mzidFLR pipeline
 
Scripts for estimating false localisation rates (FLR) using two methods, model FLR and Decoy amino acid FLR, starting from mzid files (pep.xml pipeline [here](https://github.com/PGB-LIV/PhosphoFLR)). 

# Usage


*TPP_reusable* folder contains supported module code for the analysis of Datasets using TPP PTMprophet mzid output files.

      $py TPP_comparison.py [mzid_file] [PXD] [mod:target:decoy] [optional: decoy prefix]  [optional: modification:mass(2dp)] [optional: FDR_cutoff]

To run in verbose mode, include "--verbose" as command line parameter


![Workflow_image](https://user-images.githubusercontent.com/57440286/205335117-e3eea3e7-371c-4736-9d7a-2baf0f10996f.jpg)

Seperate script (Site_based.py) generates site based files:

      $py Site_based.py [optional: PMID] [optional: SDRF file]

Generates "Site-PSM-centric.tsv" and "Site-Peptidofom-centric.tsv" 


# Pipeline description

	Main script "TPP_comparison.py" calls modules:

	Convert_mzIdentML_sax(_MSFrag).py
		Converts mzid file to csv for downstream analysis
		Input = *.mzid
		Output = *.mzid.csv
	FDR.py
		Calculate PSM level FDR (decoy/target count)
		- FDR.jpg – PSM count vs PSM level FDR (with q-value trend)
		- FDR_score.jpg – Peptide probability vs global FDR
		Calculate mass tolerance
		- PPM_error_FDR.jpg – PSM probability vs PPM error
		Input=*.mzid.csv
		Output=FDR_output.csv
	Post_analysis.py
		Filter for PSM level FDR (default = 0.01 or specified threshold)
		Expand to site based format – one site per row
		Filter for specified modification (based on mzid MS name, if unknown mod use optional mod mass mapping [optional: modification:mass(2dp)])
		Remove decoy and contaminant hits
		Calculate model FLR
		- 1-(PTMscore*PSMscore) / SiteCount
		Calculate decoy FLR
		- Ratio* DecoyCount / (SiteCount-DecoyCount)
			Ratio = TargetCount/DecoyCount (ie.STY/A) 
			DecoyCount=count of sites where peptidoform contains decoy modification, regardless of site
		Input=FDR_output.csv
		Output=Site-based_FLR.csv
	Binomial_adjustment.py
		Binomial correction 
		- P_x=(n¦x) p^x q^(n-x)
			P=sum of unique PSMs containing decoy modification/sum of unique PSMs
			X=modified site count
			N=count site seen, modified or not
		Collapse by peptidoform site
		Recalculate model FLR 
		Recalculate decoy FLR
		Input=Site-based_FLR.csv
		Output=binomial_peptidoform_collapsed_FLR.csv



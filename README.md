# mzidFLR pipeline
 
Scripts for estimating false localisation rates (FLR) using two methods, model FLR and Decoy amino acid FLR, starting from mzid files (pep.xml pipeline [here](https://github.com/PGB-LIV/PhosphoFLR)). 

# Usage


*TPP_reusable* folder contains supported module code for the analysis of Datasets using TPP PTMprophet mzid output files.

      $py TPP_comparison.py [mzid_file] [PXD] [modification:target:decoy] [optional: decoy prefix] [optional: modification:mass(2dp)] [optional: PSM FDR_cutoff] [optional: Decoy method] [optional: target species eg species_YEAST] [optional: contaminant prefix eg contam_CONTAM or if prefix is absent and needs to be added, contam_UNKNOWN]


To run in verbose mode, include "--verbose" as command line parameter

*[optional: decoy prefix]* -
As default, uses "DECOY_" as database decoy prefix. Can be changed using this parameter

*[optional: modification:mass(2dp)]* - 
Optional parameter to include modification name for specified mod mass (given to 2dp)

*[optional: FDR_cutoff]* - 
As default, uses an FDR cutoff of 0.01 but can be changed using optional parameter

*[optional: Decoy method 'peptidoform' or 'site']* -
As default, uses "peptidoform decoy method" where decoys are considered any site where there is a pA decoy on the peptidoform. Decoy method optional parameter can be changed to "site" in order to consider decoys as just those where the site is a decoy amino acid, disregarding other sites in the peptidoform. 

*[optional: target species prefix]* -
Target species can be specified using this parameter. If contaminant prefix is given as "UNKNOWN" (ie. prefix is absent and needs to be added), then this MUST be specified, to allow removal of target species from contaminants list when assigning contaminant proteins. 

*[optional: contaminant prefix]* -
As default, uses "CONTAM_" as database contaminant prefix. Can be changed using this parameter. 
If contaminant prefix absent from database, use "contam_UNKNOWN"


![Workflow_image](https://user-images.githubusercontent.com/57440286/205335117-e3eea3e7-371c-4736-9d7a-2baf0f10996f.jpg)


## FLR_counts_pipeline.py
Generates FLR counts for all searches:

	$py FLR_counts_pipeline.p [file_list.txt]
 
Where file_list.txt contains the locations of analysis files (ie. PXD/Experiment_name).

Generates "FLRcounts_noA.csv" and "FLRcounts_no_choice_noA.csv" where no-choice peptides are ignored from FLR counts.


## Site_based_format_GSB_counts_pipeline.py

Generates site based files and Gold-Silver-Bronze classification:

	$py Site_based_format_GSB_counts_pipeline.py [file_list.txt] [Gold count threshold] [Silver count threshold] [optional: meta.csv] [optional: SDRF location]

Where file_list.txt contains the locations of analysis files (ie. PXD/Experiment_name), meta.csv contains available metadata (at least PMID) and SDRF location is the location of SDRF files (ie. Meta/SDRFs)

Generates "Site-PSM-centric.tsv" and "Site-Peptidofom-centric.tsv" per experiment as well as PXD merged files.
	"GSB_0.05_protein_pos.csv" and "GSB_0.05_protein_pos_categories.png" show allocations and counts of Gold-Silver-Bronze categories




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
  		If contaminant prefix is "UNKNOWN", contaminant prefix is added, removing target species from contaminant list
		Remove decoy and contaminant hits, where all proteins are decoy or contaminant
		Calculate model FLR
		- 1-(PTMscore*PSMscore) / SiteCount
		Calculate decoy FLR
		- Ratio* DecoyCount / (SiteCount-DecoyCount)
			Ratio = TargetCount/DecoyCount (ie.STY/A) 
			DecoyCount=count of sites where peptidoform contains decoy modification, regardless of site (using Peptidoform decoy method)
   				alternatively, DecoyCount = count of sites with decoy modification, regardless of peptidoform (using Site decoy method)
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



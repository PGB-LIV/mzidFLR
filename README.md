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
Using this paramter requires a list of contaminant IDs to be provided. We have supplied the file "cRAP_contaminants.txt" which contains the IDs of cRAP database contaminants. Any further contaminants used would need to be added to this list in order to be correctly identified.


![Workflow_image](https://user-images.githubusercontent.com/57440286/205335117-e3eea3e7-371c-4736-9d7a-2baf0f10996f.jpg)


## FLR_counts_pipeline.py
Generates FLR counts for all searches:

	$py FLR_counts_pipeline.py [file_list.txt] [optional:decoy amino acid]
 
Where file_list.txt contains the locations of analysis files (ie. PXD/Experiment_name). Decoy amino acid should be specified as single letter code ie.A.
 eg. $py FLR_counts_pipeline.py filenames.txt A 

Generates "FLRcounts_noA.csv" and "FLRcounts_no_choice_noA.csv" where no-choice peptides are ignored from FLR counts.


## Site_based_format_GSB_counts_pipeline.py

Generates site based files and Gold-Silver-Bronze classification:

Can be ran in "SDRF mode", where SDRF files are present giving meta-data per sample:

	$py Site_based_format_GSB_counts_pipeline.py [file_list.txt] [meta.tsv] [SDRF location] NA [Gold count threshold] [Silver count threshold] [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]
 
 OR "Simple meta-data mode", where tsv files are given containing metadata per PXD in the given location "[simple meta location]":

 	$py Site_based_format_GSB_counts_pipeline.py [file_list.txt] NA NA [simple meta location] [Gold count threshold] [Silver count threshold] [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]

OR "No meta-data mode", where no SDRF or meta-data is available, all meta-data columns will be populated with "NA":

	$py Site_based_format_GSB_counts_pipeline.py [file_list.txt] NA NA NA [Gold count threshold] [Silver count threshold] [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]

Requires txt file with locations of analysis files (ie. PXD/Experiment_name), csv file (meta.csv) with meta data (at least PMID), sdrf file location (eg. "PRIDE/SDRFs/"), gold count and silver count - if meta or sdrf files not available, use "NA"
Alternatively, can use simple meta data files, please give location (eg. "simple_meta/") when the meta data is the same for all files in a data set
Can also accept optional decoy and contam prefix as well as modification:target:decoy, if not specified "DECOY" and "CONTAM" prefixes will be used as default 
and phospho:STY:A used as search mod.

Simple meta-data tsv files MUST contain the following column names: 
- PubMedIDs
- Sample ID
- Organism
- Organism Part
- Cell Line
- Disease Information
 

Generates the following files: 
- CSV file and png showing allocations and counts of Gold-Silver-Bronze categories, for single and all protein mapping for each decoy-method used in the pipeline. 
- CSV giving the residue counts in each category for single and all protein mapping. 
- TSV file for PSM and Peptidoform centiric formats; per experiment, merged per PXD dataset and merged all datasets. 
- Uniprot format PSM/Peptidoform centric merged all datasets


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



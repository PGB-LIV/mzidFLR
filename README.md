# mzid_pipeline
 
Scripts for estimating false localisation rates (FLR) using two methods, model FLR and Decoy amino acid FLR, starting from mzid files (pep.xml pipeline [here](https://github.com/PGB-LIV/PhosphoFLR)). 

Also generates peptidoform site based format. 

# Usage

*TPP_reusable* folder contains supported module code for the analysis of Datasets using TPP PTMprophet mzid output files.

      $py python TPP_comparison.py [mzid_file] [PXD] [mod:target:decoy] [optional: FDR_cutoff]


****Site_peptidoform_centric_format_no_pApeptidoforms.py****

Used for creating site based peptidoform format.
 - Takes _binomial_collpased_FLR.csv_ as input

# Pipeline description

Main script "TPP_comparison.py" calls modules:

	"convert_mzIdentML_sax.py" -> converts MZML to CSV
	"FDR.py" -> calculate Global FDR
		output: "FDR_output.csv"
	"Post_analysis.py" 
			-filter for FDR (default 0.01 or user specified)
			-removes decoy and contam hits
			-expands to site-based format (one row shows PSM with PTM score/position for one site on each **peptide**)
			-only keep phospho sites
		calls "pAla.py" -> calculates decoy amino acid FLR
		output: "Site-based_FLR_pAla.csv"
	"Binomial_adjusment.py" -> calculates pA FLR using Oscar's binomial method. 
		outputs:
		1) "binomial.csv"
				-binomial FLR calculated
				-all hits
		2) "binomial_collapsed_FLR.csv"
				-collapsed on peptide position

Site based format "Site_peptidoform_centric_format_no_pApeptidoforms.py":

	takes "binomial_collapsed_FLR.csv"
	-filter for 0.05FLR threshold
	-remove no choice peptides
	-remove peptidoforms with pA
	-collapse on peptidoform modification site

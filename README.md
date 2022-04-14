# mzid_pipeline
 
Scripts for estimating false localisation rates (FLR) using two methods, model FLR and Decoy amino acid FLR, starting from mzid files (pep.xml pipeline [here](https://github.com/PGB-LIV/PhosphoFLR)). 

Also generates peptidoform site based format. 

# Usage

*TPP_reusable* folder contains supported module code for the analysis of Datasets using TPP PTMprophet mzid output files.

      $py TPP_comparison params_TPP	
    The params file (params_TPP) must be given and must contain the directories of each of the search file locations. 


****Site_peptidoform_centric_format_no_pApeptidoforms.py****

Used for creating site based peptidoform format.
 - Takes _binomial_collpased_FLR.csv_ as input

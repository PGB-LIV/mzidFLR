
#Working directory
p="D:/Dropbox/PTMExchange/Rice/New_build_ID/PXD000923"
#Working directory of each search to be compared, if more than one search
TPP_wd = p  #if decoy amino acid not specified in folder name (ie pL/pLeu), will automatically use pAla method
#location of search database
database="Z:/Rice/Osativa_super_annotation_union/Osativa_super_annotation_union_noIC4R_v2_cRAP_decoy.fasta"
#PXD ID
PXD = "PXD000923"
#FDR cutoff (0.0-1.0)
FDR_cutoff = 0.01
#List all working directories (if more than one, w=[TPP_wd,TPP_wd2...])
w = [TPP_wd]
#list all search names (if more than one, s=["TPP_tryptic","TPP_semitryptic"...])
s = ["TPP_tryptic"]


#Working directory
p="D:/Dropbox/PTMExchange/Rice/New_build_ID/PXD000923"
#p="D:/Dropbox/PTMExchange/Plasmodium/PXD026474/"
#Working directory of each search to be compared
TPP_wd = p  #if decoy amino acid not specified in folder name (ie pL/pLeu), will automatically use pAla method
#location of search database
#database = "Z:/Plasmodium/PlasmoDB-51_Pfalciparum3D7_AnnotatedProteins_cRAP_PA_THISP_Level1_2020-09-01_target_decoy.fasta"
database="Z:/Rice/Osativa_super_annotation_union/Osativa_super_annotation_union_noIC4R_v2_cRAP_decoy.fasta"
#PXD ID
PXD = "PXD002756"
#FDR cutoff (0.0-1.0)
FDR_cutoff = 0.01
#List all working directories
w = [TPP_wd]
#list all search names
s = ["TPP_tryptic"]
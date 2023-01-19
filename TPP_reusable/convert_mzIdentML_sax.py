import xml.sax
import csv
import pandas as pd
import os
import sys

csv.field_size_limit(100000000)

def convert(input):
    class TPPHandler(xml.sax.ContentHandler):
        def __init__(self):
            self.CurrentData = ""
            self.Seq = ""
            self.Counter=0

        def endElement(self, tag):
            if self.CurrentData == "PeptideSequence":
                output_peptidedict.write(self.Seq)
            self.CurrentData = ""
            self.Seq=""

        def characters(self, content):
            if self.CurrentData == "PeptideSequence":
                self.Seq += content

        # Call when an element starts
        def startElement(self, tag, attributes):
            #Add to peptide dict: pepid to peptide, mod and pos
            if tag == "Peptide":
                output_peptidedict.write("\n"+attributes['id']+",")
            if tag == "Modification":
                output_peptidedict.write(":"+attributes['location']+";")
                if 'residues' in attributes:
                    output_peptidedict.write(attributes['residues']+";")
                else:
                    output_peptidedict.write("N-term;")
                if 'monoisotopicMassDelta' in attributes:
                    output_peptidedict.write(attributes['monoisotopicMassDelta']+";")
                else:
                    output_peptidedict.write(";")
            if tag == "cvParam" and attributes['cvRef']=="UNIMOD":
                output_peptidedict.write(attributes['name']+";")

            #add to protein dict: pepid to dbseqref, position, +1/-1
            if tag == "PeptideEvidence":
                output_proteindict.write("\n"+attributes['dBSequence_ref']+",")
                output_proteindict.write(attributes['peptide_ref']+",")
                output_proteindict.write(attributes['start']+",")
                try:
                    output_proteindict.write(attributes['post']+",")
                except KeyError:
                    output_proteindict.write("N/A,")
                try:
                    output_proteindict.write(attributes['pre'])
                except KeyError:
                    output_proteindict.write("N/A,")

            #add to protein_id dict: dbseqref to protein
            if tag == "DBSequence":
                output_proteinIDdict.write("\n"+attributes['id']+",")
                output_proteinIDdict.write(attributes['accession'])


            self.CurrentData = tag
            if self.CurrentData == "SpectrumIdentificationResult":
                self.Counter = 0
                output.write("\n")
            if self.CurrentData == "SpectrumIdentificationItem":
                self.Counter += 1
            # add to temp csv
            if tag == "SpectrumIdentificationResult":
                output.write(attributes["name"] + ",")
            # Only take first identification
            if self.Counter == 1:
                if tag== "SpectrumIdentificationItem":
                    output.write(attributes["chargeState"] +",")
                    output.write(attributes["calculatedMassToCharge"] + ",")
                    output.write(attributes["experimentalMassToCharge"]+",")
                    output.write(attributes["peptide_ref"] + ",")
                if tag == "cvParam" and attributes["name"]=="Comet:expectation value":
                    output.write(attributes["value"] + ",")
                if tag == "cvParam" and attributes["name"]=="PSM-level probability":
                    output.write(attributes["value"]+",")
                if tag=="cvParam" and attributes["name"]=="PTMProphet probability":
                    output.write(attributes["value"] + ";")
            if tag=="cvParam" and attributes["name"]=="retention time":
                output.write(","+attributes["value"])


    pep_dict_name='IDML_peptidedict.csv'
    protein_dict_name='IDML_proteindict.csv'
    proteinID_dict_name='IDML_proteinIDdict.csv'
    temp_file_name="IDML2CSV_temp.csv"
    final_file_name=input+".csv"

    output = open(temp_file_name, "w")
    output.close()
    output_peptidedict = open(pep_dict_name,"w")
    output_peptidedict.close()
    output_proteindict = open(protein_dict_name,"w")
    output_proteindict.close()
    output_proteinIDdict = open(proteinID_dict_name,"w")
    output_proteinIDdict.close()
    output = open(temp_file_name, "a")
    output_peptidedict = open(pep_dict_name,"a")
    output_proteindict = open(protein_dict_name,"a")
    output_proteinIDdict = open(proteinID_dict_name,"a")
    parser = xml.sax.make_parser()
    Handler = TPPHandler()
    parser.setContentHandler(Handler)
    parser.parse(input)
    output.close()
    output_peptidedict.close()
    output_proteindict.close()
    output_proteinIDdict.close()

    with open(pep_dict_name, mode='r') as infile:
        next(infile)
        reader = csv.reader(infile)
        pep_dict = {rows[0]:rows[1] for rows in reader}

    with open(protein_dict_name, mode='r') as infile:
        next(infile)
        reader = csv.reader(infile)
        protein_dict = {rows[1]:str(rows[0]+";"+rows[2]+";"+rows[3]+";"+rows[4]) for rows in reader}

    with open(proteinID_dict_name, mode='r') as infile:
        next(infile)
        reader = csv.reader(infile)
        proteinID_dict = {rows[0]:rows[1] for rows in reader}

    spectrum=[]
    z=[]
    calc_mass=[]
    mass_exp=[]
    peptide=[]
    mods=[]
    pos=[]
    mass=[]
    res=[]
    protein=[]
    protein_pos=[]
    evalue=[]
    psm_prob=[]
    ptm_prob=[]
    rt=[]

    file = open(temp_file_name, 'r')
    for line in file:
        if "," in line:
            data_row = line.strip().split(',')
            spectrum.append(data_row[0])
            z.append(data_row[1])
            calc_mass.append(data_row[2])
            mass_exp.append(data_row[3])
            peptide_temp=data_row[4]
            peptide.append(pep_dict[peptide_temp].split(":")[0])
            if ";" in pep_dict[peptide_temp]:
                mod_list=pep_dict[peptide_temp].split(":")[1:]
            else:
                mod_list=[]
            mod_names=""
            mod_pos=""
            mod_mass=""
            mod_res=""
            for i in range(0,len(mod_list)):
                mod_pos+=mod_list[i].split(";")[0]+";"
                mod_res+=mod_list[i].split(";")[1]+";"
                #if mod_list[i].split(";")[2]=="":
                #    mod_mass+="0;"
                #else:
                mod_mass+=mod_list[i].split(";")[2]+";"
                if mod_list[i].split(";")[3]=="":
                    mod_names+="unknown_mod;"
                else:
                    mod_names+=mod_list[i].split(";")[3]+";"
            mods.append(mod_names)
            pos.append(mod_pos)
            mass.append(mod_mass)
            res.append(mod_res)
            protein_temp=protein_dict[peptide_temp].split(";")[0]
            protein_pos.append(protein_dict[peptide_temp].split(";")[1])
            #protein_pre=protein_dict[peptide_temp].split(";")[2]
            #protein_post=protein_dict[peptide_temp].split(";")[3]
            protein.append(proteinID_dict[protein_temp])
            evalue.append(data_row[5])
            psm_prob.append(data_row[6])
            #start at 7, increase check for blank or ":" (ptm prophet), if blank next cell is RT, if ptm prophet result, next cell is RT
            #only taking first peptide for each spectrum
            for i in range(7,(len(data_row))):
                if ":" in data_row[i]:
                    PTM_info_temp=[]
                    PTM_info_list=""
                    for p in data_row[i][:-1].split(";"):
                        if p.split(":")[2] not in PTM_info_temp:
                            PTM_info_temp.append(p.split(":")[2])
                            PTM_info_list+=p+";"
                    ptm_prob.append(PTM_info_list[:-1])
                    rt.append(data_row[i+1])
                    break
                elif data_row[i]=="":
                    rt.append(data_row[i+1])
                    ptm_prob.append("")
                    break

    try:
        df = pd.DataFrame({"Spectrum":spectrum,"Charge":z,"Calculated mass":calc_mass,"Experimental mass":mass_exp,"Peptide":peptide,"Modifications":mods,"Positions":pos,"Modification mass":mass,"Modification residue":res,
                       "Protein":protein,"Protein position":protein_pos,"e-value":evalue,"PSM probability":psm_prob,"PTM info":ptm_prob,"Retention time":rt})

        df.to_csv(final_file_name, index=False)
        file.close()
        os.remove("IDML_peptidedict.csv")
        os.remove("IDML_proteindict.csv")
        os.remove("IDML_proteinIDdict.csv")
        os.remove("IDML2CSV_temp.csv")
    except:
        file.close()
        os.remove("IDML_peptidedict.csv")
        os.remove("IDML_proteindict.csv")
        os.remove("IDML_proteinIDdict.csv")
        os.remove("IDML2CSV_temp.csv")
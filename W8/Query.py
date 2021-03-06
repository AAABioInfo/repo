#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 23:23:27 2020

@author: Aldo Amaya
"""


#Import all of these

import argparse
import mygene
import ensembl_rest
import re
import Bio
import requests


#Import your variables 

if __name__ == "__main__":
    import_args = argparse.ArgumentParser(description='Which gene are you interested in?')
    import_args.add_argument("-g",  "--Gene", required=True, help="Place gene name here")
    import_args.add_argument("-f",  "--Fields", required=True, help="Which field")
    import_args.add_argument("-s",  "--Species", required=True, help="Place species query")
    import_args.add_argument("-e",  "--email", required=True, help="e-mail for entrez search")

    i_args = import_args.parse_args()

ref_gene= i_args.Gene
ref_field= i_args.Fields
ref_SPC= i_args.Species
email= i_args.email



mg = mygene.MyGeneInfo()

#use this query to find all matches for gene of interest
#for testing
ref_gene='MC1R'
ref_SPC='human'
ref_field='ensembl.gene'
email='aamaya@brandeis.edu'

######## This is your ensembl I
out = mg.querymany(ref_gene, scopes='symbol', fields=ref_field, species=ref_SPC)
ENSB_ID=out[0]['ensembl']['gene']
refine=mg.query(ENSB_ID)
FullName=refine['hits'][0]['name']
Symbol=FullName=refine['hits'][0]['symbol']


######## Grab more data from ensembl 

ESBL_data=ensembl_rest.sequence_id(ENSB_ID) #dictionary
ESBL_data.keys()

Molecule=ESBL_data['molecule']
if Molecule == 'dna':
    print('This is a dna sequence')
else:
    pass

FASTA=ensembl_rest.sequence_id(ENSB_ID,headers={'content-type': 'text/x-fasta'}) #dictionary

ref_gene='MC1R'
ref_SPC='human'
ref_field='ensembl.gene'


###### Lets make a fasta file
Filename=ref_gene+"_"+ref_SPC+".fasta"
f=open(Filename,'w') 
for item in FASTA:
    f.write(item)
f.close()

#####find the longest ORG

LongestORF=max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',FASTA), key = len)

#####convert from DNA to AA

from Bio.Seq import Seq
my_seq = Seq(LongestORF)
my_aa=my_seq.translate()
extractedAA = ""
for i in my_aa:
    extractedAA+=i

#### Write to Fasta File

with open(Filename, "a+") as f:
    f.write(">Amino Acid conversion from Longest ORF of %s| \n" % (ref_gene))
    f.write(extractedAA)
    f.write('"\n')
f.close()


######find homology list


ensembl_server = 'http://rest.ensembl.org'

###Support script
def do_request(server, service, *args, **kwargs):
    url_params = ''
    for a in args:
        if a is not None:
            url_params += '/' + a
    req = requests.get('%s/%s%s' % (server, service, url_params), params=kwargs, headers={'Content-Type': 'application/json'})
    if not req.ok:
        req.raise_for_status()
    return req.json()

hom_response = do_request(ensembl_server, 'homology/id', ENSB_ID, type='orthologues', sequence='none')
homologies = hom_response['data'][0]['homologies']
Homology_Ls=[]
for homology in homologies:
    temp_Orthologe=homology['target']['species']
    Homology_Ls.append(temp_Orthologe)

##### Lets print to new document

Filename2=ref_gene+'_homology_list.txt'
f=open(Filename2,'w') 
for item in Homology_Ls:
    f.write(item+"\n")
f.close()


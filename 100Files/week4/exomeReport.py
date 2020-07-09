#!/usr/bin/python

#This script will create a report for Dr. Banner
#inputs: Sequenced_20r30mm.txt '/Users/brownbear/sbx_bash/W4_AAA/exomes' (need absolute dir of where python script is present)


#Output: Will display for each speices the following: 
#Sequenced_20r30mm.txt - this will be used to extract the species from All_exomes



# Comment out for testing
#Largs= ['test.py', 'Sequenced_20r30mm.txt','/Users/brownbear/sbx_bash/W4_AAA/exomes']
# test
import sys
Largs=[]
for largs in sys.argv:
    Largs.append(largs)


import os


if os.getcwd()==Largs[-1]:
    pass
else:
    os.chdir(Largs[-1])
#used for testing directory
#os.chdir("..")
#os.getcwd()
#os.listdir()
#os.chdir(.\)


# this will make a list and dictionary from the clinical data
#there will be two outputs 

#Species=[] This will create a sorted list of species
#Clin_list2={} This will create a dictionary, will append uniq list and multiple pairs 

import csv
Species=[]
Clin_list2={}


with open(Largs[1], newline = '') as SelectedSub:                                                                                          
    	SS_reader = csv.reader(SelectedSub, delimiter='\t')
    	for Subjects in SS_reader:         
            Clin_list2[Subjects[-1]]= Subjects[0:-1]
            Species.append(Subjects[-1])

del Species[0]
Species.sort()
tagSpecies=list(range(0, len(Species)))
uniqlists = [ [] for _ in range(len(tagSpecies)) ] # this will be used for 
ubilists = [ [] for _ in range(len(tagSpecies)) ] # this will be used for 
#attach empty list to attach unique genes 
Tagged_SP=[Species,tagSpecies,uniqlists,ubilists] 

# this function uses os imports

def ProcessDir(ind_Species, group):
    """ This function will create list [containing gene, seq, tag] """  
    os.chdir(ind_Species)
    os.getcwd()
    SearchL = os.listdir() # This will create a list of whats available
    import fnmatch
    pattern='*precrispr.fasta'
    FT_matching = fnmatch.filter(SearchL, pattern)
    FT_matching  # this is the list of variables that need to be combined.
    
    
    #append data sort and determien duplicate
    combined_SFasta=[]
    gene_num_app=[]
    gene_seq_app=[]  
  
    for app in FT_matching: # this will append all precris genes for species
        with open(app, newline = '') as SelectedSubF:                                                                                          
           FAreader = csv.reader(SelectedSubF, delimiter=' ')
           for FastaR in FAreader: 
               if  ">" not in FastaR[0]:
                   #print ('sequence')
                   if FastaR[0] not in gene_seq_app:
                         gene_seq_app.append(FastaR[0])
                   else: 
                         print ('duplicate seq')
                        
               else:
                  #print ('gene')
                  if FastaR[0] not in gene_num_app:
                          gene_num_app.append(FastaR[0])
                  else: 
                         print ('duplicate gene')
                 
    #tag data               
    groupT=[group] * len(gene_num_app)
    combined_SFasta=[gene_num_app, gene_seq_app, groupT]   
    os.chdir("..")
    os.getcwd()
    # and then send as output 
    return combined_SFasta         
               
#make a container to house all sequences in order of tagged species
Container=[]

for each_Spe in Tagged_SP[0]:    
    group=Tagged_SP[1][Tagged_SP[0].index(each_Spe)]
    ind_Species=each_Spe
    
    TempVar=ProcessDir(ind_Species,group)
    Container.append(TempVar)

#combine container      
C_Container=[] # create 3 arrays that are indexed. 

All_C_genes=[]
All_C_seq=[]
All_C_marker=[]

for each_SpeCont in Container:
    counter=0
    while counter < len(each_SpeCont):
        if counter==1:
            All_C_genes.extend(each_SpeCont[0])
        elif counter==2:
            All_C_seq.extend(each_SpeCont[1])
        else:
            All_C_marker.extend(each_SpeCont[2])           
        counter=counter+1
             
C_Container=[All_C_genes, All_C_seq, All_C_marker]      
   

#determine unique genes     
DupGenes_perS=[] # this means when genes were imported from diff motifs, there was a dup

def getIndexPositions(listOfElements, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    indexPosList = []
    indexPos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            indexPos = listOfElements.index(element, indexPos)
            # Add the index position in list
            indexPosList.append(indexPos)
            indexPos += 1
        except ValueError as e:
            break
 
    return indexPosList


# this will parse gene name list
for compSeq in All_C_seq:
    rep_tempIndex=getIndexPositions(All_C_seq,compSeq) #how many times is seq present
    if len(rep_tempIndex)==1:
        tempI=rep_tempIndex[0]
        ind_tag=Tagged_SP[1].index(C_Container[2][tempI])
        Tagged_SP[2][ind_tag].append(C_Container[0][tempI])
        
    elif len(rep_tempIndex)==7: 
        print('this is shared by all species')
        tempI=rep_tempIndex[0]
        ind_tag=Tagged_SP[1].index(C_Container[2][tempI])
        Tagged_SP[-1][ind_tag].append(C_Container[0][tempI])
        
    elif 1<len(rep_tempIndex)<7:
        pass
    else:
        print('You have a duplicate gene')
        tempI=rep_tempIndex[0]
        DupGenes_perS.append(C_Container[0][tempI])

#DupGenes_perS   
#lets make a list of all unique genes with tags before we assign to dic
        
# print final report using numpy

#os.chdir("..")
#os.getcwd()
#os.listdir()

import numpy as np

f=open('Final_exomeReport.txt','ab')

with open("Final_exomeReport.txt", "ab") as f:

    for prt_spe in Tagged_SP[0]:
        f.write(b"\n")
        tempind=Tagged_SP[0].index(prt_spe)
        DISCOVERER=Clin_list2[prt_spe][0]
        DIAMETER=Clin_list2[prt_spe][2]
        environment=Clin_list2[prt_spe][3]
        location=Clin_list2[prt_spe][1]
        if len(ubilists[tempind])==0:
           Tubilists="0"
        else:
            Tubilists=ubilists[tempind]
        # Creating the type of a structure
        structuredArr=np.array([("Organism:",prt_spe),("Discovered by:",DISCOVERER),("Has a diameter(mm) of:",DIAMETER),("Found at:",environment), ("In this Location:", location), ("The following genes unique to only itself:", uniqlists[tempind]), ("The following genes ubiquitous with all species in this list:",Tubilists)])
        np.savetxt(f, structuredArr, delimiter=' ', fmt=['%s' , '%s'])   
     
f.close()
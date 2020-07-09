#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 23:23:27 2020

@author: Aldo Amaya
"""
# mtools sort -m 100M -o {name}.sorted.bam {name}.bam

#### This set of the script will take in all your arg
#### + verify pwd to correct folder

#####################used for testing
# Largs=['pipeline.py', 'harrington_clinical_data.txt','dgorgon_reference.fa', 'hawkins_pooled_sequences.fastq']
####################

#Import all of these
import os
import sys
import pandas as pd
import numpy as np
import re
import csv
import argparse
import gzip
import subprocess
import pysam
import glob



#Import your variables 

if __name__ == "__main__":
    import_args = argparse.ArgumentParser(description='import all arguments for this script')
    import_args.add_argument("-c",  "--ClinData", required=True, help="Place clinical data")
    import_args.add_argument("-rs",  "--ref_Seq", required=True, help="Place reference sequence")
    import_args.add_argument("-fq",  "--fastq", required=True, help="Place your fasq file_pooled")
    i_args = import_args.parse_args()

clinical_data= i_args.ClinData
ref_Seq= i_args.ref_Seq
Pooled_fastq= i_args.fastq

#For testing
# clinical_data=Largs[1]
# ref_Seq=Largs[2]
# Pooled_fastq=Largs[-1]


####################################Demultiplex####################################
###################################################################################
###################################################################################
#########lets crete a library with the clinical data############## 

ClinL = pd.read_csv(clinical_data,sep='\t') 
ClinLD=ClinL.set_index('Name').T.to_dict('list')
ClinLR=ClinL.set_index('Name').T.to_dict('list')

### To test whats there : 
# In:ClinL.columns
# Out:: Index(['Name', 'Color', 'Barcode'], dtype='object')

#print(ClinLD)
####################### Let's import reference ##############
GeneR=[]
SeqR=[]

with open(ref_Seq, newline = '') as SelectedSub:                                                                                          
    	SS_reader = csv.reader(SelectedSub, delimiter='\t')
    	for ref in SS_reader: 
            if re.search('>', ref[0]):
                #print("this was a Ref_gene")
                GeneR.append(ref[0])
            else:
                #print("this was a Ref_sequence")
                SeqR.append(ref[0])
                
ReferenceL=[GeneR,SeqR]    
#print(ReferenceL)
################################################
# This class was provided by Tim Song to Parse pooled fastq files
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.__next__()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        return self
     
    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)
##########################################################################

Parsed_fromPool=ParseFastQ(Pooled_fastq)

seqHeader=[]
seqStr=[]
qualHeader=[]
qualStr=[]

for fastq_obj in Parsed_fromPool:
    seqHeader.append(fastq_obj[0])
    seqStr.append(fastq_obj[1])
    qualHeader.append(fastq_obj[2])
    qualStr.append(fastq_obj[3]) 

Prs_PoolL=[seqHeader,seqStr,qualHeader,qualStr]
 
# ###############This is filtered by barcodes##############
Pooled_ID=Prs_PoolL[0]
Pooled_Seq=Prs_PoolL[1]
Pooled_QC=Prs_PoolL[3]

PbyBC_data=[]                  
for key in ClinLD:
    temp_BarCode=ClinLD[key][1]
    Temp_IDMatch=[]
    temp_Matches=[]
    Temp_QCMatch=[]
    TempContainer=[]
    for pMatch in Pooled_Seq:
        if pMatch[0:5] == temp_BarCode:
            temp_Matches.append(pMatch)
            indexL=Pooled_Seq.index(pMatch)
            Temp_QCMatch.append(Pooled_QC[indexL])
            Temp_IDMatch.append(Pooled_ID[indexL])  

        else:
            pass
    
    TempContainer=[temp_BarCode,temp_Matches, Temp_QCMatch,Temp_IDMatch]
    PbyBC_data.append(TempContainer)



# ###############clean barcodes and update qc##############

rBC_data=[]        
for setdata in PbyBC_data:
    # PbyBC_data[]
    temp_BarCode=setdata[0]
    tclearedBC=[]
    TempContainer=[]
    updated_tqcdata=[]
    delmark=len(temp_BarCode)
    InsidePdata=setdata[1]
    tqcdata=setdata[2]
    tIDs=setdata[3]
    for clearb in InsidePdata:
        TempClear=clearb.replace(temp_BarCode,"", 1)
        tclearedBC.append(TempClear)
    updated_tqcdata = [e[delmark:] for e in tqcdata] #this will erase barcode size from qc   
    TempContainer=[temp_BarCode,tclearedBC,updated_tqcdata,tIDs]
    rBC_data.append(TempContainer)


# # print(len(rBC_data[38][2][0]))
# # print(len(PbyBC_data[38][2][0]))

# ###############QC reads, after degradation##############


TQC_data=[]
TempContainer=[] 
for qcread in rBC_data:
    #set of alias qcread=rBC_data[0]
    barcode=qcread[0]
    rBCseq=qcread[1]
    qcdata=qcread[2]
    qcIDs=qcread[3]
    pMatch=['DD',"FD","FF","DF"]
    # list of new variables
    seqqc=[]
    qc_updated=[]
    index_seqdel=[]    
    TempContainer=[]
    #Create an index on qhere the quality drops per sequence
    for qc_pf in qcdata:
        # qc_pf=qcdata[0]
        temp_ecmdL=[]
        for p in pMatch:
            tempMatch=re.search(p, qc_pf)
            if tempMatch==None:
                tempMatch=[]
            else:
                 temp_ecmdL.append(tempMatch.span()[0])    
        vardel=min(temp_ecmdL)
        tqc_pf=qc_pf[:vardel]
        index_seqdel.append(vardel)
        qc_updated.append(tqc_pf)
    TempContainer=[barcode, rBCseq, index_seqdel, qc_updated, qcIDs]
    TQC_data.append(TempContainer)
 
######update the sequences for degradation
QC_data=[]
TempContainer=[]
for d in TQC_data:
    barcode=d[0]
    rBCseq=d[1]
    index_seqdel=d[2]
    qc_updated=d[3]
    qcIDs=d[4]

    # print(len(index_seqdel[2]))

    idx_c=0
    while idx_c <len(index_seqdel):
            vardel=index_seqdel[idx_c]
            tempUpdate=rBCseq[idx_c][:vardel]
            seqqc.append(tempUpdate)
            idx_c=idx_c+1    
    TempContainer=[barcode,seqqc,qc_updated,qcIDs]
    QC_data.append(TempContainer)

############## Tag new fastq to dir before export##############
###for testing
#ClinLD=ClinL.set_index('Name').T.to_dict('list')

for qcL in QC_data:
    tempBC=qcL[0]
    for key, value in ClinLD.items():
        if tempBC in value:
            ClinLD[key].append(qcL[1:])


############## iterate through clin data and make new dir+fast ##############
os.mkdir('fastqs')

for key in ClinLD:
    tempdic=ClinLD[key]
    keyname=key
    tempseq=tempdic[2][0]
    tempqc=tempdic[2][1]
    tempid=tempdic[2][2]
    tempfilename=keyname+'_trimmed.fastq'
    TempSavelocation="./fastqs/"+tempfilename
    f=open(TempSavelocation,'ab')
    icounter=0
    while icounter < len(tempid): 
        with open(TempSavelocation,'ab') as f:
            # Creating the type of a structure
            structuredArr=np.array([tempid[icounter],tempseq[icounter],"+",tempqc[icounter]])
            np.savetxt(f, structuredArr, fmt=['%s'])  
        icounter=icounter+1 
    f.close()

    
###################################################################################
###################################################################################
###################################################################################

############Perform alignment on each FASTQ to the reference sequence##############
#samtools 1.10
#Using htslib 1.10.2
#ref_Seq, Out: 'dgorgon_reference.fa'



#index your reference
RefPath='bwa index '+ref_Seq
os.system(RefPath)

#find all trimmed_fast from Fastq file and create a list
FastTrimL=os.listdir('fastqs')

# make a definition that process each file into the fastqs folder

def sam_2_sortbam (trimmed_fasq):
    Trimfile=trimmed_fasq
    name=Trimfile.split('_')[0]
    dirRefSeq='../'+ref_Seq
    # create a sam file
    pathSam='bwa mem '+dirRefSeq+ " " +Trimfile+' >'+name+'.sam'
    os.system(pathSam)  
    # create a bam file from a sam file
    pathbam='samtools view -bS '+name+'.sam > '+name+'.bam'
    os.system(pathbam)  
    # sort your bam files from a bam files
    path_sort_bam='samtools sort -m 100M -o ' +name+'.sorted.bam '+name+'.bam'
    os.system(path_sort_bam)
    # index your bam files
    path_idxsort_bam='samtools index '+name+'.sorted.bam'
    os.system(path_idxsort_bam)

# go into dir 
os.chdir('./fastqs')    
  
for FTL in FastTrimL:
    sam_2_sortbam(FTL)


#make all subforders 
path='bams'
os.makedirs(path)
####clean up
os.system('mv *sorted* bams/') 
os.system('cp -r bams ..') 
os.system('rm -r bams') 
os.system('rm *sam *bam') 
os.chdir("..")
os.makedirs('bwaindex')
os.system('mv *.fa.* bwaindex/') 
os.system('clear') 


############Variant Discovery ##############
IdxBamL = glob.glob('bams/*.bam')

#### we will make a reporting dict, pair clinical ID to the SNP, this was made at start 
###for test
#ClinLR=ClinL.set_index('Name').T.to_dict('list')
ClinLR

####define pileup(), will return a dic  FPile={int-key-read: list of nuclotides}
def pileup(idxfile):
    
    
    #test file, replaced with the sorted.bam you are using. Make sure it is indexed! (Use samtools index yourbam.sorted.bam)
    dirFile=idxfile
    samfile = pysam.AlignmentFile(dirFile, "rb")
    dictperfasta={}
    #Since our reference only has a single sequence, we're going to pile up ALL of the reads. Usually you would do it in a specific region (such as chromosome 1, position 1023 to 1050 for example)
    for pileupcolumn in samfile.pileup():
        #print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        #use a dictionary to count up the bases at each position
        # this says okay make a dictionary that temporary houses each of your reads 
        ntdict = {}
        alist=[]
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # You can uncomment the below line to see what is happening in the pileup. 
                #print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                ########## ADD ADDITIONAL CODE HERE ############# 
                alist.append(base)
        ntdict = {pileupcolumn.pos:alist} 
        dictperfasta.update(ntdict)
    samfile.close()
    
    return dictperfasta


for idxfile in IdxBamL:
    #let determine name in order to pair to clin data
    removedirID=idxfile.split("/")[1]
    name=removedirID.split(".")[0]
    #this will use pileup() to read indexed Bam files agains reference and find reads by position
    FPile=pileup(idxfile)
    # This will sort through all reads and for each will present the following:
    # 'read#_Totalreads_BaseID_Percentage.present'  
    
    BreakdownL=[]
    for read in FPile:
        WholeL=FPile[read]
        nreads=len(WholeL)
        uniqbases=sorted(set(WholeL))
        idxbases=[]
        for uniqbs in uniqbases:
            occurances=WholeL.count(uniqbs)
            percentOcr=int((occurances/len(WholeL))*100)
            tempuniq=str(read)+"_"+str(occurances)+"_"+uniqbs+"_"+str(percentOcr)
            idxbases.append(tempuniq)
        BreakdownL.append(idxbases)
    
    #This will filter the reads that do not have SNP
    AppL=[]
    for l in BreakdownL:
        if len(l) > 1:
            AppL.append(l)
        else:
            pass
        
    # this will append to ClinLR 
    ClinLR[name].append(AppL)
     

############Lets determine pool data for color ##############

##determine unique colors
All_colors=[];
for i in ClinLR:
    colorfind=ClinLR[i][0]
    All_colors.append(colorfind)
    
unqColors=sorted(set(All_colors))

#This will store the data for each color
#list Id ['Black', 'Green', 'Orange', 'Yellow']
Colordata=[ [] for _ in range(len(unqColors)) ]

    
#pool SNP reads
for ind in ClinLR:
    colorfind=ClinLR[ind][0]
    colorstats=ClinLR[ind][2][0]
    # Where to store
    savehere=unqColors.index(colorfind)
    dataparse=[] 
    temppos=[]
    tempred=[]
    tempbase=[]
    temppercent=[]
    for stats in colorstats:
        pdata=stats.split("_")
        temppos.append(int(pdata[0]))
        tempred.append(int(pdata[1]))
        tempbase.append(pdata[2])
        temppercent.append(pdata[3])

    #set the index for WT based on Reference sequence
    SNP_postion=sorted(set(temppos))[0]
    WT_from_Ref=SeqR[0][SNP_postion]
    
    if tempbase.index(WT_from_Ref)==0:
        Mutindex=1
    else:
        Mutindex=0
    #set the bases/ reads/ total/ SNP percent 
    Wild_Type_base=tempbase[tempbase.index(WT_from_Ref)]
    Mutation_base=tempbase[Mutindex]
    Wild_reads=tempred[tempbase.index(WT_from_Ref)]
    Mutation_reads=tempred[Mutindex]
    Total_Reads=Wild_reads+Mutation_reads
    SNP_percent=temppercent[Mutindex]
    #append to dictionary

    dataparse=[Total_Reads, SNP_percent, SNP_postion, Mutation_base]
    ClinLR[ind].append(dataparse)
    
    #append to color only 
    colorappend=[SNP_postion, Wild_Type_base,Wild_reads, Mutation_base, Mutation_reads]
    Colordata[savehere].append(colorappend)
       
        

############Lets Print based on sample first ##############
f=open('report.txt','ab') 
for key in ClinLR:
    
    tempdic=ClinLR[key]
    color="Color: " +ClinLR[key][0]
    keyname="Sample: "+key
    tempReportD=tempdic[3]
    Reads="Total Reads: "+ str(tempReportD[0])
    Percent="Percentage different from Wildtype: "+str(tempReportD[1]+"%")
    Location="Mutant location in respect to referance: "+str(tempReportD[2])
    MutantID="Mutation: "+tempReportD[3]
    
    
    with open('report.txt','ab') as f:
            # Creating the type of a structure
            Report=np.array([keyname,color,Reads,Location, MutantID, Percent])
            np.savetxt(f, Report, fmt=['%s'])  
            f.write(b"\n")
f.close()


############Now lets print the colors ##############

# Colordata
# unqColors

f=open('report.txt','ab') 
colori=0
while colori <len(unqColors) :
    tempColor=unqColors[colori]
    tempLocation=Colordata[colori][0][0]
    tempWT=Colordata[colori][0][1]
    tempMUT=Colordata[colori][0][3]

    color="Mold Color: "+tempColor
    Location="Caused by a mutation in position: "+str(tempLocation)
    WT="The wildtype base was: "+tempWT
    MUT="The mutation base was: "+tempMUT
    with open('report.txt','ab') as f:
            # Creating the type of a structure
            Report=np.array([color,Location,WT,MUT])
            np.savetxt(f, Report, fmt=['%s'])  
            f.write(b"\n")
    colori=colori+1
f.close()

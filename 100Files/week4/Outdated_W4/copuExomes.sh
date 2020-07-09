#!/bin/bash
# input 1: $1= clinical_data.txt
# Output: Sequenced_20r30mm.txt
#Creater Aldo Amaya
#Purpose: Take clinical data and identiy samples that are betweent 20 and 30 mm long (inclusive) and have been sequenced,


cat $1 >sourceL   
mkdir exomes
#cat clinical_data.txt >sourceL   
#Purpose 
  #Read in the clinical data file and remove all spaces from entries with spaces and save them as new file Formatted_clinical_data.txt
sed 's/ //g' sourceL> Formatted_clinical_data.txt
  #identify the samples that have a diameter between 20 and 30 mm long (inclusive) 
head -n1 Formatted_clinical_data.txt>Sequenced_20r30mm.txt 
egrep "[3][0]" Formatted_clinical_data.txt>Samples_20r30mm.txt
egrep "[2][0-9]" Formatted_clinical_data.txt>>Samples_20r30mm.txt

  #and have had their genomes sequenced. 
  
egrep "Sequenced" Samples_20r30mm.txt>>Sequenced_20r30mm.txt
  #Copy the identified exomes using the sample code names to a new directory called exomes.
  
cp Sequenced_20r30mm.txt cd exomes/   
rm Sequenced_20r30mm.txt Samples_20r30mm.txt  Formatted_clinical_data.txt  sourceL

echo "All clear, have a nice day"


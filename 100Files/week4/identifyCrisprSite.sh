#!/bin/bash
# input $1=All_exomes (target folder will need to have "All_exomes_ID")
# Output: exomename_topmotifs.fasta- this iwll do the top 3 , will also create a summary for the species species_precrispr_summary
#Creater Aldo Amaya
#Purpose: For each gene inside the exomename_topmotifs.fasta files, identify a suitable CRISPR site.


TargetSeq=$(echo $1)
#TargetSeq=$(echo All_exomes/)
cd $TargetSeq

# This is the list of names of species to examine
cat All_exomes_ID>sourceL
lenSp=$(sort sourceL| uniq|wc -l) # echo $lenSp
i=1

# go into each folder and look through the top 3 motifs to look for Crisper sites 21-bp+GG
while [ $i -le $lenSp ]
do
  tnExoname=$(head -n $i sourceL| tail -n 1)
  cd $tnExoname
  ls *fasta>tempCL
  lentempCL=$(sort tempCL| uniq|wc -l) #echo $lentempCL
  ii=1

  while [ $ii -le $lentempCL ]
  do
    tmpTop=$(head -n $ii tempCL| tail -n 1)
    RtmpTop=$(echo $tmpTop| sed 's/.fasta//g')
    awk --re-interval '/[ACTG]{21,}GG/{print a; print}{a=$0}' $tmpTop > "${RtmpTop}_precrispr".fasta
    let ii=$ii+1
  done 
  
  #clean up
  rm tempCL 
  
  cd ..
    
let i=i+1

done

rm sourceL

clear 

echo "All clear, have a nice day"

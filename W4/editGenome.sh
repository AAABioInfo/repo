#!/bin/bash
## input $1=All_exomes (target folder will need to have "Top 3 precrisper files), $2= exomes
# Output: this will create the top3 motif postcrispr with an inserted A after NGG, will also create a summary for the species
#Creater Aldo Amaya
#Purpose: From top 3 motifs each possible crispr site will insert and A before NGG site 

TargetSeq=$(echo $1)
cd $TargetSeq

TargetMov=$(echo $2)

#TargetSeq=$(echo All_exomes/)
#TargetMov=$(echo exomes/) #echo $TargetMov
# This is the list of names of species to examine
cat All_exomes_ID>sourceL
lenSp=$(sort sourceL| uniq|wc -l)
i=1

# go into each folder and look through the top 3 motifs to look for Crisper sites 21-bp+GG
while [ $i -le $lenSp ]
do
  tnExoname=$(head -n $i sourceL| tail -n 1) #echo $tnExoname
  cd $tnExoname
  ls *precrispr.fasta>tempCL #cat tempCL
  lentempCL=$(sort tempCL| uniq|wc -l)
  ii=1
  
  while [ $ii -le $lentempCL ]
  do
    tmpTop=$(head -n $ii tempCL| tail -n 1)
    RtmpTop=$(echo $tmpTop| sed 's/precrispr.fasta//g')
    sed 's/.GG/A&/' $tmpTop> "${RtmpTop}postcrispr.fasta"
    let ii=$ii+1
  done 
  
  #clean up
  rm tempCL 

  cd ..
  
  
  cp -r $tnExoname ../$TargetMov 
    
let i=i+1

done

rm sourceL All_exomes_ID
clear

echo "All clear, have a nice day"

#ls -lt|head 


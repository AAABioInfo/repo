#!/bin/bash
# input no input required; however, target folder will need to have "Top 3 precrisper files"
# Output: exomename_topmotifs.fasta- this iwll do the top 3
#Creater Aldo Amaya
#Purpose: From top 3 motifs each possible crispr site will insert and A before NGG site 


TargetSeq=$(echo All_exomes/)
cd $TargetSeq

# This is the list of names of species to examine
cat All_exomes_ID>sourceL
lenSp=$(sort sourceL| uniq|wc -l)
i=1

# go into each folder and look through the top 3 motifs to look for Crisper sites 21-bp+GG
while [ $i -le $lenSp ]
do
  tnExoname=$(head -n $i sourceL| tail -n 1)
  cd $tnExoname
  ls *.fasta>tempCL
  lentempCL=$(sort tempCL| uniq|wc -l)
  ii=1
  
  while [ $ii -le $lentempCL ]
  do
    tmpTop=$(head -n $ii tempCL| tail -n 1)
    RtmpTop=$(echo $tmpTop| sed 's/precrispr.fasta//g')
    sed 's/.GG/A&/' $tmpTop> "$RtmpTop postcrispr".fasta
    let ii=$ii+1
  done 

  #remove all the spaces in all the file names
  for f in *\ *; do mv "$f" "${f// /}"; done
  
  #clean up
  rm tempCL 

  cd ..
    
let i=i+1

done

echo "All clear, have a nice day"


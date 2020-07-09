#!/bin/bash
# input no input required; however, target folder will need to have "All_exomes_ID"
# Output: exomename_topmotifs.fasta- this iwll do the top 3
#Creater Aldo Amaya
#Purpose: For each gene inside the exomename_topmotifs.fasta files, identify a suitable CRISPR site.
          
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
  ls *fa>tempCL
  lentempCL=$(sort tempCL| uniq|wc -l)
  ii=1
  touch "$tnExoname _precrispr_summary"
  while [ $ii -le $lentempCL ]
  do
    tmpTop=$(head -n $ii tempCL| tail -n 1)
    RtmpTop=$(echo $tmpTop| sed 's/.fa//g')
    echo $RtmpTop>>"$tnExoname _precrispr_summary"          # this will be the header to determine precrisper genes
    gawk '/[ACTG]{21,}GG/{print a; print}{a=$0}' $tmpTop > "$RtmpTop _precrispr".fasta
    gawk '/[ACTG]{21,}GG/{print a}{a=$0}' $tmpTop >>"$tnExoname _precrispr_summary"
    let ii=$ii+1
  done 
  
  #remove all the spaces in all the file names
  for f in *\ *; do mv "$f" "${f// /}"; done
  
  #clean up
  rm tempCL 

  cd ..
    
let i=i+1

done

rm sourceL

echo "All clear, have a nice day"

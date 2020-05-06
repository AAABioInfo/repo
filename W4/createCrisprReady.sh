#!/bin/bash
# input 1: $1= motif_list.txt $2=All_exomes
# Output: exomename_topmotifs.fasta- this iwll do the top 3
#Creater Aldo Amaya
#Purpose: Identify the 3 highest occurring motifs in each exome inside the exomes folder. 

            # this will link up any desired list 
#TargerSeq=$(echo $2)
cat $1 >sourceL
#used to test script
#cat motif_list.txt>sourceL                # this will link up any desired list 
TargetSeq=$(echo All_exomes/)


#create a list of all exomes to test and their names to add their extensions, without calling fasta file 
cd $TargetSeq
ls *fasta > All_exomes_Test 
sed 's/.fasta//g' All_exomes_Test>All_exomes_ID
cd ..
cat motif_list.txt>sourceL     
mv sourceL cd $TargetSeq 
#set parameters for loops 

cd $TargetSeq

lenM=$(sort sourceL| uniq|wc -l)
let i=1 # when using ls it counts All_exomes_Test 
lenExo=$(sort All_exomes_Test| uniq|wc -l)

# Go through each exoname and create a directory and store high motifs 
while [ $i -le $lenExo ] # go through list of all samples
do 

# This part will grab the first exoname and create a folder. Here you will store all of the top exo Motifs_count,i.e. fox_AAAA_103.fafsa,fox_AAAA_102.fafsa, fox_AAAA_101.fafsa  
  tempF=$(head -n $i All_exomes_Test| tail -n 1)  
  tnExoname=$(head -n $i All_exomes_ID| tail -n 1)
  mkdir $tnExoname
  cat $tempF>Lfasta.fa
  
  #set params for new loop to grep each motif to fasta file
  ii=1 
  #touch "$tnExoname _motif_count_All.txt"

  
  while [ $ii -le $lenM ]
  do
  tempM=$(head -n $ii sourceL| tail -n 1)
  countN=$(grep $(echo $tempM) Lfasta.fa|wc -l) 
  echo $tnExoname $tempM $countN >>"$tnExoname _motif_count_All.txt" 
  let ii=ii+1  
  done 

#copy each of these list into exoname folder
  cat "$tnExoname _motif_count_All.txt" | sort -k3nr| head -n 3 > "$tnExoname _motif_count_TOP3.txt"
  cp "$tnExoname _motif_count_TOP3.txt" cd $tnExoname
  cp "$tnExoname _motif_count_All.txt" cd $tnExoname

  #cd ..  ls -lt |head -n 5   ls -lt cd addax
  
  # sets a temporary list  to reprocess 
  awk '{print $2}' "$tnExoname _motif_count_TOP3.txt" > TempT3
  iii=1 # this will be referencing top list 
  let lenTT3=3   #lenTT3=3 
    while [ $iii -le $lenTT3 ]
    do
      tempM3=$(head -n $iii TempT3| tail -n 1)
      countN3=$(grep $(echo $tempM3) Lfasta.fa|wc -l)
      grep -B1 $(echo $tempM3) Lfasta.fa> "$tnExoname _$tempM3 _$countN3".fa
      cp "$tnExoname _$tempM3 _$countN3".fa $tnExoname
      rm "$tnExoname _$tempM3 _$countN3".fa
      let iii=iii+1  
    done 

#remove all the spaces in all the file names
  cd $tnExoname
  for f in *\ *; do mv "$f" "${f// /}"; done
  cd ..
  
  #Clean up after each motif 
  rm "$tnExoname _motif_count_TOP3.txt" "$tnExoname _motif_count_All.txt" TempT3
   #touch $(echo $tempM).fa   
  
let i=i+1

done 



#
echo "All clear, have a nice day"

#clean up

rm All_exomes_Test All_exomes_ID sourceL Lfasta.fa 

mkdir processed_fasta
mv *fasta cd processed_fasta
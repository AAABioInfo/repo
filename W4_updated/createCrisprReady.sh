#!/bin/bash
# input 1: $1= motif_list.txt $2=All_exomes $3=$All_exomes_Test
# Output: Output: dir(species_name)
                   # files:
                   #  species_name_motif_count_All.txt
                   #   species_name_count_TOP3
                   #   species_topMotif1_count.fa
                   #   species_topMotif2_count.fa
                   #   species_topMotif3_count.fa
#Creater Aldo Amaya
#Purpose: Identify the 3 highest occurring motifs in each exome inside the exomes folder. 

            # this will link up any desired list 
#TargerSeq=$(echo $2)
cat $1 >sourceL
TargetSeq=$(echo $2)
cat $3>All_exomes_ID

#used to test script
#cat motif_list.txt>sourceL                # this will link up any desired list 
#TargetSeq=$(echo All_exomes/)
#cat All_exomes_Test>All_exomes_ID

cp sourceL $TargetSeq
cp All_exomes_ID $TargetSeq
rm sourceL All_exomes_ID
cd $TargetSeq

lenM=$(sort sourceL| uniq|wc -l)
let i=1 # when using ls it counts All_exomes_Test 
lenExo=$(sort All_exomes_ID| uniq|wc -l)

# Go through each exoname and create a directory and store high motifs 
while [ $i -le $lenExo ] # go through list of all samples
do 
# from our sublist we will create a folder to contain all files top exo Motifs_count,i.e. fox_AAAA_103.fafsa,fox_AAAA_102.fafsa, fox_AAAA_101.fafsa  
  tempF=$(head -n $i All_exomes_ID| tail -n 1)  
  mkdir $tempF
  cat "$tempF.fasta" >Lfasta.fa
  #set params for new loop to grep each motif to fasta file
  ii=1 
  #touch "$tnExoname _motif_count_All.txt"
  while [ $ii -le $lenM ]
  do
  tempM=$(head -n $ii sourceL| tail -n 1)
  countN=$(grep $(echo $tempM) Lfasta.fa|wc -l) 
  echo $tempF $tempM $countN >>"${tempF}_motif_count_All.txt" 
  let ii=ii+1  
  done 

#ls -lt| head 
#copy each of these list into exoname folder
  cat "${tempF}_motif_count_All.txt" | sort -k3nr| head -n 3 > "${tempF}_motif_count_TOP3.txt"
  cp "${tempF}_motif_count_TOP3.txt" $tempF
  cp "${tempF}_motif_count_All.txt" $tempF
  
  # sets a temporary list  to reprocess 
  awk '{print $2}' "${tempF}_motif_count_TOP3.txt" > TempT3
  iii=1 # this will be referencing top list 
  lenTT3=$(cat TempT3|wc -l)   #lenTT3=3 
    while [ $iii -le $lenTT3 ]
    do
      tempM3=$(head -n $iii TempT3| tail -n 1)
      countN3=$(grep $(echo $tempM3) Lfasta.fa|wc -l|sed 's/ //g')
      grep -B1 $(echo $tempM3) Lfasta.fa> "${tempF}_${tempM3}_${countN3}.fasta"
      cp "${tempF}_${tempM3}_${countN3}.fasta" $tempF
      rm "${tempF}_${tempM3}_${countN3}.fasta"
      let iii=iii+1  
    done 
    

  #Clean up after each motif 
  rm "${tempF}_motif_count_TOP3.txt" "${tempF}_motif_count_All.txt" TempT3
   #touch $(echo $tempM).fa   
  
let i=i+1

done 

#
clear 

echo "All clear, have a nice day"

#clean up

rm sourceL Lfasta.fa 

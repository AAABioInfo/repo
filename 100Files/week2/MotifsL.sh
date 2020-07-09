#!/bin/bash

#This code was 201RBIF-100-1DL Week 2 assignment 
#This code will open a list of different motifs and create indidual fasta files and a motif count list.
#This script requires 2 inputs ./MotifsL.sh (list of motifs=$1) (.fa file=$2) 



#Created by Aldo A. Amaya
#Date: 04/23/2020




#####    Pre-set variables   #####


### you can use any two files, identifier and reference ###
cat $1 >sourceL                # this will link up any desired list 
cat $2 >source.fa               # this will be our fasta file (current limit to one file)

#let lenM=($3)                     #user will identify motifs

lenM=$(uniq sourceL|wc -l)  #this will count how many uniq motifs or list items are available, make upper bound
let i=1 #echo $i this will serve as our counter
mkdir motifs # this will create a directory where you will move all fasta files after
touch motif_count.txt   #this will create the list we need 



#####    Prepare fasrta file for processing    #####


# you must first leaniarized fasta file before conducting grep, this will help with accurate grep
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < source.fa>Lfasta.fa




#####    Process data/ create motif.fa and generate counts    #####

while [ $i -le $lenM ] # this createa a loop to  condition that you must go through the whole motifs list
do
  tempM=$(head -n $i sourceL| tail -n 1)   # This will roll through motifs list, while only looking at tail
  touch $(echo $tempM).fa                    # This will create the fasta file where we will store all matches per motif from refence fasta file
  
  countN=$(grep $(echo $tempM) Lfasta.fa|wc -l) #This will determine how many genes had the motif present 
  grep $(echo $tempM) Lfasta.fa> $(echo $tempM).fa  #This will call our temporary motif and reference to our linear fasta and rederict outout to motif.fa
  echo $tempM $countN >>motif_count.txt  #We will then call our motifs and counts and store into a motif count list 
  let i=i+1  #increases the counter by 1
done


#####    Show results for counts    #####

echo $'Operation succesful\nYour Results are below\n Motif   Present on these genes' #this is just for display that program went well 
cat motif_count.txt  #show results



#####    Move files    #####

mv [ATCG]*fa cd motifs   #This will look at all while cards and move the fa files that contain ATCG to motif dir

 
 
 

  


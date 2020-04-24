# This will actually be a readme.md


# The reason for this script is to outline the function of MotifsL.sh

# To Run you must input 2 files

# ie.     ./MotifsL.sh (list of motifs=$1) (.fa file=$2)



##### below are the details for W2 assingment. ####
# Species: bacteria, R.bifella.
# Desired result: identify 10 motifs associated with radiation resistancez



#Recipe for sucess:


# 1.) SSH into the AWS server and create a folder called week2.----- Done-----

#a.)ssh amaya@ec2-18-224-15-142.us-east-2.compute.amazonaws.com       Use credentials to enter

#b.)mkdir week2                                                       This will create a directory to build 




# 2.) Copy the motif list and the genome of the R.bifella from ----- Done-----
#      /home/rbif/week2_assignment/

#a.)cd /home/rbif/week2_assignment/              This comand was used to access the fafta data
#b.1)scp interesting_motifs.txt /home/amaya/repo
#b.2)scp r_bifella.fasta /home/amaya/repo
#b.3) git add interesting_motifs.txt r_bifella.fasta
#b.4) git commit -m "add true files for W2-assigntmet"
#b.5) git push


#Since we are on the server you need to utilize scp to copy;/
#                                        no ip needed siince on the same file


#Notes about file: interesting_motifs.txt : 10 motifs



#Notes about file:r_bifella.fast    : 500 genes









# 3.)Create a single bash script that does the following:-----Not Done-----

# fist set your permisions: chmod 755 chmod 755 MotifsL.sh
# MotifsL.sh 
  # $0 MotifsL.sh
  # $1 interesting_motifs.tx
  # $2 r_bifella.fasta


# Start of bash Script

#Print out the number of occurrences... output to a file called motif_count.txt----- Done-----

#Create a fasta file for each motif (we know that there are 10 in total); /----- Done-----
#Must contain all genes and their corresponding sequences that have that motif./----- Done-----

#Each file should be named after the motif (ie ATG.txt) and outputted to a new dir called motifs -----Not Done-----

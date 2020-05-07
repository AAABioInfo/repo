#This script will ensure we complete all objectives.

ssh amaya@ec2-18-224-15-142.us-east-2.compute.amazonaws.com
scp copuExomes.sh aaamaya@ec2-18-224-15-142.us-east-2.compute.amazonaws.com/home/amaya/week4



copyExomes.sh     ./copyExomes.sh (input= clinical_data.txt) 
                  Output cd/exomes/Produces Sequenced_20r30mm.txt
                  
createCrisprReady.sh
                  ./copyExomes.sh (input= motif_list.txt)     
                  Output: dir(species_name)
                    files:
                      species_name_motif_count_All.txt
                      species_name_count_TOP3
                      species_topMotif1_count.fa
                      species_topMotif2_count.fa
                      species_topMotif3_count.fa
                          
identifyCrisprSite.sh
editGenome.sh
exomeReport.py





#Inputs:
/home/rbif/week4/clinical_data.txt   -copied to personal AWS and Local machine- 
/home/rbif/week4/motif_list.txt      -copied to personal AWS and Local machine-
/home/rbif/week4/exomes/*.fasta      -copied to personal AWS and Local machine-


Steps:

1. Create a directory called week4 in your home directory.        ---Done---
    mkdir week4

2. 
  i. Create a bash script called copyExomes.sh        ---Done---

  ii. Read in the clinical data file and identify the samples that ---Done---
      a. have a diameter between 20 and 30 mm long (inclusive) ---Done---
      b. genomes sequenced. ---Done---
      c.Copy the identified exomes using the sample code names to a new directory called exomes. ---Done---
   

3. 
  i. Create a bash script called createCrisprReady.sh         ---Done---

  ii. Using the motif_list.txt file, 
      a. identify the 3 highest occurring motifs in each exome inside the exomes folder. ---Done---
      b. Output the headers and corresponding sequences to a new file called exomename_topmotifs.fasta.  ---Done---
          **This means to only select for the genes that have at least one of the three top motifs and output to a new                file
      c. Rename the headers of the exomename_topmotifs.fasta files to include the motif sequences and their count.     
                                      (Example: >gene1  ===> >gene1 ATGC_1 AAAA_38) ---Done---

          Note: The sequence does not need to contain all 3 motifs, it just needs to have at least 1 of them. 
          “exomename_topmotifs.fasta” should contain the code name of the sample. For example, fox_topmotifs.fasta.

4. 
  i. Create a bash script called identifyCrisprSite.sh         ---Done---

  ii. For each gene inside the exomename_topmotifs.fasta files, 
          a. identify a suitable CRISPR site. ---Done---
          b. Find sequences that contain “NGG”, where “N” can be any base, that has at least 20 basepairs upstream. ---Done---

                  Example of upstream: ATGAACGTCTGTAAGAACTGCGGATCTGTCA 
                  (Everything left of CGG is upstream of the DNA) 
          c. Output suitable candidates (headers and sequences) to a new file called exomenames_precrispr.fasta ---Done---
          
5. 

    i. Create a bash script called editGenome.sh         ---Done---
    
    ii. Using those files, write a script that will insert the letter A right before the NGG site. 
          b. Output to a new file called exomename_postcrispr.fasta 

6. 
    i. Create a python script called exomeReport.py         ---Done---
    
    ii. Create a single report that summarizes the findings. 
        a. It should be a text file that lists 
            i. the name of the discoverer of the organism, 
            ii.the diameter, 
            iii. the code name and the environment it came from. 
            iv. include how many genes the organism has in common with all of the other organisms, 
            v. and which genes are unique to that organism.
            
      Note: 
          -You are only comparing that organism against the other organisms of the same cohort (20-30mm in diameter). 
          -You are summarizing the motif-containing exomes with just the CRISPR-ready genes, 
                - using the exomenames_precrispr.fasta 
                    (or exomenames_postcrispr.fasta since the headers should be identical between those two files)
          -There may not be any genes in common in that entire cohort. 


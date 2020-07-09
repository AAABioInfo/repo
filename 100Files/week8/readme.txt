This script was originally crafted in order to look up a gene in mygene and find 
the ensembl name. 

Witht that name, it will return the Fasta and do the following:
    1.)Create a fasta file the format Gene_species.fasta
    2.) Determine the longest RF
    3.) Convert to Amino Acids using Biopython
    4.) Create a header
    5.) Append header and Amino Acid sequence of longest ORF
    6.) Using that orthologs from Ensembl  create a list of orthologs and make a file
    
    To run, use the following code
    
    python Query.py -g "gene of interest" -s "species" -f "field for mygene"
    
    for this problem, our script was excecuted as follows"
    
    python Query.py -g MC1R -s human, -f ensembl.gene

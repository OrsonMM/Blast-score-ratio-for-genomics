# Like Blast Score Ratio (BSR)

**Find_Genomic_Islands for Vibrio parahaemolyticus** 

These script find genomic islands (GIs) that may play a role in the emergence and pathogenesis of pandemic strains. Predicted Genes with prodigal in draft assembly genomes and return a Matrix.tab with Blast score ratio for each gene in your draft genome.

The concept of the BSR you can review these paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC545078/) 

Usage:


        $ python3 like_bsr.py dir_DraftGenomes Multifasta_for_DB/*fasta #cores name_output
        
  
dir_DraftGenomes = Directory with your draft genomes or complete genomes



Multifasta_for_DB = A file multifasta with genes of your interest. check - dir_multifasta 

# Dependences

You need install all programs in your environmetal path. See instructions for each case: 

1. python3 and libraries seaborn, numpy, dataframe.
2. BLAST ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
3. Prodigal https://github.com/hyattpd/Prodigal
       

# Find_Genomic_Island_7 for Vibrio parahaemolyticus

These script find genomic islands (GIs) that may play a role in the emergence and pathogenesis of pandemic strains. Predicted Genes with prodigal in draft assembly genomes and return a Matrix.tab with Blast score ratio for each gene in your draft genome.

Usage:


        python3 like_bsr.py dir_DraftGenomes Multifasta_for_DB/*fasta #cores name_output
        
  
dir_DraftGenomes = Directory with your draft genomes or complete genomes



Multifasta_for_DB = A file multifasta with genes of your interest. check - dir_multifasta 


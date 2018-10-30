#usr/bin/python

# Author B10Comp - Orson Mestanza 
# These script find genes in draft assembly genomes and return a Matrix with Blast score ratio for each gene in your multifasta
# You need Python2
# You need BLAST install in your environment 
# You need Prodigal install in your environment
# You need Pandas python

import sys
import os
import operator
import numpy as np
import seaborn as sns
import pandas as pd
import shlex, subprocess

import matplotlib.pyplot as plt

def funcion_makeblast(multifasta):

	comand_makeblast = "makeblastdb -in %s -dbtype nucl"%multifasta
	print("##########  Step1: Make a index of data_base  ##############")

	print(comand_makeblast)
	subprocess.call(comand_makeblast, shell=True)

	print("\n##########  Done data_base  ##############")

def do_prodigal(genomes):

	subprocess.call("mkdir prodigal_fasta", shell=True)
	print("##########  Step2: Prediction of CDS using prodigal  ##############")
	

	count = 0
	for i in os.listdir(genomes):
		if i.endswith('.fasta') or i.endswith('.fa') or i.endswith('.fna'):
			print(i)
			subprocess.check_output("prodigal -c -d prodigal_fasta/pro_%s -i %s%s -m -q"%(i,genomes,i), shell=True)
			count += 1
	print("\n##########  Done %d Predictions  ##############\n"%count)

	return "prodigal_fasta"

def do_blastn(Predictions,index_blast,cores):
	

	print("##########  Step3: Do blastn   ##############")

	count_blast = 0
	for i in os.listdir(Predictions):
		if i.endswith('.fasta'):
			new_name = i.split(".fasta")[0]
			print("blastn -db %s -query %s/%s -outfmt 6 -num_alignments 5 -perc_identity 60 -qcov_hsp_perc 75 -num_threads %d -out %s/%s.txt"%(index_blast,Predictions,i,int(cores),Predictions,new_name))
			subprocess.call("blastn -db %s -query %s/%s -outfmt 6 -num_alignments 5 -perc_identity 60 -qcov_hsp_perc 75 -num_threads %d -out %s/%s.txt"%(index_blast,Predictions,i,int(cores),Predictions,new_name), shell=True)
			count_blast += 1

	print("\n##########  Done %d blastn  ##############\n"%count_blast)

def extract_bsr_reference(data_base,cores):
	
	order_ref = [] # Order of genes in reference
	bsr_ref = []   # blast score ratio of reference 
	length_ref = [] 

	file = subprocess.check_output("blastn -db %s -query %s -outfmt 6 -num_alignments 5 -perc_identity 100 -qcov_hsp_perc 100 -num_threads %d"%(data_base,data_base,int(cores)),shell=True).decode("utf-8")

	bsr = file.split("\n")

	for i in range(len(bsr)):
		new_list = bsr[i].split("\t")
		if len(new_list) > 2:
			if new_list[0] == new_list[1] and int(new_list[4]) == 0:
				order_ref.append(new_list[0])
				bsr_ref.append(int(new_list[11]))
				length_ref.append(int(new_list[3]))

	return order_ref, bsr_ref, length_ref

def do_matrix_bsr(blastn_files,order_gene,bsr_gene,len_gene):


	#df1 = pd.DataFrame(columns=order_gene)

	frames = []

	print("##########  Step3: Select best hit and do matrix_bsr  ##############")

	for i in os.listdir(blastn_files):
		if i.endswith('.txt'):			
			index = i.split(".txt")[0]
			nombres = index
			index = {}
			with open("%s/%s"%(blastn_files,i)) as f:
				for line in f:
					line = line.split("\t")
					if not line[1] in index:
						if float(line[3])/len_gene[order_gene.index(line[1])]*100 > 70:  # Si el alineamiento entre la query y el subjet es mas del 60%
							new = line[1]
							index[new] = float(line[11])/float(bsr_gene[order_gene.index(line[1])])
					else:
						if float(line[3])/len_gene[order_gene.index(line[1])]*100 > 70:
							if index[line[1]] < float(line[11])/float(bsr_gene[order_gene.index(line[1])]):
								index[line[1]] = float(line[11])/float(bsr_gene[order_gene.index(line[1])])

				
				df = pd.DataFrame(index, index=[nombres], columns=order_gene)
				print(df) 			
				frames.append(df) 
				index = {}

	results = pd.concat(frames)
	results = results.replace(np.nan, 0)
	
	results.to_csv("result_Dataframe.txt", sep='\t', encoding='utf-8')
	#print(results)
	#sns.heatmap(results, annot=False, xticklabels=1, linewidths=1, linecolor='black')
	#plt.show()


def main():
		
	print("Usage:")
	print("       python like_bsr.py dir_DraftGenomes Multifasta_for_DB\n\n")  

	file0 = sys.argv[0]
	genomes = sys.argv[1]
	data_base = sys.argv[2]
	cores = sys.argv[3]

	
	funcion_makeblast(data_base)           # These function do makeblastdb with multifasta input 

	nucle_files = do_prodigal(genomes)     # These function predicted CDS using prodigal and return the new dir with cds. 

	do_blastn(nucle_files,data_base,cores) # These functoin make a blastn cds predicted with prodigal vs your genes in multifasta.

	order_gen, bsr_gen, len_gen = extract_bsr_reference(data_base,cores) # These function return two arrays. 1. Is array with order of genes of reference. 2 is array with blast score ratio for each gene 3. Is array with the length of each gene. 

	do_matrix_bsr(nucle_files,order_gen,bsr_gen,len_gen)

if __name__ == "__main__":
	main()


#Some problems that you can solve in the future 
# 1. For example the genes is not in the same position to reference.
# 2. Is necessary plot a heatmap with these data in near future.

	

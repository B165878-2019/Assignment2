#!/usr/bin/python3
import subprocess
import os



##### Entrez utilities
taxonomy = input("which taxonomic group would you like to inspect? ")
protname = input("what is the name of your desired protein family? ")
sequences_fa = open(f"{taxonomy}.fa").read()
with open(f"{taxonomy}.fa", "w") as f:
 query_term = (f"{taxonomy}[ORGN] AND {protname}[PROT] NOT PARTIAL")
 esearch = subprocess.Popen(
	['esearch', '-db', 'protein', '-query', query_term],
	stdout=subprocess.PIPE)
 efetch = subprocess.Popen(
	['efetch', '-db', 'protein', '-format', 'fasta'],
	stdin=esearch.stdout, stdout=f)


##### Clustal
clustalo = subprocess.Popen(['clustalo', '-i', sequences_fa], stdout=subprocess.PIPE) #doesnt open file correctly
subprocess.Popen(['cons', '-outseq', 'consensus.fa'], stdin=clustalo.stdout)




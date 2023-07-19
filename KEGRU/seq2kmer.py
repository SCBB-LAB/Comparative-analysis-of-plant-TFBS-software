#!usr/bin/python
import sys,os
infile=[x.strip() for x in open(sys.argv[1]).readlines()]#Input sequence or reads file

for line in infile:
	if line.startswith(">"):
		seq_id = line.strip()
	else:		
		def seq2kmer(seq, k):

			kmer = [seq[x:x+k] for x in range(0,len(seq)+1-k,2)] # stride size 2
			
			kmers = " ".join(kmer)			
			print(kmers)
			return kmers
		seq2kmer(line,5)

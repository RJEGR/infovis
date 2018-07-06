#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# @uthor:      rgomez@cicese.edu.mx with help from Ramirez-Rafael J.
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#
# This file contains code that helps to wragling fasta/fastq files.
# use python 3.6+ to improve this code
# to use, three argument needed: sequence file (fasta/fastq format) output_name_file sequence format(fastq/fastq)
# Ex: length_seqs.py input.fasta output.lengths.csv fasta
#------------------------------
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import sys

fastafile = sys.argv[1]
out = open(sys.argv[2], 'w')
frt = sys.argv[3]
L= [[record.id,len(record)] for record in SeqIO.parse(fastafile, frt, generic_dna)]

# write the file by pipe the replace and upper function through your data content
for i in range(len(L)):
    out.write(str(L[i][0])+ " " +str(L[i][1])+'\n')
out.close()

print("Length analysis done ...")
	#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# @uthor:      rgomez@cicese.edu.mx with help from Ramirez-Rafael J.
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007

"""
    This file contains code than helps to wragling fasta/fastq files. \
    The code was developed in python 3.6+. \
    By default a set of 100 sequence is sampled and output_name_file gonna be the same than input but using the "subset" prefix, \
    This current version output an standard-format in FASTA''',
    to use, two argument needed: sequence file (fasta/fastq format) and sequence format(fastq/fastq).\
     Ex: subsample.py input.fasta fasta 
"""

from Bio import SeqIO
from random import sample
import sys



file = sys.argv[1]
#out = open(sys.argv[2], 'w')
out = open(file+'.subset', 'w')
frt = sys.argv[2]
sys.stdout = out
orig_stdout = sys.stdout


with open(file) as f:
    seqs = SeqIO.parse(f, frt)
    subset = ((seq.name, seq.seq) for seq in  sample(list(seqs),100))
    for i in subset:
        print(">{}\n{}".format(*i))

sys.stdout = orig_stdout
out.close()

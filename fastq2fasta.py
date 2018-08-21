#!/usr/bin/python3.6
# -*- coding: utf-8 -*-
#
#------------------------------
# @uthor:      rgomez@cicese.edu.mx with help from Ramirez-Rafael J.
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007

"""
    This file contains code than helps to convert fastq to fasta files. \
    The code was developed in python 3.6+. \
    to use, one argument needed: sequence file (fastq format); by default the file name is used to output (using sufix fasta)
     Ex: fastq2fasta.py input.fastq 
"""

from Bio import SeqIO
import sys

file = sys.argv[1]


SeqIO.convert(file, 'fastq', file+'.fasta', 'fasta')


print("Convertion from", file, "file was done")

exit


# python3 -c "from Bio import SeqIO;SeqIO.convert('T0-R_S1_L001_R1_001.fastq.gz.P.qtrim', 'fastq', 'T0-R_S1_L001_R1_001.fastq.gz.P.qtrim.fa', 'fasta')"

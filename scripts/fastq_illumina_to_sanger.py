#/usr/bin/env python
import sys
from Bio import SeqIO

if len(sys.argv) != 1:
    sys.exit("ERROR! Pipe me your fastq file and I'll pipe you another one!")
    
SeqIO.convert( sys.stdin , "fastq-illumina", sys.stdout, "fastq-sanger")

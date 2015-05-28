#!/usr/bin/env python

#Takes a single FASTQ file and splits to .fasta + .qual files
import sys
from Bio import SeqIO

if len(sys.argv) != 1: 
    print "This is supposed to be piped. cat your fasta, or ungzip it to stdout."
    sys.exit()

SeqIO.convert( sys.stdin , "fastq", sys.stdout, "qual")

#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
import sys
import gzip

if len(sys.argv) != 3:
    print "ERROR: Incorrect number of files"
    print "Usage:" + sys.argv[0] + " file.fasta file.qual"
    print "Fastq file will be written to stdout" 
    sys.exit()


fasta_in   = open( sys.argv[1] )
qual_in    = open( sys.argv[2] )

record_iterator = PairedFastaQualIterator( fasta_in , qual_in )
SeqIO.write(record_iterator, sys.stdout, "fastq")

 #!/bin/python
from Bio import SeqIO # Import Biopython's file reader
import sys            # Import sys to handle with arguments

# Usage: pairer.py \
#           fileIn_1
#           fileIn_2
#           fileIn_1  # Again
#           fileIn_2  # Again
#           fileOut_1
#           fileOut_2
#           fileOut_3
#           fileOut_4

if len(sys.argv) != 9:
    print "ERROR: Incorrect number of files"
    quit()

input_1      = sys.argv[1] # Pipe from _1.fastq.gz
input_2      = sys.argv[2] # Pipe from _2.fastq.gz
input_1_bis  = sys.argv[3] # Pipe from _1.fastq.gz (Again)
input_2_bis  = sys.argv[4] # Pipe from _2.fastq.gz (Again)
output_1     = sys.argv[5] # Pipe output _1
output_2     = sys.argv[6] # Pipe output _2
output_3     = sys.argv[7] # Pipe output _3
output_4     = sys.argv[8] # Pipe output _4

# readIdentifiers: opens a fastq file and returns a set made of all the identifiers inside
def readIdentifiers( fileIn ):
    identifiers = set()                                  # Initialize set
    handle_input = open( fileIn , "rU" )                 # Open connection
    for record in SeqIO.parse( handle_input , "fastq" ): # Iterate over the registries of the FASTQ
        identifiers.add(record.id)
    handle_input.close()
    return identifiers

index_1 = readIdentifiers( input_1 ) # Read identifiers from the first input
index_2 = readIdentifiers( input_2 ) # Read identifiers from the second input

index_paired = index_1 & index_2     # Find the reads that are still paired

add_to_3 = index_1 - index_paired    # Make a set with the unpaired forward reads
add_to_4 = index_2 - index_paired    # Make a set with the unpaired reverse reads

del [ index_1 , index_2 ]            # Free memory

# Writes the fastqs associated to the paired and unpaired reads
def writeFastqs ( index_paired , add_to_unpaired , input_file , output_paired , output_unpaired ):
    handle_input           = open( input_file      , "rU" ) # In (read, universal)
    handle_output_paired   = open( output_paired   , "w"  ) # Out paired (write)
    handle_output_unpaired = open( output_unpaired , "w"  ) # Out unpaired (write)
    for record in SeqIO.parse( handle_input , "fastq"):            # Iterate over the entries of the FASTQ file
        if record.id in add_to_unpaired:                           # Check belonging to the smaller set
            SeqIO.write(record, handle_output_unpaired , "fastq" )
        else:
            SeqIO.write(record, handle_output_paired   , "fastq" )
    handle_input.close()            # Close connections
    handle_output_paired.close()
    handle_output_unpaired.close()

writeFastqs( index_paired , add_to_3 , input_1_bis, output_1 , output_3 )
del add_to_3

writeFastqs( index_paired , add_to_4 , input_2_bis, output_2 , output_4 )
del add_to_4

del index_paired
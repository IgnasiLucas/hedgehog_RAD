#!/bin/bash
#
#2016-11-07
#
#Demultiplexing raw reads using Sabre (https://github.com/najoshi/sabre)
#
#Files representing raw paired-end reads: runII.1.fasta.gz(forward), runII.2.fasta.gz(reverse)
#present in the directory together with .txt file barcode

if [ ! -e Er51_436_1.fastq ]; then
   sabre pe -m 1 -c -f runII.1.fasta.gz -r runII.2.fasta.gz -b barcode.txt -u unknown_barcode1.fasta -w unknown_barcode2.fasta
fi

#Options:
#sabre (pe) pair-end reads,
#(-m 1) allows for one mismatch in the barcode sequence,
#(-c) remove barcodes from both files,
#(-f) specifying forward input file,
#(-r) specifying reverse input file, (-b) barcode file,
#(-u) output file contains unknown forward sequences,
#(-w) output file contains unknown reverse sequences.
#
#One allowed mismatch in the barcode sequences derived
#from python script hamm.py assessing minimum distance
#between words(five letters barcodes).
#Script finished run with minimum distance between words = 3.

# Here we count the frequency of the first 5-letters words in the unknown_barcode1.fasta file
if [ ! -e counts ]; then
   gawk '(NR % 4 == 2){F[substr($1,1,5)]++}END{for (f in F) print f "\t" F[f]}' unknown_barcode1.fasta | \
   sort -nrk 2,2 > counts
fi

# Now, I take the 30 most frequent words that were not identified, and measure their
# Hamming distance to the real barcodes.

if [ ! -e comparison ]; then
   if [ ! -e bad_words ]; then
      head -n 30 counts | cut -f 1 > bad_words
   fi

   if [ ! -e good_and_bad_words ]; then
      # There are 25 real codewords
      cut -f 1 -d ' ' barcode.txt > good_words
      cat good_words bad_words > good_and_bad_words
   fi

   python ./hamm.py -o comparison good_and_bad_words
   rm bad_words good_words good_and_bad_words
fi

# The file 'comparison' is a triangular matrix. There is a header row, and 25 + 30 rows.
# If I remove the first 26 rows, and the last 30 columns, I will get a 30 * 25 matrix
# with the comparison that we need. I will add the header as well.

if [ ! -e distances ]; then
   head -n  1 comparison | cut -f 1-26 > z1
   tail -n 30 comparison | cut -f 1-26 > z2
   cat z1 z2 > distances
   rm z1 z2
fi

# As expected, all distances are higher than 1. There are only two or three cases where
# the non-recognized word is at distance 2 from a real codeword and at larger distances
# from all other codewords. We can extract them from the file 'distances' like this:

if [ ! -e recovered ]; then
   gawk -v FS='\t' '(NR == 1){
      for (i=1; i<=NF; i++) {
         CODE[i] = $(i+1)
      }
fi

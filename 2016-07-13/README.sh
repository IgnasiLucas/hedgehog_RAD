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

# Now, I take the N most frequent words that were not identified, and measure their
# Hamming distance to the real barcodes.

N=100
if [ ! -e comparison ]; then
   if [ ! -e bad_words ]; then
      head -n $N counts | cut -f 1 > bad_words
   fi

   if [ ! -e good_and_bad_words ]; then
      # There are 25 real codewords
      cut -f 1 -d ' ' barcode.txt > good_words
      cat good_words bad_words > good_and_bad_words
   fi

   python ./hamm.py -o comparison good_and_bad_words
   rm bad_words good_words good_and_bad_words
fi

# The file 'comparison' is a triangular matrix. There is a header row, and 25 + N rows.
# If I remove the first 26 rows, and the last N columns, I will get a N * 25 matrix
# with the comparison that we need. I will add the header as well.

if [ ! -e distances ]; then
   head -n  1 comparison | cut -f 1-26 > z1
   tail -n $N comparison | cut -f 1-26 > z2
   cat z1 z2 > distances
   rm z1 z2
fi

# As expected, all distances are higher than 1. There are only a few cases where
# the non-recognized word is at distance 2 from a real codeword and at larger distances
# from all other codewords. We can extract them from the file 'distances' like this:

if [ ! -e closest ]; then
   gawk -v FS='\t' '(NR == 1){
      for (i=2; i<=NF; i++) {
         CODE[i] = $i
      }
   }(NR > 1){
      DIST2 = ""
      HOWMANY = 0
      for (i=2; i<=NF; i++) {
         if ($i == 2) {
            HOWMANY++
            DIST2 = CODE[i] DIST2
         }
      }
      if (HOWMANY == 1) {
         print $1 "\t" DIST2
      }
   }' distances > closest
fi

# Below, I paste the log from sabre, and I add manually the number of
# reads that we could recover from those that were not identified. To
# recover them, we should create another barcode file with the new words
# that we trust, and run sabre again on the unclassified reads, without
# allowing any mismatch.
#
#   FastQ records for barcode TAGCA: 3176888 (1588444 pairs) + 77492
#   FastQ records for barcode TACGT: 5522138 (2761069 pairs)
#   FastQ records for barcode TAATG: 5134546 (2567273 pairs)
#   FastQ records for barcode GTTGT: 4006542 (2003271 pairs)
#   FastQ records for barcode GTGTG: 7282104 (3641052 pairs)
#   FastQ records for barcode GTCAC: 6310898 (3155449 pairs) + 26639
#   FastQ records for barcode GTACA: 5576506 (2788253 pairs) +  6009
#   FastQ records for barcode GGTTC: 1733294 (866647 pairs)  + 23845
#   FastQ records for barcode GGGGA: 2113556 (1056778 pairs)
#   FastQ records for barcode GGCCT: 8847316 (4423658 pairs)
#   FastQ records for barcode GGAAG: 10170846 (5085423 pairs)+  6018
#   FastQ records for barcode GCTAA: 4293996 (2146998 pairs) + 82715
#   FastQ records for barcode GCGCC: 5311396 (2655698 pairs) + 23065
#   FastQ records for barcode GCCGG: 311268 (155634 pairs)
#   FastQ records for barcode GCATT: 558808 (279404 pairs)   + 51046
#   FastQ records for barcode GATCG: 3206390 (1603195 pairs) +  5984
#   FastQ records for barcode GAGAT: 477280 (238640 pairs)   + 37361
#   FastQ records for barcode GACTA: 1794746 (897373 pairs)  +  6458
#   FastQ records for barcode GAAGC: 2957740 (1478870 pairs)
#   FastQ records for barcode CTTCC: 7307936 (3653968 pairs)
#   FastQ records for barcode CTCTT: 6518346 (3259173 pairs)
#   FastQ records for barcode CTAGG: 7054932 (3527466 pairs)
#   FastQ records for barcode CGTAT: 7334390 (3667195 pairs)
#   FastQ records for barcode CGGCG: 3569452 (1784726 pairs)
#   FastQ records for barcode CGCGC: 9164044 (4582022 pairs)

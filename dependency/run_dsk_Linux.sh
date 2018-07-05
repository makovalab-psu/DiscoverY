#!/bin/sh

#uses dsk2.2.0
#run as ./run_dsk.sh FASTQ_file_to_be_kmerized kmer_size

if [ $# -ne 2 ]
then
	echo "Usage: $0 <FASTQ_file_to_be_kmerized> <kmer_size>"
	exit 1
fi


R1_fsY_reads=$1
kmer_size=$2

echo "Counting k-mers from the following file : " $R1_fsY_reads

dsk-v2.2.0-bin-Linux/bin/dsk -file $R1_fsY_reads -abundance-min 0 -kmer-size $kmer_size -out R1_dsk -verbose 0

dsk-v2.2.0-bin-Linux/bin/dsk2ascii -file R1_dsk -out kmers_from_reads -verbose 0

rm R1_dsk.h5

rm -rf dsk_output
mkdir dsk_output
mv kmers_from_reads dsk_output/

echo "Counting finished! K-mer counts are in the folder './dsk_output'"

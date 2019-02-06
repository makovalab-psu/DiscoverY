from collections import defaultdict
from Bio import SeqIO
from pybloomfilter import BloomFilter
import sys


def exit_gracefully():
    sys.exit("^DiscoverY exited with error, please see the message above.^")


def reverse_complement(seq):
    complement_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq[::-1].translate(complement_map)


def make_bloom_from_kmer_abundance(ip_kmer_set, kmer_size, bf_size, bf_filename, ignore_rc=True):
    bf = BloomFilter(bf_size, 0.001, bf_filename)
    bf_count = 0
    for kmer in ip_kmer_set :
        bf.add(kmer)
        if not ignore_rc :
            bf.add(reverse_complement(kmer))
        bf_count+=1
    return bf


def test_valid_kmer_format(sample_line, kmer_size):
    sample_kmer_abun = sample_line.split(" ")
    if len(sample_kmer_abun) != 2:
        print("Format of kmer file is incorrect. Need a (kmer, count) pair separated by whitespace.")
        exit_gracefully()
    if len(sample_kmer_abun[0]) != kmer_size:
        print("Format of kmer file does not match kmer size specified for DiscoverY")
        exit_gracefully()
    try:
        int_test_line_kmer_abundance = int(sample_kmer_abun[1])
    except ValueError:
        print("Format of kmer file is incorrect. K-mer abundance in kmers_from_reads is not an int")
        exit_gracefully()


def make_dict_from_kmer_abundance (ip_file, kmer_size):
    kmer_dicts = defaultdict(int)
    with open(ip_file,'r') as file_handle1:
        first_line = file_handle1.readline()
        test_valid_kmer_format(first_line, kmer_size)
        for line in file_handle1:
            current_abundance = int(line.split(' ')[1])
            kmer_dicts[line[:kmer_size]] = current_abundance
    return kmer_dicts



# make a FASTA file with annotated contigs
def write_annotated_contigs_to_fasta (annotated_records) :
    with open("proportion_annotated_contigs.fasta","w") as output_handle :
        SeqIO.write(annotated_records, output_handle, "fasta")
    return 1


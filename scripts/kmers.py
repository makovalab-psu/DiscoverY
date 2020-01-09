from collections import defaultdict
import sys


def exit_gracefully():
    sys.exit("^DiscoverY exited with error, please see the message above.^")

try:
	# for python 3, where maketrans() is a function
	from string import maketrans
	complement_map = maketrans("ACGTNacgtn", "TGCANtgcan")
	def reverse_complement(seq):
	    return seq[::-1].translate(complement_map)
except ImportError:
	# for python 3, where maketrans() is a class function
	def reverse_complement(seq):
	    complement_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")
	    return seq[::-1].translate(complement_map)


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
        file_handle1.seek(0)
        for line in file_handle1:
            current_abundance = int(line.split(' ')[1])
            kmer_dicts[line[:kmer_size]] = current_abundance
    return kmer_dicts





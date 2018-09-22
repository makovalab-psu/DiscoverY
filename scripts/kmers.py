from collections import defaultdict
from Bio import SeqIO
from pybloomfilter import BloomFilter
#from string import maketrans
import sys


def exit_gracefully():
    sys.exit("^DiscoverY exited with error, please see the message above.^")


def reverse_complement(seq):
    complement_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq[::-1].translate(complement_map)


def kmerize(ip_string, kmer_size, step_size=1, offset=0):
    """
    function that kmerizes an input_string and returns a list of kmers
    """
    
    return [str(ip_string[i:i + kmer_size]) for i in range(offset, len(ip_string) - kmer_size + 1, step_size)
            if str(ip_string[i:i + kmer_size]) != "N"*kmer_size]
    


def make_bloom_from_kmer_abundance(ip_kmer_set, kmer_size, bf_size, bf_filename, ignore_rc=True):
    bf = BloomFilter(bf_size, 0.001, bf_filename)
    bf_count = 0
    for kmer in ip_kmer_set :
        bf.add(kmer)
        if not ignore_rc :
            bf.add(reverse_complement(kmer))
        bf_count+=1
        #if bf_count % 10000 == 0 :
        #    print "Number of kmers add to BF so far is : ", bf_count
    #bf.update(ip_kmer_set)
    return bf


# this one is for pybloom_live
def legacy_make_bloom_from_kmer_abundance(ip_file, kmer_size, ignore_rc=True):
    """
    function that bloomifies a large file with kmers and count
    file looks like : ATGCT 101
    """
    op_bloom = pybloom_live.ScalableBloomFilter(mode=pybloom_live.ScalableBloomFilter.SMALL_SET_GROWTH)
    with open(ip_file, 'r') as kmers_fp:
        for line in kmers_fp:
            op_bloom.add(line[:kmer_size])
            if ignore_rc == False:
                op_bloom.add(reverse_complement(line[:kmer_size]))
    return op_bloom


def make_set_from_kmer_abundance(ip_file, kmer_size, ignore_rc=False):
    """
    function that settifies a large file with kmers and count
    file looks like : ATGCT 101
    """
    kmers = []
    with open(ip_file, 'r') as kmers_fp:
        for line in kmers_fp:
            kmers.append(line[:kmer_size])
            if ignore_rc == False :
                kmers.append(reverse_complement(line[:kmer_size]))
    return set(kmers)


def make_dict_from_kmer_abundance (ip_file, kmer_size, ignore_rc=False):
    kmer_dicts = defaultdict(int)
    with open(ip_file,'r') as file_handle1 :
        for line in file_handle1 :
            current_abundance = int(line.split(' ')[1])
            kmer_dicts[line[:kmer_size]] = current_abundance
            # kmer_dicts[reverse_complement(line[:kmer_size])] = current_abundance
    return kmer_dicts


# function that uses BioPython to convert a fasta file to a dict
def parse_ctgs_fasta (ip_fasta_file) :
    with open(ip_fasta_file, "r") as handle :
        contigs_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    return contigs_dict


def test_parse_ctgs () :
    ip_fasta_file = "test/contigs.fasta"
    contigs = parse_ctgs_fasta(ip_fasta_file)
    for k, _ in contigs.iteritems() :
        print(k, contigs[k].seq)


# make a FASTA file with annotated contigs
def write_annotated_contigs_to_fasta (annotated_records) :
    with open("proportion_annotated_contigs.fasta","w") as output_handle :
        SeqIO.write(annotated_records, output_handle, "fasta")
    return 1


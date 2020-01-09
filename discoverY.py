import argparse
from scripts import classify_ctgs, kmers


def main():
    """
    The following input files are required in ./data folder
    1. contigs.fasta : Male contigs fasta file
    2. kmers_from_reads : Kmers from male reads with their counts
    3. female.bloom or female.fasta or female_kmers
    """
    parser = argparse.ArgumentParser(description='DiscoverY selects Y-specific contigs from a male assembly')
    parser.add_argument('--female_bloom', help='Use if female BloomFilter provided (defaults to False)',
                        action='store_true', required=False)
    parser.add_argument('--female_kmers_set', help='Use if female kmers set provided (defaults to False)',
                        action='store_true', required=False)
    parser.add_argument('--kmer_size', help='Set kmer size (defaults to 25)', required=False)
    parser.add_argument("--female_bloom_capacity", type=int, help="If female BloomFilter is not provided, \
    This sets the capacity of the bloom filter. This should be greater than the number of distinct kmers \
    in female.fasta or female_kmers.")

    parser.add_argument('--mode', help="Set to run in femal_only mode or female+male mode. \
    In female_only mode, the proportion shared between each contig with a female reference is computed.\
    In female+male mode, both the proportion of k-mers shared with the female and the median depth of coverage of each male contig is measure", required=True)
    
    args = vars(parser.parse_args())

    k_size = 25
    bloom_filt = False
    female_kmers = False

    # set kmer_size from argument or using default here
    if not args['kmer_size']:
        k_size = 25
    else:
        try:
            k_size = int(args['kmer_size'])
        except ValueError:
            print("Error : kmer_size provided is not an integer")
            kmers.exit_gracefully()

    if args['female_bloom']:
        bloom_filt = True
    elif args['female_bloom_capacity']:
        try:
            bf_capacity = int(args['female_bloom_capacity'])
        except ValueError:
            print("Error : female_bloom_capacity provided is not an integer")
            kmers.exit_gracefully()
    else:
        bf_capacity = 3 * 1000 * 1000 * 1000

    print("Started DiscoverY")

    if args['female_kmers_set']:
        female_kmers = True

    mode = args["mode"]
    print("Mode", mode)
    # declare defaults
    print("Using default of k=25 and input folder='data'")
    print("Shortlisting Y-contigs")

    classify_ctgs.classify_ctgs(k_size, bloom_filt, bf_capacity, female_kmers, mode)

    print("DiscoverY completed successfully")

if __name__ == "__main__":
    main()


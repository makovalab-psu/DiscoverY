from __future__ import division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scripts import kmers
from pybloomfilter import BloomFilter
from collections import defaultdict
from numpy import median


def classify_ctgs(kmer_size, bf, fem_kmers):
    contigs_fasta_file = "data/contigs.fasta"
    reads_kmers = "data/kmers_from_male_reads"
    # output
    annotated_records = []

    # make Bloom Filter of kmers from female
    if bf:
        print("Opening Bloom Filter of k-mers from female")
        female_kmers_bf = BloomFilter.open("data/female.bloom")
        print("Done")
    else:
        print("Need to make Bloom Filter of k-mers from female")
        print("Reading female reference one record at a time and k-merizing each record...")
        if fem_kmers: # if female kmers file exist
            female_kmers_file = "data/female_kmers"
            female_kmers_set = kmers.make_set_from_kmer_abundance(female_kmers_file, kmer_size)
        else :
            female_reference_file = "data/female.fasta"
            all_female_kmers = []
            for record in SeqIO.parse(female_reference_file,"fasta"):
                curr_seq = record.seq
                to_kmerize_fwd = str(curr_seq).upper()
                length = len(to_kmerize_fwd)
                for i in range(0, length-kmer_size+1):
                    all_female_kmers.append(to_kmerize_fwd[i:i+kmer_size])

            print("All kmerizing done, now converting to a set")
            female_kmers_set = set(all_female_kmers)

        print("Here is a sample female kmer in set : ", next(iter(female_kmers_set)))
        print("Cardinality of the female kmer set is : ", len(female_kmers_set))

        print("Generating Bloom Filter of k-mers from female reference")
        female_kmers_bf_size = 3000000000
        male_and_absent_in_female_bf_filename = "data/female.bloom"
        female_kmers_bf = kmers.make_bloom_from_kmer_abundance\
            (female_kmers_set, kmer_size, female_kmers_bf_size,
             male_and_absent_in_female_bf_filename, ignore_rc=True)
        print("Done")

    # make a dict of all kmers from male
    print("Generating a dictionary from kmers in kmers_from_male_reads")
    kmer_abundance_dict_from_male = defaultdict(int)
    kmer_abundance_dict_from_male = kmers.make_dict_from_kmer_abundance(reads_kmers, kmer_size)

    # store all contigs in a dict
    print("Parsing contigs into a SeqIO dictionary...")
    all_male_ctgs = kmers.parse_ctgs_fasta(contigs_fasta_file)
    print("Done")

    ctg_count = 0
    # for each contig in contigs.fasta, calculate proportion of kmers that are present in female data set
    for contig in all_male_ctgs.keys():
        ctg_count += 1
        print("No. of contigs seen so far: ", ctg_count)
        print("Current contig ID is : ", contig)

        # convert the SeqRecord object's .seq into a string that needs to be kmerized
        to_kmerize_fwd = str(all_male_ctgs[contig].seq).upper()

        length = len(to_kmerize_fwd)
        reverse = kmers.reverse_complement(to_kmerize_fwd)
        count_of_male_kmers = 0
        count_of_male_kmers_not_shared_with_female = 0
        curr_kmer_abundances = []
        for i in range(0, length-kmer_size+1):
            kmer = to_kmerize_fwd[i:i+kmer_size]
            rev_kmer = reverse[i:i+kmer_size]
            count_of_male_kmers += 1
            # feature 1 : compute proportion of kmers shared with female
            if kmer not in female_kmers_bf and rev_kmer not in female_kmers_bf:
                    count_of_male_kmers_not_shared_with_female += 1
            # feature 2 : find median of kmer abundances from reads
            if kmer in kmer_abundance_dict_from_male :
                curr_kmer_abundances.append(kmer_abundance_dict_from_male[kmer])
            elif rev_kmer in kmer_abundance_dict_from_male :
                curr_kmer_abundances.append(kmer_abundance_dict_from_male[rev_kmer])
            else:
                continue
        if curr_kmer_abundances:
            curr_median = median(curr_kmer_abundances)
        else:
            curr_median = 0

        print("Total No. of k-mers from this contig: ", count_of_male_kmers)
        print("No. of k-mers not shared with female: ", count_of_male_kmers_not_shared_with_female)

        if count_of_male_kmers != 0 :
            proportion = count_of_male_kmers_not_shared_with_female / count_of_male_kmers
        else:
            proportion = -1

        print("Proportion is: ", proportion)
        print("Median is: ", curr_median)

        annotated_records.append(SeqRecord((all_male_ctgs[contig].seq), id=all_male_ctgs[contig].id, description=
        str(str(len(to_kmerize_fwd)) + " " + str(proportion) + " " + str(curr_median))))

    kmers.write_annotated_contigs_to_fasta(annotated_records)
    return 1

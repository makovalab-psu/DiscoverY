from __future__ import division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scripts import kmers
from pybloomfilter import BloomFilter
from collections import defaultdict
from numpy import median



def getbloomFilter(bf, fem_kmers, kmer_size):
    if bf:
        print("Opening Bloom Filter of k-mers from female")
        female_kmers_bf = BloomFilter.open("data/female.bloom")
        print("Done")
    else:
        print("Need to make Bloom Filter of k-mers from female")
        print("Reading female reference one record at a time and k-merizing each record...")
        bf_size = 3 * 1000 * 1000 * 1000
        bf_filename = "data/female.bloom"
        female_kmers_bf = BloomFilter(bf_size, .001, bf_filename)

        if fem_kmers: # if female kmers file exist
            female_kmers_file = "data/female_kmers"
            with open(female_kmers_file, 'r') as fm_kmers:
                #assumes kmers are uppercase 
                for line in fm_kmers:
                    female_kmers_bf.add(line[:kmer_size])
        else :
            female_reference_file = "data/female.fasta"
            n_kmers = "N"*kmer_size
            for record in SeqIO.parse(female_reference_file,"fasta"):
                to_kmerize_fwd = str(record.seq).upper()
                length = len(to_kmerize_fwd)
                for i in range(0, length-kmer_size+1):
                    female_kmer = to_kmerize_fwd[i:i+kmer_size]
                    if female_kmer != n_kmers:
                        female_kmers_bf.add(to_kmerize_fwd[i:i+kmer_size])

        print("Done creating bloom filter")
    return female_kmers_bf


def classify_fm_male_mode(kmer_size, female_kmers_bf):
    male_contigs_file = "data/male_contigs.fasta"
    reads_kmers = "data/kmers_from_male_reads"
    annotated_contigs = open("proportion_annotated_contigs.fasta", "w")
    
    print("Generating a dictionary from kmers in kmers_from_male_reads")
    kmer_abundance_dict_from_male = kmers.make_dict_from_kmer_abundance(reads_kmers, kmer_size)
    
    ctg_count = 0
    n_kmers = "N"*kmer_size
    #for each contig, calculate proportion of kmers that are present in female data set
    for record in SeqIO.parse(male_contigs_file, "fasta"):
        ctg_count += 1
        print("No. of contigs seen so far: ", ctg_count)
        print("Current contig ID is : ", record.id)
        to_kmerize_fwd = str(record.seq).upper()
        length = len(to_kmerize_fwd)
        reverse = kmers.reverse_complement(to_kmerize_fwd)
        count_of_male_kmers = 0
        count_of_male_kmers_not_shared_with_female = 0
        curr_kmer_abundances = []
        for i in range(0, length-kmer_size+1):
            kmer = to_kmerize_fwd[i:i+kmer_size]
            if kmer == n_kmers:
                continue
            rev_kmer = reverse[length-kmer_size-i:length-i]
            count_of_male_kmers += 1
            # feature 1 : compute proportion of kmers shared with female
            if kmer not in female_kmers_bf and rev_kmer not in female_kmers_bf:
                count_of_male_kmers_not_shared_with_female += 1
            # feature 2 : find median of kmer abundances from reads
            if kmer in kmer_abundance_dict_from_male:
                curr_kmer_abundances.append(kmer_abundance_dict_from_male[kmer])
            elif rev_kmer in kmer_abundance_dict_from_male :
                curr_kmer_abundances.append(kmer_abundance_dict_from_male[rev_kmer])
            else:
                continue
        if curr_kmer_abundances:
            curr_median = median(curr_kmer_abundances)
            print("Median is: ", curr_median)
        else:
            curr_median = 0
        print("Total No. of k-mers from this contig: ", count_of_male_kmers)
        print("No. of k-mers not shared with female: ", count_of_male_kmers_not_shared_with_female)

        if count_of_male_kmers != 0:
            proportion = count_of_male_kmers_not_shared_with_female / count_of_male_kmers
        else:
            proportion = 0
        print("Proportion is: ", proportion)

        SeqIO.write(SeqRecord(record.seq, id=record.id, description=str(length) + " " + str(proportion) + " " + str(curr_median)), annotated_contigs, "fasta")
    
    return


def classify_fm_mode(kmer_size, female_kmers_bf):
    male_contigs_file = "data/male_contigs.fasta"
    annotated_contigs = open("proportion_annotated_contigs.fasta", "w")
    ctg_count = 0
    n_kmers = "N"*kmer_size
    
    for record in SeqIO.parse(male_contigs_file, "fasta"):
        ctg_count += 1
        print("No. of contigs seen so far: ", ctg_count)
        print("Current contig ID is : ", record.id)

        to_kmerize_fwd = str(record.seq).upper()
        length = len(to_kmerize_fwd)
        reverse = kmers.reverse_complement(to_kmerize_fwd)
        count_of_male_kmers = 0
        count_of_male_kmers_not_shared_with_female = 0
        curr_kmer_abundances = []
        for i in range(0, length-kmer_size+1):
            kmer = to_kmerize_fwd[i:i+kmer_size]
            if kmer == n_kmers:
                continue
            count_of_male_kmers += 1
            if kmer not in female_kmers_bf:
                count_of_male_kmers_not_shared_with_female += 1
            else:
                rev_kmer = reverse[length-kmer_size-i:length-i]
                if rev_kmer not in female_kmers_bf:
                    count_of_male_kmers_not_shared_with_female += 1

        print("Total No. of k-mers from this contig: ", count_of_male_kmers)
        print("No. of k-mers not shared with female: ", count_of_male_kmers_not_shared_with_female)

        try:
            proportion = count_of_male_kmers_not_shared_with_female / count_of_male_kmers
        except ZeroDivisionError:
            proportion = 0
        print("Proportion is: ", proportion)

        SeqIO.write(SeqRecord(record.seq, id=record.id, description=str(length) + " " + str(proportion)), annotated_contigs, "fasta")
    annotated_contigs.close()
    return

def classify_ctgs(kmer_size, bf, fem_kmers, mode):
    female_kmers_bf = getbloomFilter(bf, fem_kmers, kmer_size)
    if mode == "female+male":
        classify_fm_male_mode(kmer_size, female_kmers_bf)
    elif mode == "female_only":
        classify_fm_mode(kmer_size, female_kmers_bf)
    else:
        print("Incorrect mode entered")
        kmers.exit_gracefully()


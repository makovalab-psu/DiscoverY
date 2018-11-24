from __future__ import division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scripts import kmers
from pybloomfilter import BloomFilter
from collections import defaultdict
from numpy import median


def classify_ctgs(kmer_size, bf, fem_kmers, mode):
    male_contigs_file = "data/male_contigs.fasta"
    reads_kmers = "data/kmers_from_male_reads"
    # output

    # make Bloom Filter of kmers from female
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
                for line in fm_kmers:
                    female_kmers_bf.add(line[:kmer_size])
        else :
            female_reference_file = "data/female.fasta"    
            for record in SeqIO.parse(female_reference_file,"fasta"):
                to_kmerize_fwd = str(record.seq).upper()
                length = len(to_kmerize_fwd)
                for i in range(0, length-kmer_size+1):
                    female_kmers_bf.add(to_kmerize_fwd[i:i+kmer_size])

        print("Done creating bloom filter")


    annotated_records = []
    if mode == "female+male":
            # make a dict of all kmers from male
        print("Generating a dictionary from kmers in kmers_from_male_reads")
        kmer_abundance_dict_from_male = defaultdict(int)
        kmer_abundance_dict_from_male = kmers.make_dict_from_kmer_abundance(reads_kmers, kmer_size)
    ctg_count = 0
           # for each contig in contigs.fasta, calculate proportion of kmers that are present in female data set
    n_kmers = "N"*kmer_size
    for record in SeqIO.parse(male_contigs_file, "fasta"):
        ctg_count += 1
        print("No. of contigs seen so far: ", ctg_count)
        print("Current contig ID is : ", record.id)    
        # convert the SeqRecord object's .seq into a string that needs to be kmerized
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
            rev_kmer = reverse[i:i+kmer_size]
            if rev_kmer == n_kmers:
                continue
            count_of_male_kmers += 1
            # feature 1 : compute proportion of kmers shared with female
            if kmer not in female_kmers_bf and rev_kmer not in female_kmers_bf:
                count_of_male_kmers_not_shared_with_female += 1
            if mode=="female+male":
                # feature 2 : find median of kmer abundances from reads
                if kmer in kmer_abundance_dict_from_male :
                    curr_kmer_abundances.append(kmer_abundance_dict_from_male[kmer])
                elif rev_kmer in kmer_abundance_dict_from_male :
                    curr_kmer_abundances.append(kmer_abundance_dict_from_male[rev_kmer])
                else:
                    continue
        if mode=="female+male":        
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
            proportion = -1
            print("Proportion is: ", proportion)

        if mode=="female+male":
            annotated_records.append(SeqRecord(record.seq, id=record.id, description=str(length) + " " + str(proportion) + " " + str(curr_median)))
        else:
            annotated_records.append(SeqRecord(record.seq, id=record.id, description=str(length) + " " + str(proportion)))

    with open("proportion_annotated_contigs.fasta", "w") as annotated_contigs_file:
        SeqIO.write(annotated_records, annotated_contigs_file, "fasta") 
    return 1
    

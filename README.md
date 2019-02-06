# DiscoverY
DiscoverY is a toool to shortlist Y-specific contigs from an assembly of male whole genome sequencing data, based on exact _k-mer_ matches with a female. This tool is platform agnostic and has been tested on male assemblies from Illumina, PacBio and 10x Genomics. The female can be a reference assembly or a low coverage raw reads dataset. DiscoverY can be ran in two different modes: female_only or female+male.  In the female_only mode, the proportion shared between each contig with a female reference is computed; female+male mode uses both proportion of k-mers of each contig shared with a reference female and the kmer counts from the male raw reads to estimate each contig's depth-of-coverage.


## Usage 

Before running DiscoverY, the following input files are required in ./data folder. 
Note that currently the names of the "data" folder and of the files are hardcoded into DiscoverY.

	male_contigs.fasta : contigs from WGS assembly, which will be annotated by running DiscoverY
	kmers_from_male_reads : kmers from raw male reads used for assembly, used for computing coverage
	[optional] female.bloom : Bloom Filter of kmers from female if available
	[optional] female_kmers : kmers from female if Bloom Filter has not been constructed
	[optional] female.fasta : Female reference in FASTA format if BF and kmer set both unavailable
	

[See below for more info on generating these files.](https://github.com/md5sam/DiscoverY/blob/master/README.md#input-files) 

A typical run of DiscoverY looks like this. 

	python discoverY.py --female_bloom --mode female+male
	

DiscoverY accepts the following parameters. 

	python discoverY.py [--help] 

- --help: print usage information.

- --female_bloom

- --female_kmers_set

- --kmers_size

- --mode

The output of DiscoverY is an annotated file with : ```proportion_annotated_contigs.fastq``` in which the fasta headers have information about the proportion shared with female. 

## Installation 

To download, 

	git clone https://github.com/makovalab-psu/DiscoverY
	
DiscoverY also requires the numpy and biopython python packages in order to run the classifier plots.
These packages can be installed on many  systems as follows:

    pip install numpy
    pip install biopython
    pip install matplotlib
    pip install seaborn

DiscoverY also uses the k-mer counter DSK. The latest DSK binaries (v2.2.0 for Linux 64 bit and v2.2.0 for Mac OSX) are provided in the dependency folder. Thus, if you are using either of these operating systems, DSK need not be installed, and you may use the binaries as provided. For other operating systems, or if alternate versions or functionality of DSK is desired, see https://gatb.inria.fr/software/dsk/.

## Example

An example is provided in the /data folder.



## Input files


### Generating ```kmers_from_male_reads```
The ```kmers_from_male_reads``` is a file which includes a line for every distinct kmer in the male reads dataset, and each line contains the kmer sequence followed by whitespace followed by the number of times it occurs. For example, the first three lines could look like this.

	AAAAAAAAAAAAAAAAGAAAAACAA 5
	AAAAAAAAAAAAAAACAAGCTGAAT 8
	AAAAAAAAAAAAAAAGAAAAACAAA 3

The ```kmers_from_male_reads``` file can be generated by the k-mer counter DSK. The ./dependency folder contains DSK binaries for Linux 64 bit and Mac OSX. Usage is as follows (example shown below is for a Linux system) :

    cd dependency
    ./run_dsk_Linux.sh <FASTQ_file> <kmer_size>


If the k-mer counts file for male is not already provided, the user may need to generate k-mer counts manually using DSK. To generate k-mer counts with DSK, the following steps are needed : 

    cd dependency 
    ln -s ../data/female.fasta  # make sure the correct reads file is provided to DSK
    ./run_dsk_Linux.sh r1.fastq 25  


The kmer\_counts table will be generated in :

    dependency/dsk_output/kmers_from_female


This file can be copied or linked to the data folder so that DiscoverY can use it : 

    cd ../data
    ln -s ../dependency/kmers_from_female 



## Miscellaneous

The following scripts are included with this distribution of DiscoverY, and are automatically run by discovery.py as part of the pipeline. Users may consider them separately for custom needs if required. 

	
**kmers.py** 
	
	a set of general purpose functions to work with kmers

**classify_ctgs.py**
	
	input : all contigs from WGS assembly, male read kmers and female kmers.
	output : annotated contigs with proportions shared with female

	

## License
This program is released under the MIT License. Please see LICENSE.md for details


## Citation
If you use DiscoverY in your research, please cite this repository. 



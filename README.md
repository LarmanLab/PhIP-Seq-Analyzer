

#########################
PhIP-seq Analyzer

Demo version
Date: 20170524
#########################


apply for a server node

#interact -p parallel -n 24 -t 24:0:0


################
Data preparation
################

1. It is supposed that the home folder is ./, and the sequencing files determined by sequencing analyzer were put in the fold ./raw/. so they are
fastq file: ./raw/PhIPseq_R1.fastq.gz
index file: ./raw/PhIPseq_I1.fastq.gz
barcode file: ./raw/Sample-Barcode.txt

2. The splitted fastq files would be stored into the directory ~/phip/PhIP-seq_Analyzer/test_rawdata/

3. The running mode would be phipseq analysis mode, so the output folders would be ./test_human and ./test_virus. So the last parameter should be -y ./test

4. enter the folder known as PhIP-seaq_Analyzer


#########################
1: Demultiplex FASTQ file
#########################

Example 1: Mostly used function involves demultiplexing and generating sample_info.csv and variables.txt run the command line like the below:

$ python3 ./bin/bioTreatFASTQ.py -i ./raw/PhIPseq_I1.fastq.gz -f ./raw/PhIPseq_R1.fastq.gz -b ./raw/Sample-Barcode.txt -o ./test_rawdata/ -y ./test

Once the proceduce is done, some folders and files would be created: There would be many fastq files in ./test_rawdata/, and the file names would be determined by the barcode file.
The folders known as ./test_human and ./test_virus would be created, of which has sample_info.csv and variables.txt
Note: the pipeline only accept the absolute path of a directory or file


Example 2:  trim nucleotides when demultiplexing, let say sequencing cycles is 50nt
trim 10nt from 3-end or keep the first 40nt:                -r 40
remove 10nt from 5-end and keep the left:                   -t 10
trim 5nt from 5-end and keep 40nt and discard the left: 	-t 5 -r 35


Example 3: The length of exported reads in FASTQ should be kept equal determined by the Sequencing Analyzer. In some cases, They are not because no quality filtering apply. We could specify -l 100nt. That option would discard all reads shorter than 100nt when demultiplexing.


Example 4: Demultiplexing step can be skipped, and directly get sample_info.csv and variables.txt. Here *.fastq files were stored at ./test_rawdata/

$ python3 ./bin/bioTreatFASTQ.py -o ./test_rawdata/ -y ./test


Example 5: skip demultiplexing step, but trim fastq reads only. Here, all fastq files were put in the ./test_rawdata/, and the trimmed fastq files (remove 40nt from 3-end of each read with 100nt) were saved into ./trim_rawdata/

$ python3 ./bin/bioTreatFASTQ.py -r 60 -x ./test_rawdata/ -o./trim_rawdata/ -y ./test




Example 5: The default reference peptide libraries are human and virus. There are additional peptide library known as allergome(allergic peptides) and PE(public epitopes) 
 -c human,virus 
 -c virus,allergome,PE

###################
2: phipseq analysis
###################
run the command line like the below:

#human library
$ python3 ./bin/bioPHIPseq.py ./test_human/variables.txt
#virus library
$ python3 ./bin/bioPHIPseq.py ./test_virus/variables.txt


###################
Requirements
###################

The barcode file:
1. barcode file should be *.txt seperated by tab. The first and second columns should be barcode sequences and sample names, respectively.
2. Regarding sample names, avoid some characters namely slash(/ or \), asterisk(*), at sign(@), any brackets or white space. And the characters dash(-), underscore(_), or dot(.) are acceptable.
3. No while line is allowed.

The reads file and index file
1. FASTQ format
2. Both of the files should be matched and applied together. 
3. support compressed format with *.gz

Running environments
1. Linux
2. Python 3.4.0 above


######################
ERROR Handling
######################
ERROR 1: Mac OS X: ValueError: unknown locale: UTF-8 in Python

Resolution: If you have faced the error on MacOS X, here's the quick fix - add these lines to your ~/.bash_profile:
	export LC_ALL=en_US.UTF-8
	export LANG=en_US.UTF-8
end then reload bash_profile: 
	# source ~/.bash_profile

ERROR 2: Traceback (most recent call last):
  File "./bin/bioTreatFASTQ.py", line 141, in <module>
    myGenome.genome(par['fq_file']).demultiplex_fq(par)
  File "/home/yuan/phip/PhIP-Seq_Analyzer/bin/myGenome.py", line 222, in demultiplex_fq
    for L1,La, L2,Lb, L3,Lc, L4,Ld in itertools.zip_longest(*[F1,F2]*4):
AttributeError: 'module' object has no attribute 'zip_longest'

Resolution: python3 instead of python2. 



####
#end







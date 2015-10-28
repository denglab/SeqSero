# SeqSero
Salmonella serotyping from genome sequencing data


# Introduction 
SeqSero is a pipeline for Salmonella serotype determination from raw sequencing reads or genome assemblies. A web app is available at www.denglab.info/SeqSero 

# Dependencies 
SeqSero depends on: 
1. Python 2.7 and [Biopython 1.65](http://biopython.org/wiki/Download) 
2. [Burrows-Wheeler Aligner](http://sourceforge.net/projects/bio-bwa/files/) 
3. [Samtools](http://sourceforge.net/projects/samtools/files/samtools/) 
4. [NCBI BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/) 
5. [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software ) 

Make sure all the executables are added to your search path. SeqSero also uses the in silico PCR program (isPcr) written by [Jim Kent](http://hgwdev.cse.ucsc.edu/~kent/exe/linux/), which is supplied by the SeqSero package. 

# Executing the code 
Usage: SeqSero.py 
-m <int> (input data type, '1' for interleaved paired-end reads , '2' for separated paired-end reads, '3' for single reads, '4' for genome assembly) 
-i <file> (/path/to/input/file) 
-b <string> (algorithms for bwa mapping; 'mem' for mem, 'sam' for samse/sampe; default=sam; optional) 

# Output 
Upon executing the command, a directory named 'SeqSero_result_<time_you_run_SeqSero>' will be created. Your result will be stored in 'Seqsero_result.txt' in that directory

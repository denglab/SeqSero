# SeqSero
Salmonella serotyping from genome sequencing data


# Introduction 
SeqSero is a pipeline for Salmonella serotype determination from raw sequencing reads or genome assemblies. A web app is available at www.denglab.info/SeqSero 

# Dependencies 
SeqSero depends on:

1. Python 2.7 and [Biopython 1.65](http://biopython.org/wiki/Download); 

2. [Burrows-Wheeler Aligner](http://sourceforge.net/projects/bio-bwa/files/); 

3. [Samtools](http://sourceforge.net/projects/samtools/files/samtools/);

4. [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download);

5. [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software);

6. [isPcr](http://hgwdev.cse.ucsc.edu/~kent/exe/linux/) written by Jim Kent. 

# Executing the code 
Usage: SeqSero.py 

-m <int> (input data type, '1' for interleaved paired-end reads , '2' for separated paired-end reads, '3' for single reads, '4' for genome assembly) 

-i <file> (/path/to/input/file) 

-b <string> (algorithms for bwa mapping; 'mem' for mem, 'sam' for samse/sampe; default=sam; optional) 

# Output 
Upon executing the command, a directory named 'SeqSero_result_<time_you_run_SeqSero>' will be created. Your result will be stored in 'Seqsero_result.txt' in that directory

# Citation
_S. Zhang, Y. Yin, M.B. Jones, Z. Zhang, B.L. Deatherage Kaiser, B.A. Dinsmore, C. Fitzgerald, P.I. Fields, X. Deng. Salmonella serotype determination utilizing high-throughput genome sequencing data.J Clin Microbiol, 53 (2015), pp. 16-1692_ ([paper source](http://jcm.asm.org/content/early/2015/03/05/JCM.00323-15))


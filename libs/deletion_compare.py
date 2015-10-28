
import os
from Bio import SeqIO
import sys
from Initial_functions import Uniq
from Bio.Blast import NCBIXML


target=sys.argv[1] #should be sra format
test_gene=sys.argv[2]
mapping_mode=sys.argv[3]
if sys.argv[4] not in ("1","2","3"):
  additional_file=sys.argv[4]
  file_mode=sys.argv[5]
else:
  additional_file=""
  file_mode=sys.argv[4]




def Copenhagen(sra_name,additional_file,mapping_mode,file_mode):
  if file_mode=="1":#interleaved
    if sra_name[-3:]=="sra":
      del_fastq=1
      for_fq=sra_name.replace(".sra","_1.fastq")
      rev_fq=sra_name.replace(".sra","_2.fastq")
      for_sai=sra_name.replace(".sra","_1.sai")
      rev_sai=sra_name.replace(".sra","_2.sai")
      sam=sra_name.replace(".sra",".sam")
      bam=sra_name.replace(".sra",".bam")
    else:
      del_fastq=0
      core_id=sra_name.split(".fastq")[0]
      for_fq=core_id+"-read1.fastq"
      rev_fq=core_id+"-read2.fastq"
      for_sai=core_id+"_1.sai"
      rev_sai=core_id+"_2.sai"
      sam=core_id+".sam"
      bam=core_id+".bam"
  elif file_mode=="2":#seperated
    forword_seq=sra_name
    reverse_seq=additional_file
    for_core_id=forword_seq.split(".fastq")[0]
    re_core_id=reverse_seq.split(".fastq")[0]
    for_fq=for_core_id+".fastq"
    rev_fq=re_core_id+".fastq"
    for_sai=for_core_id+".sai"
    rev_sai=re_core_id+".sai"
    sam=for_core_id+".sam"
    bam=sam.replace(".sam",".bam")
  elif file_mode=="3":#single-end
    if sra_name[-3:]=="sra":
        del_fastq=1
        for_fq=sra_name.replace(".sra","_1.fastq")
        rev_fq=sra_name.replace(".sra","_2.fastq")
        for_sai=sra_name.replace(".sra","_1.sai")
        rev_sai=sra_name.replace(".sra","_2.sai")
        sam=sra_name.replace(".sra",".sam")
        bam=sra_name.replace(".sra",".bam")
    else:
        del_fastq=0
        core_id=sra_name.split(".fastq")[0]
        for_fq=core_id+".fastq"
        rev_fq=core_id+".fastq"
        for_sai=core_id+"_1.sai"
        rev_sai=core_id+"_2.sai"
        sam=core_id+".sam"
        bam=core_id+".bam"
  
  database="complete_oafA.fasta"
  os.system("bwa index database/"+database)###01/28/2015
  if mapping_mode=="mem":
    os.system("bwa mem database/"+database+" "+for_fq+" "+rev_fq+" > "+sam) #2014/12/23
  elif mapping_mode=="sam":
    os.system("bwa aln database/"+database+" "+for_fq+" > "+for_sai)
    os.system("bwa aln database/"+database+" "+rev_fq+" > "+rev_sai)
    os.system("bwa sampe database/"+database+" "+for_sai+" "+ rev_sai+" "+for_fq+" "+rev_fq+" > "+sam)
  os.system("samtools view -F 4 -Sbh "+sam+" > "+bam)
  os.system("samtools view -h -o "+sam+" "+bam)
  os.system("cat "+sam+"|awk '{if ($5>0) {print $10}}'>"+sam+"_seq.txt")
  os.system("cat "+sam+"|awk '{if ($5>0) {print $1}}'>"+sam+"_title.txt")
  file1=open(sam+"_title.txt","r")
  file2=open(sam+"_seq.txt","r")
  file1=file1.readlines()
  file2=file2.readlines()
  file=open(sam+".fasta","w")
  for i in range(len(file1)):
    title=">"+file1[i]
    seq=file2[i]
    if len(seq)>40 and (len(title)>5 or ("@" not in title)):
      file.write(title)
      file.write(seq)
  file.close()
  database2="oafA_of_O4_O5.fasta"
  os.system('makeblastdb -in database/'+database2+' -out '+database2+'_db '+'-dbtype nucl')
  os.system("blastn -query "+sam+".fasta"+" -db "+database2+"_db -out "+sam+"_vs_O45.xml -outfmt 5")
  handle=open(sam+"_vs_O45.xml")
  handle=NCBIXML.parse(handle)
  handle=list(handle)
  O9_bigger=0
  O2_bigger=0
  for x in handle:
    O9_score=0
    O2_score=0
    try:
      if 'O-4_full' in x.alignments[0].hit_def:
        O9_score=x.alignments[0].hsps[0].bits
        O2_score=x.alignments[1].hsps[0].bits
      elif 'O-4_5-' in x.alignments[0].hit_def:
        O9_score=x.alignments[1].hsps[0].bits
        O2_score=x.alignments[0].hsps[0].bits
      if O9_score>O2_score:
        O9_bigger+=1
      if O9_score<O2_score:
        O2_bigger+=1
    except:
      continue
  print "$$$Genome:",sra_name
  if O9_bigger>O2_bigger:
    print "$$$Typhimurium"
  elif O9_bigger<O2_bigger:
    print "$$$Typhimurium_O5-"
  else:
    print "$$$Typhimurium, even no 7 bases difference"
  print "O-4 number is:",O9_bigger
  print "O-4_5- number is:",O2_bigger
  os.system("rm "+sam+"_title.txt")###01/28/2015
  os.system("rm "+sam+"_seq.txt")###01/28/2015
  os.system("rm "+sam+".fasta")###01/28/2015
  os.system("rm "+database2+"_db.*")###01/28/2015
  os.system("rm "+sam+"_vs_O45.xml")###01/28/2015

if test_gene=="oafA":
  Copenhagen(target,additional_file,mapping_mode,file_mode)
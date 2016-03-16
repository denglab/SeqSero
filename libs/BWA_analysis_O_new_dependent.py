#!/usr/bin/env python
#tyr_of_O2_O9.fasta should be in the same directory, in it, O9 should be first then O2

import os
from Bio import SeqIO
import sys
from Initial_functions import Uniq
from Bio.Blast import NCBIXML

def BWA_O_analysis(sra_name,additional_file,database,mapping_mode,file_mode):
  if file_mode=="1":#interleaved
    if sra_name[-3:]=="sra":
      os.system("fastq-dump --split-files "+sra_name)
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
      try:
        os.system("gunzip "+sra_name)
      except:
        pass
      dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
      os.system("perl "+dirpath+"/split_interleaved_fastq.pl --input "+core_id+".fastq --output "+core_id.replace(".","_")+".fastq")#######03152016
      ori_size=os.path.getsize(core_id+".fastq")#######03152016
      os.system("mv "+core_id.replace(".","_")+"-read1.fastq"+" "+core_id+"-read1.fastq")#######03152016
      os.system("mv "+core_id.replace(".","_")+"-read2.fastq"+" "+core_id+"-read2.fastq")#######03152016
      for_fq=core_id+"-read1.fastq"#######03152016
      rev_fq=core_id+"-read2.fastq"#######03152016
      if float(os.path.getsize(for_fq))/ori_size<=0.1 or float(os.path.getsize(rev_fq))/ori_size<=0.1:#09092015#12292015#######03152016
        os.system("echo haha")#09092015
        os.system("perl "+dirpath+"/splitPairedEndReads.pl "+core_id+".fastq")#09092015
        os.system("mv "+core_id+".fastq_1 "+for_fq)##09092015
        os.system("mv "+core_id+".fastq_2 "+rev_fq)##09092015
      else:#09092015
        os.system("echo hehe")#09092015
      for_sai=core_id+"_1.sai"
      rev_sai=core_id+"_2.sai"
      sam=core_id+".sam"
      bam=core_id+".bam"
  elif file_mode=="2":#seperated
    forword_seq=sra_name
    reverse_seq=additional_file
    try:
        os.system("gunzip "+forword_seq)
    except:
        pass
    try:
        os.system("gunzip "+reverse_seq)
    except:
        pass
    for_core_id=forword_seq.split(".fastq")[0]
    re_core_id=reverse_seq.split(".fastq")[0]
    for_fq=for_core_id+".fastq"
    rev_fq=re_core_id+".fastq"
    dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))#######03152016
    print "check fastq id and make them in accordance with each other...please wait..."
    os.system("python "+dirpath+"/compare_and_change_two_fastq_id.py "+for_fq+" "+rev_fq)#######03152016
    for_sai=for_core_id+".sai"
    rev_sai=re_core_id+".sai"
    sam=for_core_id+".sam"
    bam=sam.replace(".sam",".bam")
  elif file_mode=="3":#single-end
    if sra_name[-3:]=="sra":
      os.system("fastq-dump --split-files "+sra_name)###01/28/2015
      del_fastq=1
      for_fq=sra_name.replace(".sra","_1.fastq")
      for_sai=sra_name.replace(".sra","_1.sai")
      sam=sra_name.replace(".sra",".sam")
      bam=sra_name.replace(".sra",".bam")
    else:
      del_fastq=0
      core_id=sra_name.split(".fastq")[0]
      try:
        os.system("gunzip "+sra_name)
      except:
        pass  
      for_fq=core_id+".fastq"
      for_sai=core_id+"_1.sai"
      sam=core_id+".sam"
      bam=core_id+".bam"
  
  os.system("bwa index "+database)
  if file_mode!="3":
    if mapping_mode=="sam":
      os.system("bwa aln "+database+" "+for_fq+" > "+for_sai)
      os.system("bwa aln "+database+" "+rev_fq+" > "+rev_sai)
      os.system("bwa sampe "+database+" "+for_sai+" "+ rev_sai+" "+for_fq+" "+rev_fq+" > "+sam)
    elif mapping_mode=="mem":
      os.system("bwa mem "+database+" "+for_fq+" "+rev_fq+" > "+sam) #2014/12/23
  else:
    if mapping_mode=="mem":
      os.system("bwa mem "+database+" "+for_fq+" > "+sam) #2014/12/23
    elif mapping_mode=="sam":
      os.system("bwa aln "+database+" "+for_fq+" > "+for_sai)
      os.system("bwa samse "+database+" "+for_sai+" "+for_fq+" > "+sam)
  os.system("samtools view -F 4 -Sbh "+sam+" > "+bam)
  os.system("samtools view -h -o "+sam+" "+bam)

  file=open(sam,"r")
  handle=file.readlines()
  name_list=[]
  for line in handle:
    if len(line)>300:
      name_list.append(line.split("\t")[2])
  a,b=Uniq(name_list)
  c=dict(zip(a,b))
  final_O=sorted(c.iteritems(), key=lambda d:d[1], reverse = True) #order from frequency high to low, but tuple while not list
  Sero_list_O=[]
  print "Final_Otype_list:"
  print final_O
  num_1=0#new inserted
  O9_wbav=0
  O310_wzx=0
  O946_wzy=0
  if len(final_O)>0: #new inserted
    for x in final_O:#new inserted
      num_1=num_1+x[1]#new inserted
      if "O-9,46_wbaV" in x[0]:
        O9_wbaV=x[1]
      if "O-3,10_wzx" in x[0]:
        O310_wzx=x[1]
      if "O-9,46_wzy" in x[0]:
        O946_wzy=x[1]
      if "O-3,10_not_in_1,3,19" in x[0]:
        O310_no_1319=x[1]
      if "O-9,46,27_partial_wzy" in x[0]:
        O94627=x[1]        
  O_list=[]
  O_choice=""


  print "$$$Genome:",sra_name
  if len(final_O)==0:
    print "$$$No Otype, due to no hit"
  else:
    if final_O[0][1]<8:
      print "$$$No Otype, due to the hit reads number is small."
    else:
      for x in final_O:
        if x[1]>5:
          O_list.append(x[0])
      qq=1#
      for x in final_O:#
        if "sdf" in x[0] and x[1]>3:#
          qq=0#
          print "$$$",x[0],"got a hit, reads:",x[1]#
      if qq!=0:#
        print "$$$No sdf exists"#

      if "O-9,46_wbaV" in O_list and float(O9_wbaV)/float(num_1) > 0.1:
        if "O-9,46_wzy" in O_list and float(O946_wzy)/float(num_1) > 0.1:
          O_choice="O-9,46"
          print "$$$Most possilble Otype:  O-9,46"
        elif "O-9,46,27_partial_wzy" in O_list and float(O94627)/float(num_1) > 0.1:
          O_choice="O-9,46,27"
          print "$$$Most possilble Otype:  O-9,46,27"
        else:
          O_choice="O-9"
          if file_mode=="3":
            rev_fq=""
            rev_sai=""
            assembly(sra_name,O_choice,for_fq,rev_fq,for_sai,rev_sai,sam,bam,mapping_mode)
          else:
            assembly(sra_name,O_choice,for_fq,rev_fq,for_sai,rev_sai,sam,bam,mapping_mode)
      elif ("O-3,10_wzx" in O_list) and ("O-9,46_wzy" in O_list) and float(O310_wzx)/float(num_1) > 0.1 and float(O946_wzy)/float(num_1) > 0.1:
        if "O-3,10_not_in_1,3,19" in O_list and float(O310_no_1319)/float(num_1) > 0.1:
          O_choice="O-3,10"
          print "$$$Most possilble Otype:  O-3,10"
        else:
          O_choice="O-1,3,19"
          print "$$$Most possilble Otype:  O-1,3,19"
      else:
        try: 
          O_choice=final_O[0][0].split("_")[0]
          if O_choice=="O-1,3,19":
            O_choice=final_O[1][0].split("_")[0]
          print "$$$Most possilble Otype: ",O_choice
        except:
          print "$$$No suitable Otype, or failure of mapping (please check the quality of raw reads)"


def assembly(sra_name,potential_choice,for_fq,rev_fq,for_sai,rev_sai,sam,bam,mapping_mode):
  database="ParaA_rfb.fasta"
  os.system("bwa index database/"+database)###01/28/2015
  if rev_fq=="":#2015/09/09
    if mapping_mode=="mem":
      os.system("bwa mem database/"+database+" "+for_fq+" > "+sam) #2014/12/23
    elif mapping_mode=="sam":
      os.system("bwa aln database/"+database+" "+for_fq+" > "+for_sai)
      os.system("bwa samse database/"+database+" "+for_sai+" "+for_fq+" > "+sam)
  else:
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
    if len(seq)>=50 and len(title)>6:#generally,can be adjusted
      file.write(title)
      file.write(seq)
  file.close()
  database2="tyr_of_O2_O9.fasta"
  os.system('makeblastdb -in database/'+database2+' -out '+database2+'_db '+'-dbtype nucl')
  os.system("blastn -query "+sam+".fasta"+" -db "+database2+"_db -out "+sam+"_vs_O29.xml -outfmt 5")
  handle=open(sam+"_vs_O29.xml")
  handle=NCBIXML.parse(handle)
  handle=list(handle)
  O9_bigger=0
  O2_bigger=0
  for x in handle:
    O9_score=0
    O2_score=0
    try:
      if 'O-9' in x.alignments[0].hit_def:
        O9_score=x.alignments[0].hsps[0].bits
        O2_score=x.alignments[1].hsps[0].bits
      elif 'O-2' in x.alignments[0].hit_def:
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
    print "$$$Most possible Otype is O-9"
  elif O9_bigger<O2_bigger:
    print "$$$Most possible Otype is O-2"
  else:
    print "$$$No suitable one, because can't distinct it's O-9 or O-2, but ",potential_choice," has a more possibility."
  print "O-9 number is:",O9_bigger
  print "O-2 number is:",O2_bigger

  os.system("rm "+sam+"_title.txt")###01/28/2015
  os.system("rm "+sam+"_seq.txt")###01/28/2015
  os.system("rm "+sam+".fasta")###01/28/2015
  os.system("rm "+database2+"_db.*")###01/28/2015
  os.system("rm "+sam+"_vs_O29.xml")###01/28/2015

 
target=sys.argv[1] #should be sra format
data_base=sys.argv[2]
mapping_mode=sys.argv[3]
if sys.argv[4] not in ("1","2","3"):
  additional_file=sys.argv[4]
  file_mode=sys.argv[5]
else:
  additional_file=""
  file_mode=sys.argv[4]

BWA_O_analysis(target,additional_file,data_base,mapping_mode,file_mode)

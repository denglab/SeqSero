#!/usr/bin/env python

import os
from Bio import SeqIO
import sys
from Initial_functions import Uniq
from Bio.Blast import NCBIXML

def BWA_analysis(sra_name,additional_file,database,mapping_mode,file_mode,z):
  global fliC_option
  global fljB_option
  global family_c_list
  global family_b_list
  global list_c_length
  global list_b_length
  family_c_list=[]   #will catch families on the family test mapping
  family_b_list=[]
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

  database=database.split("/")[-1]##########1/27/2015
  os.system("bwa index database/"+database)
  if file_mode!="3":
    if mapping_mode=="mem":
      os.system("bwa mem database/"+database+" "+for_fq+" "+rev_fq+" > "+sam) #2014/12/23
    elif mapping_mode=="sam":
      os.system("bwa aln database/"+database+" "+for_fq+" > "+for_sai)
      os.system("bwa aln database/"+database+" "+rev_fq+" > "+rev_sai)
      os.system("bwa sampe database/"+database+" "+for_sai+" "+ rev_sai+" "+for_fq+" "+rev_fq+" > "+sam)
    elif mapping_mode=="nanopore":
      os.system("bwa mem -x ont2d database/"+database+" "+for_fq+" "+rev_fq+" > "+sam)
  else:
    if mapping_mode=="mem":
      os.system("bwa mem database/"+database+" "+for_fq+" > "+sam) #2014/12/23
    elif mapping_mode=="sam":
      os.system("bwa aln database/"+database+" "+for_fq+" > "+for_sai)
      os.system("bwa samse database/"+database+" "+for_sai+" "+for_fq+" > "+sam)
    elif mapping_mode=="nanopore":
      os.system("bwa mem -x ont2d database/"+database+" "+for_fq+" > "+sam) 
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
  Sero_list_C=[]
  Sero_list_B=[]
  description=[] #for storage of whole desription to extract the multifasta file for second BWA mapping
  fliC={}
  fljB={}
  for x in c:
    if x[:4]=="fliC":
      fliC[x]=c[x]
  for x in c:
    if x[:4]=="fljB":
      fljB[x]=c[x]
  final_fliC=sorted(fliC.iteritems(), key=lambda d:d[1], reverse = True) #order from frequency high to low, but tuple while not list
  final_fljB=sorted(fljB.iteritems(), key=lambda d:d[1], reverse = True) #order from frequency high to low, but tuple while not list
  print "Final_filC_list:"
  print final_fliC
  num_1=0#new inserted
  num_2=0#new inserted
  if len(final_fliC)>0: #new inserted
    for x in final_fliC:#new inserted
      num_1=num_1+x[1]#new inserted
  print "Final_fliC_number_together: ",num_1#new inserted
  print "Final_fljB_list:"
  print final_fljB
  if len(final_fljB)>0: #new inserted
    for x in final_fljB: #new inserted
      num_2=num_2+x[1] #new inserted
  print "Final_fljB_number_together: ",num_2#new inserted
  print "$$Genome:",sra_name
  try:
    fliC_option=final_fliC[0][0].split("_")[1]
  except:
    fliC_option="-"
  try:
    fljB_option=final_fljB[0][0].split("_")[1]
  except:
    fljB_option="-"
   
  if z==0:
    if len(final_fliC)==0 or num_1<=10:
      print "$$$No fliC, due to no hit"
    else:
      if final_fliC[0][1]<=1 and z==1:
        print "$$$No fliC, due to the hit reads number is small."
      else:
        try:
          family=final_fliC[0][0].split("_")[-1]
          Sero_list_C.append(family)
          description.append(final_fliC[0][0])
          print "$$Most possilble fliC family: ",Sero_list_C[0]," Number: ",final_fliC[0][1]
          i=0
          for x in final_fliC:
            if x[0].split("_")[-1] not in Sero_list_C:
              if i<=x[1]:
                i=x[1]
                sec_choice=x[0].split("_")[-1]
                des=x[0]
                number=x[1]
          if locals().has_key('sec_choice'):
            Sero_list_C.append(sec_choice)
            description.append(des)
            print "$$Sec possilble fliC family: ",sec_choice," Number: ",number
            j=0
            for x in final_fliC:
              if x[0].split("_")[-1] not in Sero_list_C:
                if j<=x[1]:
                  j=x[1]
                  third_choice=x[0].split("_")[-1]
                  des=x[0]
                  number=x[1]
            if locals().has_key('third_choice'):
              Sero_list_C.append(third_choice)
              description.append(des)
              print "$$Third possilble fliC family: ",third_choice," Number: ",number 
        except:
          print "$$$No fliC, or failure of mapping"
    try:
      ratio=float(num_2)/float(num_1)
    except:
      ratio=0

    if len(final_fljB)==0 or num_2<=5 or ratio<0.15:
      print "$$$No fljB, due to no hit"
    else:
      if final_fljB[0][1]<=1 and z==1:
        print "$$$No fljB, due to the hit reads number is small."
      else:
        try:
          family=final_fljB[0][0].split("_")[-1]
          Sero_list_B.append(family)
          description.append(final_fljB[0][0])
          print "$$Most possilble fljB family: ",Sero_list_B[0]," Number: ",final_fljB[0][1]
          i=0
          for x in final_fljB:
            if x[0].split("_")[-1] not in Sero_list_B:
              if i<=x[1]:
                i=x[1]
                B_sec_choice=x[0].split("_")[-1]
                des=x[0]
                number=x[1]
          if locals().has_key('B_sec_choice'):
            Sero_list_B.append(B_sec_choice)
            description.append(des)
            print "$$Sec possilble fljB: ",B_sec_choice," Number: ",number   
            j=0
            for x in final_fljB:
              if x[0].split("_")[-1] not in Sero_list_B:
                if j<=x[1]:
                  j=x[1]
                  B_third_choice=x[0].split("_")[-1]
                  des=x[0]
                  number=x[1]
            if locals().has_key('B_third_choice'):
              Sero_list_B.append(B_third_choice)
              description.append(des)
              print "$$Third possilble fljB: ",B_third_choice," Number: ",number 
        except:
          print "$$$No fljB, or failure of mapping"
    if len(description)==0:    #used for the case which fljB and fliC both has no hit, it will directly cease the function
      return
    handle=SeqIO.parse("database/"+database,"fasta")########1/27/2015
    handle=list(handle)
    file1=open("temp"+sra_name+".fasta","w")
    for x in handle:
      for y in description:
        if x.description==y:
          title=">"+x.description
          seq=str(x.seq)
          file1.write(title+"\n")
          file1.write(seq+"\n")
    file1.close()
    data_base2="temp"+sra_name+".fasta"
    os.system("mv "+data_base2+" database")######1/27/2015
    z=1
    BWA_analysis(sra_name,additional_file,data_base2,mapping_mode,file_mode,z) #z=1 this time, z is the pointer for the seqeunce of run
    return

  if z==1:
    list_c_length=len(final_fliC)
    list_b_length=len(final_fljB)
    family_test(final_fliC,"fliC",type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam)
    family_test(final_fljB,"fljB",type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam)

  #os.system("rm "+for_fq)###01/28/2015#$$$$$$$
  #os.system("rm "+rev_fq)###01/28/2015
  #os.system("rm "+for_sai)
  #os.system("rm "+rev_sai)
  #os.system("rm "+bam)
  #os.system("rm "+sam)###01/28/2015
  #os.system("rm temp.txt")###01/28/2015
  #os.system("rm "+sam+".fasta"+"_db.*")###01/28/2015
  #os.system("rm temp"+sra_name+".fasta")###01/28/2015  

def family_test(fliC_fljB_list,listtype,type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam):
  Sero_list_C=[]
  Sero_list_B=[]
  if listtype=="fliC":
    if len(fliC_fljB_list)==0:
      print "$$No fliC, due to no hit" #because the only possible situation for len(final_fliC)==0 is above (z=0) len(final_fliC)==0, so there is no need to use "$$$" here
    else:
      if fliC_fljB_list[0][1]<=1:
        print "$$No fliC, due to the hit reads number is small." #similiar with above, no "$$$"
      else:
        if fliC_fljB_list[0][0].split("_")[-1]=="g,m":
          type="fliC"
          database="FliC_f_g_s_whole.fasta"
          database2="FliC_Family_g_special_genes.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_c_length,mapping_mode)
        elif fliC_fljB_list[0][0].split("_")[-1]=="l,v":
          type="fliC"
          database="fliC_l_z13_whole.fasta"
          database2="FliC_Family_l,v_special_genes.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_c_length,mapping_mode)
        elif fliC_fljB_list[0][0].split("_")[-1]=="r,i":
          type="fliC"
          database="fliC_r_whole.fasta"
          database2="FliC_Family_r,i_special_genes_short.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_c_length,mapping_mode)
        elif fliC_fljB_list[0][0].split("_")[-1]=="k,z":
          type="fliC"
          database="fliC_k,z_whole.fasta"
          database2="FliC_Family_k,z_special_genes.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_c_length,mapping_mode)
        elif fliC_fljB_list[0][0].split("_")[-1]=="b,d,j":
          type="fliC"
          database="fliC_b_whole.fasta"
          database2="FliC_Family_b,d,j_special_genes.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_c_length,mapping_mode)
        elif fliC_fljB_list[0][0].split("_")[-1]=="z4,z23":
          type="fliC"
          database="fliC_z4z23_whole.fasta"
          database2="fliC_z4,z23_family.fasta" #"FliC_Family_z4z23_special_genes.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_c_length,mapping_mode)
        elif fliC_fljB_list[0][0].split("_")[-1]=="z36,z38":
          type="fliC"
          database="fliC_z36,z38_whole.fasta"
          database2="FliC_Family_z36z38_special_genes.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_c_length,mapping_mode)
        else:
          try: 
            Sero_list_C.append(fliC_fljB_list[0][0].split("_")[1])
            print "$$$Most possilble fliC: ",Sero_list_C[0]," Number: ",fliC_fljB_list[0][1]
            i=0
            for x in fliC_fljB_list:
              if x[0].split("_")[1] not in Sero_list_C:
                if i<=x[1]:
                  i=x[1]
                  sec_choice=x[0].split("_")[1]
                  number=x[1]
            if locals().has_key('sec_choice'):
              Sero_list_C.append(sec_choice)
              print "$$$Sec possilble fliC: ",sec_choice," Number: ",number
              j=0
              for x in fliC_fljB_list:
                if x[0].split("_")[1] not in Sero_list_C:
                  if j<=x[1]:
                    j=x[1]
                    third_choice=x[0].split("_")[1]
                    number=x[1]
              if locals().has_key('third_choice'):
                Sero_list_C.append(third_choice)
                print "$$$Third possilble fliC: ",third_choice," Number: ",number 
          except:
            print "$$$No fliC, or failure of mapping (second run)"

  if listtype=="fljB":
    if len(fliC_fljB_list)==0:
      print "$$No fljB, due to no hit" #similiar with above, no "$$$"
    else:
      if fliC_fljB_list[0][1]<=1:
        print "$$No fljB, due to the hit reads number is small." #similiar with above, no "$$$"
      else:
        if fliC_fljB_list[0][0].split("_")[-1]=="1":
          type="fljB"
          database="FljB_1_2_7_whole.fasta"
          database2="FljB_Family_1_special_genes_all.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_b_length,mapping_mode)
        elif fliC_fljB_list[0][0].split("_")[-1]=="e":
          type="fljB"
          database="fljB_e,n,z15_whole.fasta"
          database2="FljB_Family_e_special_genes.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_b_length,mapping_mode)
        elif fliC_fljB_list[0][0].split("_")[-1]=="l,v":
          type="fljB"
          database="FljB_l,z13,z28_whole.fasta"
          database2="FljB_Family_l,v_special_genes.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_b_length,mapping_mode)
        elif fliC_fljB_list[0][0].split("_")[-1]=="k,z":
          type="fljB"
          database="FljB_z6_whole.fasta"
          database2="FljB_Family_k,z_special_genes.fasta"
          assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_b_length,mapping_mode)
        else:
          try: 
            Sero_list_B.append(fliC_fljB_list[0][0].split("_")[1])
            print "$$$Most possilble fljB: ",Sero_list_B[0]," Number: ",fliC_fljB_list[0][1]
            i=0
            for x in fliC_fljB_list:
              if x[0].split("_")[1] not in Sero_list_B:
                if i<=x[1]:
                  i=x[1]
                  B_sec_choice=x[0].split("_")[1]
                  number=x[1]
            if locals().has_key('B_sec_choice'):
              Sero_list_B.append(B_sec_choice)
              print "$$$Sec possilble fljB: ",B_sec_choice," Number: ",number   
              j=0
              for x in fliC_fljB_list:
                if x[0].split("_")[1] not in Sero_list_B:
                  if j<=x[1]:
                    j=x[1]
                    B_third_choice=x[0].split("_")[1]
                    number=x[1]
              if locals().has_key('B_third_choice'):
                Sero_list_B.append(B_third_choice)
                print "$$$Third possilble fljB: ",B_third_choice," Number: ",number 
          except:
            print "$$$No fljB, or failure of mapping (second run)"


def assembly(type,sra_name,for_fq,rev_fq,for_sai,rev_sai,sam,bam,database,database2,list_length,mapping_mode):
  os.system("bwa index database/"+database)
  if file_mode!="3":
    if mapping_mode=="mem":
      os.system("bwa mem database/"+database+" "+for_fq+" "+rev_fq+" > "+sam) #2014/12/23
    elif mapping_mode=="sam":
      os.system("bwa aln database/"+database+" "+for_fq+" > "+for_sai)
      os.system("bwa aln database/"+database+" "+rev_fq+" > "+rev_sai)
      os.system("bwa sampe database/"+database+" "+for_sai+" "+ rev_sai+" "+for_fq+" "+rev_fq+" > "+sam)
    elif mapping_mode=="nanopore":
      os.system("bwa mem -x ont2d database/"+database+" "+for_fq+" "+rev_fq+" > "+sam) 
  else:
    if mapping_mode=="mem":
      os.system("bwa mem database/"+database+" "+for_fq+" > "+sam) #2014/12/23
    elif mapping_mode=="sam":
      os.system("bwa aln database/"+database+" "+for_fq+" > "+for_sai)
      os.system("bwa samse database/"+database+" "+for_sai+" "+for_fq+" > "+sam)
    elif mapping_mode=="nanopore":
      os.system("bwa mem -x ont2d database/"+database+" "+for_fq+" > "+sam)
  
  os.system("samtools view -F 4 -Sbh "+sam+" > "+bam)
  os.system("samtools view -h -o "+sam+" "+bam)
  os.system("cat "+sam+"|awk '{if ($5>=0) {print $10}}'>"+database+sam+"_seq.txt")
  os.system("cat "+sam+"|awk '{if ($5>=0) {print $1}}'>"+database+sam+"_title.txt")
  file1=open(database+sam+"_title.txt","r")
  file2=open(database+sam+"_seq.txt","r")
  file1=file1.readlines()
  file2=file2.readlines()
  file=open(database+"_"+sam+".fasta","w")
  for i in range(len(file1)):
    title=">"+file1[i]
    seq=file2[i]
    if len(seq)>=50 and len(title)>6:#generally, can be adjusted with different situations
      file.write(title)
      file.write(seq)
  file.close()
  os.system("mv "+database+"_"+sam+".fasta"+" database")######1/27/2015
  os.system('makeblastdb -in database/'+database+"_"+sam+".fasta"+' -out '+sam+".fasta"+'_db '+'-dbtype nucl >temp.txt') #temp.txt is to forbid the blast result interrupt the output of our program###1/27/2015
  os.system("blastn -query database/"+database2+" -db "+sam+".fasta"+"_db -out "+database2+"_vs_"+sam+".xml -outfmt 5 >temp.txt")###1/27/2015
  handle=open(database2+"_vs_"+sam+".xml")
  handle=NCBIXML.parse(handle)
  handle=list(handle)
  List=[]
  List_score=[]
  for i in range(len(handle)):
    if len(handle[i].alignments)>0:
      List.append(handle[i].query)
      score=0
      for j in range(len(handle[i].alignments)):
        for z in range(len(handle[i].alignments[j].hsps)):
            score+=handle[i].alignments[j].hsps[z].bits
      List_score.append(score)
  temp=dict(zip(List,List_score))
  Final_list=sorted(temp.iteritems(), key=lambda d:d[1], reverse = True)
  family=database2.split("_")[2]
  try:
    Final_list[0][0].split("_")[1] # or it will always print "$$$Genome...."(next line)
    print "$$$Genome:",sra_name
    print "$$$Most possilble "+type+": ",Final_list[0][0].split("_")[1]," Score(due_to_special_test, number changed to score): ",Final_list[0][1]
    print Final_list
  except:
    if type=="fliC":
      print "$$$There may be no hit for "+type+"_"+family+" family due to the reads not covering core seqeunce, but just based on reads hit number, the most possible one is: ",fliC_option
    if type=="fljB":
      print "$$$There may be no hit for "+type+"_"+family+" family due to the reads not covering core seqeunce, but just based on reads hit number, the most possible one is: ",fljB_option
  os.system("rm "+database2+"_vs_"+sam+".xml")###01/28/2015
  os.system("rm "+database+sam+"_seq.txt")###01/28/2015
  os.system("rm "+database+sam+"_title.txt")###01/28/2015
  os.system("rm temp.txt")###01/28/2015
  os.system("rm "+sam+".fasta"+"_db.*")###01/28/2015

z=0
target=sys.argv[1] #should be sra format
data_base=sys.argv[2]
mapping_mode=sys.argv[3]
if sys.argv[4] not in ("1","2","3"):
  additional_file=sys.argv[4]
  file_mode=sys.argv[5]
else:
  additional_file=""
  file_mode=sys.argv[4]
BWA_analysis(target,additional_file,data_base,mapping_mode,file_mode,z)



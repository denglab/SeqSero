#!/usr/bin/env python



import os
from Bio import SeqIO
import sys
import itertools
from Initial_Conditions import phase1
from Initial_Conditions import phase2
from Initial_Conditions import phaseO
from Initial_Conditions import sero
import time
import multiprocessing
import string

#m=string.atoi(sys.argv[1])
m=1 #temperorily, m can be set as one, because we just need one core to deal with it
file_name=sys.argv[1]

def Combine(b,c):
  fliC_combinations=[]
  fliC_combinations.append(",".join(c))
  temp_combinations=[]
  for i in range(len(b)):
    for x in itertools.combinations(b,i+1):
      temp_combinations.append(",".join(x))
  for x in temp_combinations:
    temp=[]
    for y in c:
      temp.append(y)
    temp.append(x)
    temp=",".join(temp)
    temp=temp.split(",")
    temp.sort()
    temp=",".join(temp)
    fliC_combinations.append(temp)
  return fliC_combinations


def Test(file1,z,q):
  fliC="?"
  fljB="?"
  Otype="?"
  oafA=""#$$$$
  O3_10=""
  O1_3_19=""
  file2=file1.replace(' ','_').replace(":","__").replace("[","").replace("]","")
  try:
    os.rename(file1, file2)
    real_file=file2
  except:
    real_file=file1
  #print "###The genome name:",file1
  dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))###01/27/2015
  os.system('touch result.txt')
  database_path="database"###01/27/2015
  os.system('python '+dirpath+'/Otype_determine_analysis.py '+database_path+'/Typhimurium_LT2_gnd_galF.fasta '+real_file+' '+database_path+'/new_Oserotype.fasta >temp_result_'+str(q)+'O.txt')
  os.system('cat temp_result_'+str(q)+'O.txt>>data_log.txt')
  handle=open('temp_result_'+str(q)+'O.txt',"r")
  handle=handle.readlines()
  for line in handle:
    if "$$$ Most" in line and "O_type" in line:
      Otype=line.split("O-")[1].split("_")[0].split(" ")[0]
      Otype=Otype.replace("\n","").strip()
      #print line,
    elif "$$$ No" in line:
      Otype="-"
      if "O-9" in line:
        Otype="9"
      #print line,
    elif "$$$O5-" in line:#$$$
      oafA="-"
    elif "$$$O3,10 more possible" in line:#$$$
      O3_10="+"
    elif "$$$O1,3,19 more possible" in line:
      O1_3_19="+"
  if Otype=="1,3,19" or Otype=="3,10":#$$$judge O3,10 before formula forms
    if O3_10=="+":
      Otype="3,10"
    elif O1_3_19=="+":
      Otype="1,3,19"
    else:
      print "No_O3,10_O1,3,19_spe_sequences"
  os.system('python '+dirpath+'/H_combination_output_analysis.py '+real_file+' '+database_path+'/H_new_fliC_protein_database.fasta '+database_path+'/H_new_fljB_protein_database.fasta >temp_result_'+str(q)+'H.txt')
  os.system('cat temp_result_'+str(q)+'H.txt>>data_log.txt')
  handle2=open('temp_result_'+str(q)+'H.txt',"r")
  handle2=handle2.readlines()
  suspect="no" #for the first choice doesn't hit core sequence
  for line in handle2:
    if "$$$ Most" in line and "fliC" in line:
      #print line,
      fliC=line.split("fliC_type:  ")[1].split("_")[0].strip()
      if fliC=="g,m,p,s":
        fliC="g,m,s"
    elif "$$$ No" in line and "fliC" in line:
      fliC="-"
      #print line,
    elif "$$$ Most" in line and "fljB" in line:
      #print line,
      fljB=line.split("fljB_type:  ")[1].split("_")[0].strip()
    elif "$$$ No" in line and "fljB" in line:
      fljB="-"
      #print line,
  if Otype=="9" and fliC=="g,m" and fljB=="-":
    os.system('python '+dirpath+'/special_gene_test_assemblies.py '+database_path+'/specific_genes.fasta '+real_file+' sdf >temp_result_'+str(q)+'sdf.txt')
    os.system('cat temp_result_'+str(q)+'sdf.txt>>data_log.txt')
    handle3=open('temp_result_'+str(q)+'sdf.txt',"r")
    sdf=""
    for line in handle3:
      if "$$$" in line and "got a hit" in line:
        #print line,
        sdf="+"
    if sdf!="+":
        sdf="-"

  seronames=[]
  for i in range(len(phase1)):
    fliC_combine=[]
    fljB_combine=[]
    if phaseO[i]==Otype:
      if phase1[i].count("[")==0:
        fliC_combine.append(phase1[i])
      elif phase1[i].count("[")>=1:
        c=[]
        b=[]
        if phase1[i][0]=="[" and phase1[i][-1]=="]" and phase1[i].count("[")==1:
          content=phase1[i].replace("[","").replace("]","")
          fliC_combine.append(content)
          fliC_combine.append("-")
        else:
          for x in phase1[i].split(","):
            if "[" in x:
              b.append(x.replace("[","").replace("]",""))
            else:
              c.append(x)
          fliC_combine=Combine(b,c) #Combine will offer every possible combinations of the formula, like f,[g],t: f,t  f,g,t
      if phase2[i].count("[")==0:
        fljB_combine.append(phase2[i])
      elif phase2[i].count("[")>=1:
        d=[]
        e=[]
        if phase2[i][0]=="[" and phase2[i][-1]=="]" and phase2[i].count("[")==1:
          content=phase2[i].replace("[","").replace("]","")
          fljB_combine.append(content)
          fljB_combine.append("-")
        else:
          for x in phase2[i].split(","):
            if "[" in x:
              d.append(x.replace("[","").replace("]",""))
            else:
              e.append(x)
          fljB_combine=Combine(d,e) 
      new_fliC=fliC.split(",") #because some antigen like r,[i] not follow alphabetical order, so use this one to judge and can avoid missings
      new_fliC.sort()
      new_fliC=",".join(new_fliC)
      new_fljB=fljB.split(",")
      new_fljB.sort()
      new_fljB=",".join(new_fljB)
      if (new_fliC in fliC_combine or fliC in fliC_combine) and (new_fljB in fljB_combine or fljB in fljB_combine):
        seronames.append(sero[i])
  if len(seronames)==0:
    seronames=["N/A (The predicted antigenic profile does not exist in the White-Kauffmann-Le Minor scheme)"]
  star=""
  star_line=""
  if len(seronames)>1:
    star="*"
    star_line="The predicted serotypes share the same general formula:\t"+Otype+":"+fliC+":"+fljB+"\n"##
  #print "$$$The most possible formula is: (by the order O:H1:H2) ",Otype,":",fliC,":",fljB
  #print "$$$The possible serotyes are:",seronames
  m=0
  for y in seronames:
    if y in file1:
      #print "$$$ Is the judgement true? Answer:Yes!"       #here we use file1, because we want ":", while file2 turned it to "__"
      answer="Yes"
      m=1
  if m==0:
    #print "$$$ Is the judgement true? Answer: Need to check the records and file names"
    answer="Not sure"
  print "\n","\n"
  predict_form=Otype+":"+fliC+":"+fljB
  predict_sero=(" or ").join(seronames)
  if predict_form=="9:g,m:-":#
    predict_form=predict_form+"\nSdf prediction:"+sdf #
    if sdf=="-":#
      star="*"#
      star_line="Additional characterization is necessary to assign a serotype to this strain.  Commonly circulating strains of serotype Enteritidis are sdf+, although sdf- strains of serotype Enteritidis are known to exist. Serotype Gallinarum is typically sdf- but should be quite rare. Sdf- strains of serotype Enteritidis and serotype Gallinarum can be differentiated by phenotypic profile or genetic criteria.\n"##
      predict_sero="See comments below"#
  elif predict_form=="4:i:-":#03252016#
    predict_sero="potential monophasic variant of Typhimurium"#03252016#
  elif predict_form=="4:r:-":#03252016#
    predict_sero="potential monophasic variant of Heidelberg"#03252016# 
  elif predict_form=="4:b:-":#03252016#
    predict_sero="potential monophasic variant of Paratyphi B"#03252016# 
  elif predict_form=="8:e,h:1,2":#03282016#
    predict_sero="Newport"#03282016#
    star="*"##03282016#
    star_line="Serotype Bardo shares the same antigenic profile with Newport, but Bardo is exceedingly rare."#03282016#
  claim="The serotype(s) is/are the only serotype(s) with the indicated antigenic profile currently recognized in the Kauffmann White Scheme.  New serotypes can emerge and the possibility exists that this antigenic profile may emerge in a different subspecies.  Identification of strains to the subspecies level should accompany serotype determination; the same antigenic profile in different subspecies is considered different serotypes."##
  if "N/A" in predict_sero:###added after standalone version, 2015/2/3
    claim=""###added after standalone version, 2015/2/3  
  '''
  new_file=open(file2+".txt","w")
  new_file.write(file2+"\t"+"O-"+Otype+"\t"+fliC+"\t"+fljB+"\t"+Otype+":"+fliC+":"+fljB+"\t"+(" or ").join(seronames)+"\t"+answer+"\t"+suspect+"\n")
  new_file.close()
  '''
  if "Typhimurium" in predict_sero and oafA=="-":#$$$$#03252016#
    predict_sero=predict_sero.strip()+"(O5-)"#03252016#
    star="*"#
    star_line="Detected the deletion of O5-."
  new_file=open("Seqsero_result.txt","w")
  new_file.write("Input files:\t"+file2+"\n"+"O antigen prediction:\t"+"O-"+Otype+"\n"+"H1 antigen prediction(fliC):\t"+fliC+"\n"+"H2 antigen prediction(fljB):\t"+fljB+"\n"+"Predicted antigenic profile:\t"+predict_form+"\n"+"Predicted serotype(s):\t"+predict_sero+star+"\n"+star+star_line+claim+"\n")##
  new_file.close()
  os.system("rm temp_result_"+str(q)+"*.txt")###01/28/2015
  os.system("rm result.txt")###01/28/2015
  #os.system("rm -rf database")###01/28/2015 
  os.system("rm *.fasta *.xml *.fa")###01/28/2015

    
def main():
  files1=[]
  files1.append(file_name)
  file_names=[]
  fastq_names=[]
  for file1 in files1:
    if file1[-6:]=='.fasta' or file1[-4:]=='.fna' or file1[-3:]=='.fa' or file1[-4:]=='.fsa':
      file_names.append(file1)
    if file1[-9:]==".fastq.gz" or file1[-6:]==".fastq":
      core_name=file1[:8]
      fastq_names.append(core_name)
  fastq_names=list(set(fastq_names))
  file_names=file_names+fastq_names
  for i in range(0,len(file_names),m):
    jobs=[]
    txt_names=[]
    if len(file_names)>=i+m:
      for j in range(m):
        p = multiprocessing.Process(target=Test,args=(file_names[j+i],i+j+1,i+j,))
        jobs.append(p)
        p.start()
    else:
      t=m+i-len(file_names)
      for j in range(m-t):
        p = multiprocessing.Process(target=Test,args=(file_names[j+i],i+j+1,i+j,))
        jobs.append(p)
        p.start()
        '''
    for j in xrange(len(jobs)):
      jobs[j].join()
      txt_names.append(file_names[j+i].replace(' ','_').replace(":","__").replace("[","").replace("]","")+".txt")
    print txt_names
    for j in xrange(len(txt_names)):
      print i,"and",j
      print i+j+1
      file=open(txt_names[j],"r")
      handle=list(file)
      b=handle[0].split("\t")
      print b
      sheet.write(i+j+1,0,b[0])
      sheet.write(i+j+1,1,b[1])
      sheet.write(i+j+1,2,b[2])
      sheet.write(i+j+1,3,b[3])
      sheet.write(i+j+1,4,b[4])
      sheet.write(i+j+1,5,b[5])
      sheet.write(i+j+1,6,b[6])
      sheet.write(i+j+1,7,b[7])

  print "End time,",time.time()
  file3.save("Seqsero_result2.xls")
  '''


if __name__ == '__main__':
  main()









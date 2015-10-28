#just an possible use, we can use it to replace H**.py, treat fliC and fljB as the target genes? 

from __future__ import division
import sys
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO

def special_gene(target_fie,database,gene_list):
  database=database.split("/")[-1]##########1/27/2015
  os.system('makeblastdb -in database/'+database+' -out '+database+'_db -dbtype nucl')##########1/28/2015
  os.system('blastn -query '+target_file+' -db '+database+'_db -out '+database+'_vs_'+target_file+'.xml '+'-outfmt 5')##########1/28/2015
  xml_file=database+'_vs_'+target_file+'.xml'
  result_handle=open(xml_file)
  blast_record=NCBIXML.parse(result_handle)
  records=list(blast_record)
  E_thresh=1e-10
  for x in gene_list:
    handle=SeqIO.parse("database/"+database,"fasta")##########1/28/2015
    length_list=[]
    for y in handle:
      if x in y.description:
        length_x=len(y.seq)
        length_list.append(length_x)
    aver_len=float(sum(length_list))/len(length_list)
    hspbit=[]
    alignmentlist=[]
    for record in records:
      for alignment in record.alignments:
        if x in alignment.hit_def: #multi gene database, so...
          print x,"got a hit, evaluating the hit quality..."
          score=0
          for hsp in alignment.hsps:
            if hsp.expect<E_thresh:
              score+=hsp.bits
          alignment=alignment.hit_def+':'+str(score)
          hspbit.append(score)
          alignmentlist.append(alignment)
    scorelist=dict(zip(alignmentlist,hspbit))
    score=0
    for Htype in scorelist:
      if scorelist[Htype]>score:
        First_Choice=Htype
        score=scorelist[Htype]
    if float(score)>=0.1*aver_len:
      print "$$$",First_Choice,"got a hit, score:",score
    else:
      print "$$$No ",x,"exists"
  os.system("rm "+database+"_db.*")##########1/28/2015
  os.system("rm "+xml_file)##########1/28/2015

database=sys.argv[1]
target_file=sys.argv[2]
gene_list=[]
a=1
i=3
while a==1:
  try:
    gene_list.append(sys.argv[i])
    i+=1
  except:
    a=0
special_gene(target_file,database,gene_list)
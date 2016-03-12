#!/usr/bin/env python
# "H_combination_output_analysis.py target.fasta fliCdatabase.fasta fljBdatabase.fasta"
# must have ispcr and primers of fliC and fljB at the same directory


import os
from Bio import SeqIO
import sys
from Bio.Blast import NCBIXML
from Initial_Conditions import phase1
from Initial_Conditions import phase2



target=sys.argv[1]
database_fliC=sys.argv[2]
database_fljB=sys.argv[3]
output=target.split('.')[0]+'_out.fasta'

dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))###01/27/2015
database_path="database"###01/27/2015,database_path=dirpath+"/database"
os.system(dirpath+'/isPcr maxSize=3000 tileSize=7 minPerfect=7 minGood=7 '+target+' '+dirpath+'/../primers/seq_primer_fliC.txt '+target+'_fliC.fa')
os.system(dirpath+'/isPcr maxSize=3000 tileSize=7 minPerfect=7 minGood=7 '+target+' '+dirpath+'/../primers/seq_primer_fljB.txt '+target+'_fljB.fa')
fliC=target+'_fliC.fa'
fljB=target+'_fljB.fa'

if os.path.getsize(fliC)>10:
  os.system('makeblastdb -in '+database_fliC+' -out '+database_fliC+'_db '+'-dbtype prot')###01/28/2015,no need to add fljB address,  because input is abs address already
  os.system('blastx -seg=no -query '+fliC+' -db '+database_fliC+'_db '+'-out '+'FliC_Htype_'+target+'.xml '+'-outfmt 5')
  print target
  fliC_XML='FliC_Htype_'+target+'.xml'
  fliC_handle=open(fliC_XML)
  records=NCBIXML.parse(fliC_handle)
  fliC_records=list(records)
  E_thresh=1e-10
  hspbit=[]
  alignmentlist=[]
  for record in fliC_records:
    for alignment in record.alignments:
      hsp_bit_score=0
      startlist=[]#the percentage algorithm don't consider one situation, the new hsp cover old hsp
      endlist=[]#
      for hsp in alignment.hsps:
        start=hsp.query_start#
        end=hsp.query_end#
        leng=abs(start-end)#
        if hsp.expect<E_thresh:#
          if start>end:#
            temp=start#
            start=end#
            end=start#
          if len(startlist)==0:#
            hsp_bit_score=hsp_bit_score+hsp.bits#
            startlist.append(start)#
            endlist.append(end)#
          else:#
            for i in range(len(startlist)):#
              if startlist[i]<start<endlist[i]:#
                start=endlist[i]+1#
              if startlist[i]<end<endlist[i]:#03112016
                end=startlist[i]-1#
            if end<start:#the new hsp was included in old hsp#
              percentage=0#
            else:
              percentage=float(end-start)/leng#
              startlist.append(start)#
              endlist.append(end)#
            hsp_bit_score=hsp_bit_score+percentage*hsp.bits#          
      alignment=alignment.hit_def+':'+str(hsp_bit_score)
      hspbit.append(hsp_bit_score)
      alignmentlist.append(alignment)
  scorelist=dict(zip(alignmentlist,hspbit))
  score=0
  serotype=[]
  seroscore=[]
  for Htype in scorelist:
    if scorelist[Htype]>score:
      First_Choice=Htype
      score=scorelist[Htype]
  if  locals().has_key('First_Choice'):
    serotype.append(First_Choice.split("__")[0])
    seroscore.append(score)
    secscore=0
    for Htype in scorelist:

      if scorelist[Htype]>secscore and (Htype.split("__")[0] not in serotype):
        Sec_Choice=Htype
        secscore=scorelist[Htype]
    if locals().has_key('Sec_Choice'):  
      serotype.append(Sec_Choice.split("__")[0])
      seroscore.append(secscore)
      thirdscore=0
      for Htype in scorelist:
        if scorelist[Htype]>thirdscore and (Htype.split("__")[0] not in serotype):
          Third_Choice=Htype
          thirdscore=scorelist[Htype]
      if locals().has_key('Third_Choice'):
        serotype.append(Third_Choice.split("__")[0])
        seroscore.append(thirdscore)
  print serotype,seroscore
  if score>100:
    print '#',target,'$$$ Most possible H_fliC_type: ',First_Choice,'\n'
    print '$$$ bit_score:',score,'\n'
    if locals().has_key('secscore'):
      if secscore>100:
        print '#',target,'$$$ Second possible H_fliC_type: ',Sec_Choice,'\n'
        print '$$$ Second bit_score:',secscore,'\n'
        if locals().has_key('thirdscore'):
          if thirdscore>100:
            print '#',target,'$$$ Third possible H_fliC_type: ',Third_Choice,'\n'
            print '$$$ Third bit_score:',thirdscore,'\n'
  else:
    print '$$$ No fliC in',target
else:
  score=1
  print '$$$ No fliC (no file created) in',target



if os.path.getsize(fljB)>10:
  os.system('makeblastdb -in '+database_fljB+' -out '+database_fljB+'_db '+'-dbtype prot')###01/28/2015,no need to add fljB address,  because input is abs address already
  os.system('blastx -query '+fljB+' -db '+database_fljB+'_db '+'-out '+'FljB_Htype_'+target+'.xml '+'-outfmt 5')
  print target
  fljB_XML='FljB_Htype_'+target+'.xml'
  fljB_handle=open(fljB_XML)
  records=NCBIXML.parse(fljB_handle)
  fljB_records=list(records)
  E_thresh=1e-10
  hspbit=[]
  alignmentlist=[]
  for record in fljB_records:
    for alignment in record.alignments:
      hsp_bit_score=0
      startlist=[]#
      endlist=[]#
      for hsp in alignment.hsps:
        start=hsp.query_start#
        end=hsp.query_end#
        leng=abs(start-end)#
        if hsp.expect<E_thresh:#
          if start>end:#
            temp=start#
            start=end#
            end=start#
          if len(startlist)==0:#
            hsp_bit_score=hsp_bit_score+hsp.bits#
            startlist.append(start)#
            endlist.append(end)#
          else:#
            for i in range(len(startlist)):#
              if startlist[i]<start<endlist[i]:#
                start=endlist[i]+1#
              if startlist[i]<end<endlist[i]:#03112016
                end=startlist[i]-1#
            if end<start:#the new hsp was included in old hsp#
              percentage=0#
            else:
              percentage=float(end-start)/leng#
              startlist.append(start)#
              endlist.append(end)#
            hsp_bit_score=hsp_bit_score+percentage*hsp.bits#   
      alignment=alignment.hit_def+':'+str(hsp_bit_score)
      hspbit.append(hsp_bit_score)
      alignmentlist.append(alignment)
  fljB_scorelist=dict(zip(alignmentlist,hspbit))

  fljB_score=0
  fljB_serotype=[]
  fljB_seroscore=[]
  for Htype in fljB_scorelist:
    if fljB_scorelist[Htype]>fljB_score:
      fljB_First_Choice=Htype
      fljB_score=fljB_scorelist[Htype]
  if  locals().has_key('fljB_First_Choice'):
    fljB_serotype.append(fljB_First_Choice.split("__")[0])
    fljB_seroscore.append(fljB_score)
    fljB_secscore=0
    for Htype in fljB_scorelist:
      if fljB_scorelist[Htype]>fljB_secscore and (Htype.split("__")[0] not in fljB_serotype):
        fljB_Sec_Choice=Htype
        fljB_secscore=fljB_scorelist[Htype]
    if locals().has_key('fljB_Sec_Choice'):  
      fljB_serotype.append(fljB_Sec_Choice.split("__")[0])
      fljB_seroscore.append(fljB_secscore)
      fljB_thirdscore=0
      for Htype in fljB_scorelist:
        if fljB_scorelist[Htype]>fljB_thirdscore and (Htype.split("__")[0] not in fljB_serotype):
          fljB_Third_Choice=Htype
          fljB_thirdscore=fljB_scorelist[Htype]
      if locals().has_key('fljB_Third_Choice'):
        fljB_serotype.append(fljB_Third_Choice.split("__")[0])
        fljB_seroscore.append(fljB_thirdscore)

  if fljB_score>100:
    print '#',target,'$$$ Most possible H_fljB_type: ',fljB_First_Choice,'\n'
    print '$$$ Most bit_score:',fljB_score,'\n'
    if locals().has_key('fljB_secscore'):
      if fljB_secscore>100:
        print '#',target,'$$$ Second possible H_fljB_type: ',fljB_Sec_Choice,'\n'
        print '$$$ Second bit_score:',fljB_secscore,'\n'
        if locals().has_key('fljB_thirdscore'):
          if fljB_thirdscore>100:
            print '#',target,'$$$ Third possible H_fljB_type: ',fljB_Third_Choice,'\n'
            print '$$$ Third bit_score:',fljB_thirdscore,'\n'
  else:
    print '$$$ No fljB in',target
else:
  fljB_score=1
  print '$$$ No fljB (no file created) in',target


if score>100 and fljB_score>100:
  fliC_sero=dict(zip(serotype,seroscore))
  fljB_sero=dict(zip(fljB_serotype,fljB_seroscore))
  combination=[]
  combination_score=[]
  for seroname in fliC_sero:
    for fljB_seroname in fljB_sero:
      for i in range(len(phase1)):
        if phase1[i]==seroname and phase2[i]==fljB_seroname:
          name=seroname+"_"+fljB_seroname
          score=fliC_sero[seroname]+fljB_sero[fljB_seroname]
          combination.append(name)
          combination_score.append(score)
  combinationlist=dict(zip(combination,combination_score))  #we can do the filteration here
  final_dict=sorted(combinationlist.iteritems(), key=lambda d:d[1], reverse = True)
  print "$$_H:Order:",final_dict
elif score>100 and fljB_score<100:
  print "$$_H:No fljB, only fliC, and its order:",First_Choice,Sec_Choice,Third_Choice
elif score<100 and fljB_score>100:
  print "$$_H:No fliC, only fljB, and its order:",fljB_First_Choice,fljB_Sec_Choice,fljB_Third_Choice
elif score==1 and fljB_score>100:
  print "$$_H:No fliC (file) existed, only fljB, and its order:",fljB_First_Choice,fljB_Sec_Choice,fljB_Third_Choice
elif score==1 and fljB_score<100:
  print "$$_H:No fliC (file) existed, and no fljB"
elif score>100 and fljB_score==1:
  print "$$_H:No fljB (file) existed, only fliC, and its order:",First_Choice,Sec_Choice,Third_Choice
elif score<100 and fljB_score==1:
  print "$$_H:No fljB (file) existed, and no fliC"
else:
  print "$$_H:No fliC and fljB"


'''
E_thresh=1e-10
hspbit=[]
alignmentlist=[]
for record in fliC_records:
  for alignment in record.alignments:
    hsp_bit_score=0
    for hsp in alignment.hsps:
      if hsp.expect<E_thresh:
        hsp_bit_score=hsp_bit_score+hsp.bits
    alignment=alignment.hit_def+':'+str(hsp_bit_score)
    hspbit.append(hsp_bit_score)
    alignmentlist.append(alignment)

scorelist=dict(zip(alignmentlist,hspbit))
score=0
for Htype in scorelist:
  if scorelist[Htype]>score:
    First_Choice=Htype
    score=scorelist[Htype]
if score>100:
  print '#',target,'Most possible H_fliC_type: ',First_Choice,'\n'
  print '#bit_score:',score,'\n'
else:
  print '#No fliC in',target


E_thresh=1e-10
hspbit=[]
alignmentlist=[]
for record in fljB_records:
  for alignment in record.alignments:
    hsp_bit_score=0
    for hsp in alignment.hsps:
      if hsp.expect<E_thresh:
        hsp_bit_score=hsp_bit_score+hsp.bits
    alignment=alignment.hit_def+':'+str(hsp_bit_score)
    hspbit.append(hsp_bit_score)
    alignmentlist.append(alignment)

scorelist=dict(zip(alignmentlist,hspbit))
fljB_score=0
for Htype in scorelist:
  if scorelist[Htype]>fljB_score:
    First_Choice=Htype
    fljB_score=scorelist[Htype]
if fljB_score>100:
  print '#',target,'Most possible H_fljB_type: ',First_Choice,'\n'
  print '#bit_score:',fljB_score,'\n'
else:
  print '#No fljB in',target
'''


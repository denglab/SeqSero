#!/usr/bin/env python

'''
revised at 11/24/2013
note:
The .py file mainly designed to deal with 5 different situations: 1.the gnd and galF on same contig; 2.the two genes on different contigs;
3.the two genes on different contigs, and gnd gene is splitted by two contigs; 4.the two genes on different contigs, and galF gene is splitted by two contigs;
5.the two genes on different contigs, and the two genes are both splitted by two contigs respectively;

input:
1) <gnd_galF_sequence>: a fasta file containing galF and gnd sequences
2) <target>: a fasta file containing the genome (may be seperated contigs), mainly as the blast database for <gnd_galF_sequence> to extract rfb region
3) <database>: a fasta file containing rfb regions of 46 different O-serotypes

output: the region between the galF and gnd in the target(blatdb) sequebce (the rfb region, assuming the rfb region is flanked by the two genes)
synopsis: ./Otype_determine_analysis.py <gnd_galF_sequence> <target> <database>
'''
from __future__ import division
import sys
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO


def test_O29(subdatabase):
  first_hspbit=[]
  first_alignmentlist=[]
  global Choice,Choice_score,O_9_score,O_2_score
  O_9_score=0
  O_2_score=0
  os.system('makeblastdb -in database/'+subdatabase+' -out '+subdatabase+'_db -dbtype nucl')###01/28/2015
  os.system('blastn -query '+target+' -db '+subdatabase+'_db -out '+subdatabase+'_vs_'+target+'.xml '+'-outfmt 5')###01/28/2015
  xml_file=subdatabase+'_vs_'+target+'.xml'
  result_handle=open(xml_file)
  blast_record=NCBIXML.parse(result_handle)
  records=list(blast_record)
  if subdatabase=="tyr_of_O2_O9.fasta":
    for record in records:    #there are many records (i.e. the '>' in query file), so change another method
      for alignment in record.alignments:
        if 'O-9' in alignment.hit_def:
          for hsp in alignment.hsps:
            if hsp.expect<E_thresh:
              O_9_score=O_9_score+hsp.bits
        elif 'O-2' in alignment.hit_def:
          for hsp in alignment.hsps:
            if hsp.expect<E_thresh:
              O_2_score=O_2_score+hsp.bits
    if O_9_score>100:
      if O_9_score>O_2_score:
        print '$$$ Most possible O_type: O-9','\n'
        print '$$$ longest_bit_score:',O_9_score,'\n'
      else:
        print '$$$ Most possible O_type: O-2','\n'
        print '$$$ longest_bit_score:',O_2_score,'\n'
    else:
      print "Assumpition wrong, no O2 or O9, return to re-analysis"
      if ('O-2_' not in Sec_Choice) and ('O-9_' not in Sec_Choice):
        print '$$$ Most possible O_type Choice (no tyr difference):',Sec_Choice
      if (('O-2_' in Sec_Choice) or ('O-9_' in Sec_Choice)) and ('O-2_' not in Third_Choice) and ('O-9_' not in Third_Choice):
          print '$$$ Most possible O_type Choice (no tyr difference):',Third_Choice
    os.system("rm tyr_of_O2_O9.fasta_db.*")###01/28/2015
    os.system("rm "+xml_file)###01/28/2015
  if subdatabase=="oafA_of_O4_O5.fasta":
    try:
      for record in records:    #there are many records (i.e. the '>' in query file), so change another method
        for alignment in record.alignments:
          if 'O-4_full' in alignment.hit_def:
            for hsp in alignment.hsps:
              if hsp.expect<E_thresh:
                O_9_score=O_9_score+hsp.bits
          elif 'O-4_5-' in alignment.hit_def:
            for hsp in alignment.hsps:
              if hsp.expect<E_thresh:
                O_2_score=O_2_score+hsp.bits
      if O_9_score>100:
        if O_9_score>O_2_score:
          print '$$$O5_none_7_base_deletion','\n'
          print '$$$ longest_bit_score:',O_9_score,'\n'
        else:
          print '$$$O5-','\n'
          print '$$$ longest_bit_score:',O_2_score,'\n'
      else:
        print '$$$O5_none_7_base_deletion,unsure','\n'
      os.system("rm oafA_of_O4_O5.fasta_db.*")###01/28/2015
      os.system("rm "+xml_file)###01/28/2015
    except:
      print "No oafA genes"
  if subdatabase=="O_3,10_and_1,3,19_spe.fasta":
    try:
      for record in records:    #there are many records (i.e. the '>' in query file), so change another method
        for alignment in record.alignments:
          if 'O_3,10_not_in_1,3,19' in alignment.hit_def:
            for hsp in alignment.hsps:
              if hsp.expect<E_thresh:
                O_9_score=O_9_score+hsp.bits
          elif 'O_1,3,19_not_in_3,10' in alignment.hit_def:
            for hsp in alignment.hsps:
              if hsp.expect<E_thresh:
                O_2_score=O_2_score+hsp.bits
      if O_9_score>200:
        if O_9_score>O_2_score:
          print '$$$O3,10 more possible','\n'
          print '$$$ longest_bit_score:',O_9_score,'\n'
      else:
        if O_2_score>100:
          print '$$$O1,3,19 more possible','\n'
          print '$$$ longest_bit_score:',O_2_score,'\n'
      os.system("rm O_3,10_and_1,3,19_spe.fasta_db.*")###01/28/2015
      os.system("rm "+xml_file)###01/28/2015
    except:
      print "No O3,10_and_O1,3,19 spe sequences"


def show_result():
        global First_Choice
        global Sec_Choice
        global Third_Choice
        hspbit=[]
        alignmentlist=[]
        for record in records:
          for alignment in record.alignments:
            hsp_bit_score=0
            startlist=[]
            endlist=[]
            for hsp in alignment.hsps:
              if hsp.expect<E_thresh:
                start=hsp.query_start
                end=hsp.query_end
                if len(startlist)==0:
                  startlist.append(start)
                  endlist.append(end)
                  hsp_bit_score=hsp_bit_score+hsp.bits
                elif len(startlist)==1:
                  if end<=startlist[0] or start>=endlist[0]:
                    startlist.append(start)
                    endlist.append(end)
                    hsp_bit_score=hsp_bit_score+hsp.bits
                  if start>=startlist[0] and end<=endlist[0]:
                    print hsp_bit_score,'Fully overlapped hsp, ignore',hsp_bit_score
                    continue
                  if start<startlist[0] and startlist[0]<end<endlist[0]:
                    print hsp_bit_score,'Partially overlapped hsp, trying to extract unoverlap part',hsp_bit_score
                    newend=startlist[0]-1
                    percentage=(newend-start)/(end-start)
                    startlist.append(start)
                    endlist.append(newend)
                    hsp_bit_score=hsp_bit_score+percentage*hsp.bits
                  if end>endlist[0] and startlist[0]<start<endlist[0]:
                    print hsp_bit_score,'Partially overlapped hsp, trying to extract unoverlap part',hsp_bit_score
                    newstart=endlist[0]+1
                    percentage=(end-newstart)/(end-start)
                    startlist.append(newstart)
                    endlist.append(end)
                    hsp_bit_score=hsp_bit_score+percentage*hsp.bits
                elif len(startlist)>1:
                  startlist.sort()
                  endlist.sort()
                  length=len(startlist)
                  if end<=startlist[0] or start>=endlist[length-1]:
                    startlist.append(start)
                    endlist.append(end)
                    hsp_bit_score=hsp_bit_score+hsp.bits
                  if start<startlist[0] and startlist[0]<end<endlist[0]:
                    print hsp_bit_score,'Partially overlapped hsp, trying to extract unoverlap part',hsp_bit_score
                    newend=startlist[0]-1
                    percentage=(newend-start)/(end-start)
                    startlist.append(start)
                    endlist.append(newend)
                    hsp_bit_score=hsp_bit_score+hsp.bits*percentage
                  if end>endlist[length-1] and startlist[length-1]<start<endlist[length-1]:
                    print hsp_bit_score,'Partially overlapped hsp, trying to extract unoverlap part',hsp_bit_score
                    newstart=endlist[length-1]+1
                    percentage=(end-newstart)/(end-start)
                    startlist.append(newstart)
                    endlist.append(end)
                    hsp_bit_score=hsp_bit_score+hsp.bits*percentage
                  else:
                    for i in range(0,length-1):         #use two i and i+1 as basic combinations
                      if (start>=startlist[i] and end<=endlist[i]) or (start>=startlist[i+1] and end<=endlist[i+1]):
                        print hsp_bit_score,'Fully overlapped hsp, ignore',hsp_bit_score
                        continue
                      if endlist[i]<=start<end<=startlist[i+1]:
                        startlist.append(start)
                        endlist.append(end)
                        hsp_bit_score=hsp_bit_score+hsp.bits
                      if startlist[i]<=start<=endlist[i] and endlist[i]<end<startlist[i+1]:
                        print hsp_bit_score,'Partially overlapped hsp, trying to extract unoverlap part',hsp_bit_score 
                        newstart=endlist[i]+1
                        percentage=(end-newstart)/(end-start)
                        startlist.append(newstart)
                        endlist.append(end)
                        hsp_bit_score=hsp_bit_score+hsp.bits*percentage
                      if endlist[i]<start<startlist[i+1] and startlist[i+1]<=end<=endlist[i+1]:
                        print hsp_bit_score,'Partially overlapped hsp, trying to extract unoverlap part',hsp_bit_score
                        newend=startlist[i+1]-1
                        percentage=(newend-start)/(end-start)
                        startlist.append(start)
                        endlist.append(newend)
                        hsp_bit_score=hsp_bit_score+hsp.bits*percentage
            alignment=alignment.hit_def
            hspbit.append(hsp_bit_score)
            alignmentlist.append(alignment)

        scorelist=dict(zip(alignmentlist,hspbit))
        score=0
        for Otype in scorelist:
          if scorelist[Otype]>score:
            First_Choice=Otype
            score=scorelist[Otype]
                
        secscore=0
        for Otype in scorelist:
          if scorelist[Otype]>secscore and Otype!=First_Choice:
            Sec_Choice=Otype
            secscore=scorelist[Otype] 
        
        thirdscore=0
        for Otype in scorelist:
          if scorelist[Otype]>thirdscore and Otype!=First_Choice and Otype!=Sec_Choice:
            Third_Choice=Otype
            thirdscore=scorelist[Otype]

        if thirdscore==0:
          names=First_Choice+Sec_Choice
        else:
          names=First_Choice+Sec_Choice+Third_Choice
        if secscore==0:
          names=First_Choice
        else:
          names=First_Choice+Sec_Choice
        if score==0:
          print "$$$ No O_type, due to no hit of rfb"
          names=""
        if 'O-2_' in names and 'O-9_' in names and ('O-2_' in First_Choice or 'O-9_' in First_Choice):
          print '#Contain O2 and O9, so change to special test'
          test_O29("tyr_of_O2_O9.fasta")
        else:
          if score>0:
            print '$$$ Most possible O_type: ',First_Choice,'\n'
            print '$$$ Most bit_score:',score,'\n'
            if "O-4_" in First_Choice:#$$$$$$$
              test_O29("oafA_of_O4_O5.fasta")#$$$$$$$
            if "O-1,3,19" in First_Choice or "O-3,10" in First_Choice:
              test_O29("O_3,10_and_1,3,19_spe.fasta")#$$$$$$$
          if secscore>0:
            print '$$$ Second possible O_type: ',Sec_Choice,'\n'
            print '$$$ Second bit_score:',secscore,'\n'
          if thirdscore>0:
            print '$$$ Third possible O_type: ',Third_Choice,'\n'
            print '$$$ Third bit_score:',thirdscore,'\n'







queries=sys.argv[1]
target=sys.argv[2]
database=sys.argv[3]
output=target.split('.')[0]+'_out.fa'
print "$$:",target


os.system('makeblastdb -in '+target+' -out '+target+'_db '+'-dbtype nucl')###01/28/2015
os.system('blastn -query '+queries+' -db '+target+'_db '+'-out '+queries+'_vs_'+target+'.xml '+'-outfmt 5')###01/28/2015, since it's abs address for "run_auto*.py", so no need to change "query" address this time
xml_file=queries+'_vs_'+target+'.xml'
print '\n'
result_handle=open(xml_file)
blast_record=NCBIXML.parse(result_handle)
blast_record=list(blast_record)
target_seq=SeqIO.parse(target,'fasta')
target_seq=list(target_seq)
#os.system("rm "+xml_file)###01/28/2015
#os.system("rm "+target+'_db.*')###01/28/2015

if len(blast_record)==2:
  print 'Hits have been got'+'\n'
  if len(blast_record[0].alignments)==1 and len(blast_record[1].alignments)==1:
    print 'Checking the number of alignments:  2 alignments obtained'+'\n'
    if len(blast_record[0].alignments[0].hsps)==1 and len(blast_record[1].alignments[0].hsps)==1:
      print 'Checking the number of hsps:  each alignment has 1 hsp'+'\n'
      if blast_record[0].alignments[0].hit_def==blast_record[1].alignments[0].hit_def:
        print 'Checking locations of hits:  Both hits are located in '+'"'+str(blast_record[0].alignments[0].hit_def)+'"'+'...'+'\n'
               
        hit_1_start=blast_record[0].alignments[0].hsps[0].sbjct_start
        hit_1_end=blast_record[0].alignments[0].hsps[0].sbjct_end

        hit_2_start=blast_record[1].alignments[0].hsps[0].sbjct_start
        hit_2_end=blast_record[1].alignments[0].hsps[0].sbjct_end

        if hit_1_start>hit_1_end:
          buffer=hit_1_start
          hit_1_start=hit_1_end
          hit_1_end=buffer

        if hit_2_start>hit_2_end:
          buffer=hit_2_start
          hit_2_start=hit_2_end
          hit_2_end=buffer

        print 'hit_1_start: '+str(hit_1_start)
        print 'hit_1_end: '+str(hit_1_end)

        print 'hit_2_start: '+str(hit_2_start)
        print 'hit_2_end: '+str(hit_2_end)

    
        if hit_1_end<hit_2_start:                    
          extract_start=hit_1_end+1                     
          extract_end=hit_2_start-1                        
        else:                                       
          extract_start=hit_2_end+1
          extract_end=hit_1_start-1
                  
        print 'start: '+str(extract_start), 'end: '+str(extract_end)+'\n'

        for contig in target_seq:
          if (contig.description==blast_record[0].alignments[0]) or (contig.description.replace(" ","")==blast_record[0].alignments[0].hit_def.replace(" ","")):
            target_contig=contig

        rfb_region=target_contig[extract_start:extract_end]
        print 'Extracted rfb region length:  '+str(len(rfb_region.seq.tostring()))+'\n'
        print 'Extracted rfb region saved in:  '+output+'\n'

        outfile=open(output,'w')
        title='>'+target.split('.')[0]+' rfb region:'+blast_record[0].alignments[0].hit_def+':'+str(extract_start)+' to '+str(extract_end)+'_'+str(len(rfb_region.seq.tostring()))+'bp'+')'
        outfile.write(title)  
        outfile.write('\n')
        outfile.write(rfb_region.seq.tostring())
        outfile.close()

               
        os.system('makeblastdb -in '+database+' -out '+database+'_db '+'-dbtype nucl')
        os.system('blastn -query '+output+' -db '+database+'_db '+'-out '+'Blast_Otype_'+target+'.xml '+'-outfmt 5')
        xml_file='Blast_Otype_'+target+'.xml'
        print '\n'


        filehandle=open(xml_file)
        records=NCBIXML.parse(filehandle)
        records=list(records)
        target_hspbit=0      #give a original value to target_hspbit
        second_hspbit=0
        third_hspbit=0
        E_thresh=1e-10
        show_result()

      else:
        print 'Checking locations of hits:  the two hits are not located in same contig......'+'\n'
        hit_1_start=blast_record[0].alignments[0].hsps[0].sbjct_start
        hit_1_end=blast_record[0].alignments[0].hsps[0].sbjct_end
        hit_2_start=blast_record[1].alignments[0].hsps[0].sbjct_start
        hit_2_end=blast_record[1].alignments[0].hsps[0].sbjct_end

        if hit_1_start>hit_1_end:
          buffer=hit_1_start
          hit_1_start=hit_1_end
          hit_1_end=buffer

        if hit_2_start>hit_2_end:
          buffer=hit_2_start
          hit_2_start=hit_2_end
          hit_2_end=buffer

        outfile=open(output,'w')
        for contig in target_seq:
          if contig.description==blast_record[0].alignments[0].hit_def:
            potentialsequence1=contig.seq[0:hit_1_start-1].tostring()
            potentialsequence2=contig.seq[hit_1_end+1:].tostring()
            title1='>'+target.split('.')[0]+' potential rfb region:'+blast_record[0].alignments[0].hit_def+':0 to '+str(hit_1_start-1)+'_total:'+str(len(potentialsequence1))+'bp'
            outfile.write(title1)  
            outfile.write('\n')
            outfile.write(potentialsequence1)
            outfile.write('\n')
            title2='>'+target.split('.')[0]+' potential rfb region:'+blast_record[0].alignments[0].hit_def+':'+str(hit_1_end+1)+' to contig end'+'_total:'+str(len(potentialsequence2))+'bp'
            outfile.write(title2)  
            outfile.write('\n')
            outfile.write(potentialsequence2)
            outfile.write('\n')  
          elif contig.description==blast_record[1].alignments[0].hit_def:
            potentialsequence1=contig.seq[0:hit_2_start-1].tostring()
            potentialsequence2=contig.seq[hit_2_end+1:].tostring()
            title1='>'+target.split('.')[0]+' potential rfb region:'+blast_record[1].alignments[0].hit_def+':0 to '+str(hit_2_start-1)+'_total:'+str(len(potentialsequence1))+'bp'
            outfile.write(title1)  
            outfile.write('\n')
            outfile.write(potentialsequence1)
            outfile.write('\n')
            title2='>'+target.split('.')[0]+' potential rfb region:'+blast_record[1].alignments[0].hit_def+':'+str(hit_2_end+1)+' to contig end'+'_total:'+str(len(potentialsequence2))+'bp'
            outfile.write(title2)  
            outfile.write('\n')
            outfile.write(potentialsequence2)
            outfile.write('\n')

    
        outfile.close()

        os.system('makeblastdb -in '+database+' -out '+database+'_db '+'-dbtype nucl')
        os.system('blastn -query '+output+' -db '+database+'_db '+'-out '+'Blast_Otype_'+target+'.xml '+'-outfmt 5')
        xml_file='Blast_Otype_'+target+'.xml'
        print '\n'


        filehandle=open(xml_file)
        records=NCBIXML.parse(filehandle)
        records=list(records)
        realrecord1=records[0]
        if len(records[1].alignments)>len(records[0].alignments):
          realrecord1=records[1]
        realrecord2=records[2]
        if len(records[3].alignments)>len(records[2].alignments):
          realrecord2=records[3]
        title='>'+records[0].query.split(':')[0]+'_combined_potential_sequences'
        filecontent=SeqIO.parse(output,'fasta')
        filecontent=list(filecontent)
        for contig in filecontent:
          if contig.description==realrecord1.query:
            sequence=contig.seq.tostring()

        for contig in filecontent:
          if contig.description==realrecord2.query:
            sequence=sequence+contig.seq.tostring()

        outfile=open('combined_sequence.fasta','w')
        outfile.write(title)
        outfile.write('\n')
        outfile.write(sequence)
        outfile.close()
                
        os.system('blastn -query combined_sequence.fasta'+' -db '+database+'_db '+'-out '+'Combined_seq_blast_'+target+'.xml '+'-outfmt 5')
        xml_file='Combined_seq_blast_'+target+'.xml'
        print '\n'

        filehandle=open(xml_file)
        records=NCBIXML.parse(filehandle)
        records=list(records)
        E_thresh=1e-10
                
        show_result()
    
        

    else:
      print '$$$ No O_type result, please check the number of hsps:  some alignment have more than 1 hsp (galF or gnd sequences has one more hits in tested genome), that\'s unusual for for our short sequence gnd and galF, please check your submited sequence'+'\n'
     
          
  elif len(blast_record[0].alignments)>1 and len(blast_record[1].alignments)==1:
    print 'The gnd gene is splited on different contigs of your submitted sequence' +'\n'
    for record in blast_record:
      for alignment in record.alignments:
        if len(alignment.hsps)!=1:
          print '$$$ No O_type result, please check the number of hsp:  some alignment have more than 1 hsp (galF or gnd sequences has one more hits in tested genome), that\'s unusual for our short sequence gnd and galF, please check your submited sequence'+'\n'
          break
        
    print 'Each alignment has one hsps'+'\n'
     

                
    target_seq=SeqIO.parse(target,'fasta')
    target_seq=list(target_seq)
    outfile=open(output,'w')
        
    for alignment in blast_record[0].alignments:
      for contig in target_seq:
        if contig.description==alignment.hit_def and len(alignment.hsps[0].sbjct)!=len(contig.seq):
          hitstart=alignment.hsps[0].sbjct_start
          hitend=alignment.hsps[0].sbjct_end
          if hitstart>hitend:
            buffer=hitstart
            hitstart=hitend
            hitend=buffer
          potential1=contig.seq[0:hitstart-1].tostring()
          title1='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':0 to '+str(hitstart-1)+'_total:'+str(len(potential1))+'bp'
          potential2=contig.seq[hitend+1:].tostring()
          title2='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':'+str(hitend+1)+' to contig end'+'_total:'+str(len(potential2))+'bp'
          aim_sequence=potential1
          title=title1
          if len(potential1)<len(potential2):
            aim_sequence=potential2
            title=title2
          outfile.write(title)
          outfile.write('\n')
          outfile.write(aim_sequence)
          outfile.write('\n')
            
            
    for alignment in blast_record[1].alignments:
      for contig in target_seq:
        if contig.description==alignment.hit_def:
          hitstart=alignment.hsps[0].sbjct_start
          hitend=alignment.hsps[0].sbjct_end
          if hitstart>hitend:
            buffer=hitstart
            hitstart=hitend
            hitend=buffer
          potential1=contig.seq[0:hitstart-1].tostring()
          title1='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':0 to '+str(hitstart-1)+'_total:'+str(len(potential1))+'bp'
          potential2=contig.seq[hitend+1:].tostring()
          title2='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':'+str(hitend+1)+' to contig end'+'_total:'+str(len(potential2))+'bp'
          aim_sequence=potential1
          outfile.write(title1)
          outfile.write('\n')
          outfile.write(potential1)
          outfile.write('\n')
          outfile.write(title2)
          outfile.write('\n')
          outfile.write(potential2)
          outfile.write('\n')



    outfile.close()

    os.system('makeblastdb -in '+database+' -out '+database+'_db '+'-dbtype nucl')
    os.system('blastn -query '+output+' -db '+database+'_db '+'-out '+'Blast_Otype_'+target+'.xml '+'-outfmt 5')
    print '\n'
    xml_file2='Blast_Otype_'+target+'.xml'
    filehandle=open(xml_file2)
    records=NCBIXML.parse(filehandle)
    records=list(records)
    print len(records)
    realrecord1=records[0]
    if len(records[1].alignments)>len(records[0].alignments):
      realrecord1=records[1]
    realrecord2=records[2]
    if len(records[3].alignments)>len(records[2].alignments):
      realrecord2=records[3]
    title='>'+records[0].query.split(':')[0]+'_combined_potential_sequences'
    filecontent=SeqIO.parse(output,'fasta')
    filecontent=list(filecontent)
            
    for contig in filecontent:
      if contig.description==realrecord1.query:
        sequence=contig.seq.tostring()

    for contig in filecontent:
      if contig.description==realrecord2.query:
        sequence=sequence+contig.seq.tostring()
         

    outfile=open('combined_sequence.fasta','w')
    outfile.write(title)
    outfile.write('\n')
    outfile.write(sequence)
    outfile.close()
                
    os.system('blastn -query combined_sequence.fasta'+' -db '+database+'_db '+'-out '+'Combined_seq_blast_'+target+'.xml '+'-outfmt 5')
    xml_file='Combined_seq_blast_'+target+'.xml'
    print '\n'


    filehandle=open(xml_file)
    records=NCBIXML.parse(filehandle)
    records=list(records)
    E_thresh=1e-10  
    show_result()



  elif len(blast_record[0].alignments)==1 and len(blast_record[1].alignments)>1:
    print 'The galF gene is splited on different contigs of your submitted sequence' +'\n'
    for record in blast_record:
      for alignment in record.alignments:
        if len(alignment.hsps)!=1:
          print '$$$ No O_type result, please check the number of hsps:  some alignment have more than 1 hsp (galF or gnd sequences has one more hits in tested genome), that\'s unusual for our short sequence gnd and galF, please check your submited sequence'+'\n'
          break
        
    print 'Each alignment has one hsp'+'\n'                 
    outfile=open(output,'w')

    for alignment in blast_record[0].alignments:
      for contig in target_seq:
        if contig.description==alignment.hit_def:
          hitstart=alignment.hsps[0].sbjct_start
          hitend=alignment.hsps[0].sbjct_end
          if hitstart>hitend:
            buffer=hitstart
            hitstart=hitend
            hitend=buffer
          potential1=contig.seq[0:hitstart-1].tostring()
          title1='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':0 to '+str(hitstart-1)+'_total:'+str(len(potential1))+'bp'
          potential2=contig.seq[hitend+1:].tostring()
          title2='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':'+str(hitend+1)+' to contig end'+'_total:'+str(len(potential2))+'bp'
          aim_sequence=potential1
          outfile.write(title1)
          outfile.write('\n')
          outfile.write(potential1)
          outfile.write('\n')
          outfile.write(title2)
          outfile.write('\n')
          outfile.write(potential2)
          outfile.write('\n')

        
    for alignment in blast_record[1].alignments:
      for contig in target_seq:
        if contig.description==alignment.hit_def and len(alignment.hsps[0].sbjct)!=len(contig.seq):
          hitstart=alignment.hsps[0].sbjct_start
          hitend=alignment.hsps[0].sbjct_end
          if hitstart>hitend:
            buffer=hitstart
            hitstart=hitend
            hitend=buffer
          potential1=contig.seq[0:hitstart-1].tostring()
          title1='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':0 to '+str(hitstart-1)+'_total:'+str(len(potential1))+'bp'
          potential2=contig.seq[hitend+1:].tostring()
          title2='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':'+str(hitend+1)+' to contig end'+'_total:'+str(len(potential2))+'bp'
          aim_sequence=potential1
          title=title1
          if len(potential1)<len(potential2):
            aim_sequence=potential2
            title=title2
          outfile.write(title)
          outfile.write('\n')
          outfile.write(aim_sequence)
          outfile.write('\n')

         
    outfile.close()

    os.system('makeblastdb -in '+database+' -out '+database+'_db '+'-dbtype nucl')
    os.system('blastn -query '+output+' -db '+database+'_db '+'-out '+'Blast_Otype_'+target+'.xml '+'-outfmt 5')
    print '\n'
    xml_file='Blast_Otype_'+target+'.xml'
    filehandle=open(xml_file)
    records=NCBIXML.parse(filehandle)
    records=list(records)
    realrecord1=records[0]
    if len(records[1].alignments)>len(records[0].alignments):
      realrecord1=records[1]
    realrecord2=records[2]
    if len(records[3].alignments)>len(records[2].alignments):
      realrecord2=records[3]
    title='>'+records[0].query.split(':')[0]+'_combined_potential_sequences'
    filecontent=SeqIO.parse(output,'fasta')
    filecontent=list(filecontent)
            
    for contig in filecontent:
      if contig.description==realrecord1.query:
        sequence=contig.seq.tostring()

    for contig in filecontent:
      if contig.description==realrecord2.query:
        sequence=sequence+contig.seq.tostring()
         

    outfile=open('combined_sequence.fasta','w')
    outfile.write(title)
    outfile.write('\n')
    outfile.write(sequence)
    outfile.close()
                
    os.system('blastn -query combined_sequence.fasta'+' -db '+database+'_db '+'-out '+'Combined_seq_blast_'+target+'.xml '+'-outfmt 5')
    xml_file='Combined_seq_blast_'+target+'.xml'
    print '\n'


    filehandle=open(xml_file)
    records=NCBIXML.parse(filehandle)
    records=list(records)
    E_thresh=1e-10
    show_result()
     

  elif len(blast_record[0].alignments)>1 and len(blast_record[1].alignments)>1:
    print 'The gnd and galF gene are both splited on different contigs of your submitted sequence' +'\n'
    for record in blast_record:
      for alignment in record.alignments:
        if len(alignment.hsps)!=1:
          print '$$$ No O_type result, please check the number of hsp:  some alignment have more than 1 hsp (galF or gnd sequences has one more hits in tested genome), that\'s unusual for our short sequence gnd and galF, please check your submited sequence'+'\n'
          break
        
    print 'Each alignment has one hsps'+'\n'
    outfile=open(output,'w')

    for alignment in blast_record[0].alignments:
      for contig in target_seq:
        if contig.description==alignment.hit_def and len(alignment.hsps[0].sbjct)!=len(contig.seq):
          hitstart=alignment.hsps[0].sbjct_start
          hitend=alignment.hsps[0].sbjct_end
          if hitstart>hitend:
            buffer=hitstart
            hitstart=hitend
            hitend=buffer
          potential1=contig.seq[0:hitstart-1].tostring()
          title1='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':0 to '+str(hitstart-1)+'_total:'+str(len(potential1))+'bp'
          potential2=contig.seq[hitend+1:].tostring()
          title2='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':'+str(hitend+1)+' to contig end'+'_total:'+str(len(potential2))+'bp'
          aim_sequence=potential1
          title=title1
          if len(potential1)<len(potential2):
            aim_sequence=potential2
            title=title2
          outfile.write(title)
          outfile.write('\n')
          outfile.write(aim_sequence)
          outfile.write('\n')

        
    for alignment in blast_record[1].alignments:
      for contig in target_seq:
        if contig.description==alignment.hit_def and len(alignment.hsps[0].sbjct)!=len(contig.seq):
          hitstart=alignment.hsps[0].sbjct_start
          hitend=alignment.hsps[0].sbjct_end
          if hitstart>hitend:
            buffer=hitstart
            hitstart=hitend
            hitend=buffer
          potential1=contig.seq[0:hitstart-1].tostring()
          title1='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':0 to '+str(hitstart-1)+'_total:'+str(len(potential1))+'bp'
          potential2=contig.seq[hitend+1:].tostring()
          title2='>'+target.split('.')[0]+' potential rfb region:'+alignment.hit_def+':'+str(hitend+1)+' to contig end'+'_total:'+str(len(potential2))+'bp'
          aim_sequence=potential1
          title=title1
          if len(potential1)<len(potential2):
            aim_sequence=potential2
            title=title2
          outfile.write(title)
          outfile.write('\n')
          outfile.write(aim_sequence)
          outfile.write('\n')

         
    outfile.close()

    os.system('makeblastdb -in '+database+' -out '+database+'_db '+'-dbtype nucl')
    os.system('blastn -query '+output+' -db '+database+'_db '+'-out '+'Blast_Otype_'+target+'.xml '+'-outfmt 5')
    xml_file='Blast_Otype_'+target+'.xml'
    print '\n'
    filehandle=open(xml_file)
    records=NCBIXML.parse(filehandle)
    records=list(records)
    realrecord1=records[0]
    if len(records[1].alignments)>len(records[0].alignments):
      realrecord1=records[1]
    realrecord2=records[2]
    if len(records[3].alignments)>len(records[2].alignments):
      realrecord2=records[3]
    title='>'+records[0].query.split(':')[0]+'_combined_potential_sequences'
    filecontent=SeqIO.parse(output,'fasta')
    filecontent=list(filecontent)
            
    for contig in filecontent:
      if contig.description==realrecord1.query:
        sequence=contig.seq.tostring()

    for contig in filecontent:
      if contig.description==realrecord2.query:
        sequence=sequence+contig.seq.tostring()
         

    outfile=open('combined_sequence.fasta','w')
    outfile.write(title)
    outfile.write('\n')
    outfile.write(sequence)
    outfile.close()
                
    os.system('blastn -query combined_sequence.fasta'+' -db '+database+'_db '+'-out '+'Combined_seq_blast_'+target+'.xml '+'-outfmt 5')
    xml_file='Combined_seq_blast_'+target+'.xml'
    print '\n'


    filehandle=open(xml_file)
    records=NCBIXML.parse(filehandle)
    records=list(records)
    E_thresh=1e-10
    show_result()



else:
  print '$$$ $$$ No O_type result, Attention: unusual number of hits, no hits for galF or gnd! Check blast output...'+'\n'


os.system('rm '+target+'_db.'+'*')



 

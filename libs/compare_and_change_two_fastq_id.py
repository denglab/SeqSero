#!/usr/bin/env python

import os,sys
file1=sys.argv[1]
file2=sys.argv[2]


def compare_and_change_two_fastq_id(file1,file2):
  a=os.popen("head "+file1).read().split("\n")
  b=os.popen("head "+file2).read().split("\n")
  for x in a:
    if x.startswith("@"):
      a_title=x.split(" ")[0]
  for x in b:
    if x.startswith("@"):
      b_title=x.split(" ")[0]
  if a_title==b_title:
    pass
  else:
    print "changing the title of two seperated fastq files..."
    print a_title,b_title
    os.system("sed "+"-i 's/.1 / /g' "+file1)
    print "finished file1"
    os.system("sed "+"-i 's/.2 / /g' "+file2)
    print "finished file2"

compare_and_change_two_fastq_id(file1,file2) 
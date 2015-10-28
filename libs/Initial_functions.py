#!/usr/bin/env python


def To_list(L):
  import string
  New_list=[]
  for x in L:
    x1=x[:-1]
    x1=string.atoi(x1)
    New_list.append(x1)
  return New_list


def Uniq(L): #return the uniq list and the count number
  Old=L
  L.sort()
  L = [L[i] for i in range(len(L)) if L[i] not in L[:i]]
  count=[]
  for j in range(len(L)):
    y=0
    for x in Old:
      if L[j]==x:
        y+=1
    count.append(y)
  return (L,count)


def Parse_seros_in_genome_trakr(L): #return the sero names in the sra_result.xlsx, the next step is usually "Uniq" in above 
  names2=[]
  for x in L:
    if "serovar" in x:
      key_word=x.split("serovar")[1].split("_")[1]  #the seronames
      if key_word!="str." and key_word!="group": #to eliminate some "serovar_str." and "serovar_group"
        if key_word in ["I","II","III","IIIa","IIIb","IV","VI","B"]: #the serovar is behind those letters
          if x.split("serovar")[1].split("_")[2]!="str.":
            names2.append(x.split("serovar")[1].split("_")[2])
        else:
          names2.append(x.split("serovar")[1].split("_")[1])
  return names2




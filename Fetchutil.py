from Bio import Entrez
import subprocess,os
Entrez.email=os.environ['ENTREZ_EMAIL']

from socket import error as SocketError
import errno
import time

def seqfetch(acc):
  "prints the sequence in fasta format to file"
  if os.path.exists('Proteomes/orgmap_all'):
    org=orgfetch(acc)
    binom=org[0]
    bspl=binom.split()
    four=bspl[0][0]+bspl[1][:3]
    four.lower()
    ##DB params
    db='Proteomes/'+four
    out='Orthos/'+acc+'.fas'
    subprocess.call(['blastdbcmd','-db',db,'-out',out,'-entry',acc])
    new_out=[]
    with open(out) as out_fil:
      new_out=out_fil.readlines()
      new_out[0]='>'+org[0]+': '+org[1]+' '+acc+'\n'
    with open(out,'w') as new_fil:
      new_fil.write(''.join(new_out))
    return out
  else:
    hand=Entrez.efetch(db='protein',id=str(acc),rettype='fasta')
    string=''
    hand.readline()
    org=orgfetch(acc)
    string+='>'+org[0]+': '+org[1]+' '+acc+'\n'
    while hand:
      inpout=hand.readline()
      string+=inpout
      if inpout=='':
        if not os.path.exists('Orthos'):
          subprocess.call(['mkdir','Orthos'])
          open('Orthos/'+acc+'.fasta','w').write(string)
          return 'Orthos/'+acc+'.fasta'

def addseq(oldseq, newseq):
  "Transfers the sequence from newseq to oldseq"
  old=open(oldseq,'a')
  new=open(newseq)
  while new:
    add=new.readline()
    old.write(add)
    if add=='':
      return
def namefetch(acc):
  "Returns the description of acc"
  hand=Entrez.esummary(db='protein',id=acc)
  res=Entrez.read(hand)[0]
  hand.close()
  return res['Title']
#def namefetch(acc):
#    "Fetchs definition of specified gene product"
#    hand=Entrez.efetch(db="protein", id=acc ,rettype='gp')
#    record=''
#    while hand:
#        read=hand.readline()
#        record+=read.strip()+'  '
#        if read=='':
#            break
#    rec=record.split('  ')
#    recun=[]
#    for s in rec:
#        if s!='':
#            recun.append(s)
#    i=0
#    name=''
#    while i<len(recun):
#        if recun[i].endswith('DEFINITION'):
#            while i:
#                i+=1
#                if recun[i].startswith('ACCESSION'):
#                    return name
#                name+=recun[i]+' '
#        i+=1

def orgfetch(acc):
  "Fetches source organism for acc"
  ret=[]
  if os.path.exists('Proteomes/orgmap_all'):
    dict_all={}
    with open('Proteomes/orgmap_all') as infile:
      for line in infile:
        spl=line.split()
        dict_all.update({spl[0]:spl[1:]})
    taxid=dict_all[acc][0]
    taxall={}
    with open('Proteomes/all_tax') as infile:
      for line in infile:
        spl=line.split('\t')
        taxall.update({spl[0]:spl[1]})
      binom=taxall[taxid]
      ret.append(binom)
      bspl=binom.split()
      four=bspl[0][0]+bspl[1][:3]
      four.lower()
      with open(four+'_tax_fetch') as taxf:
        read=Entrez.read(taxf)
        ret.append(read[0]['LineageEx'][1]['ScientificName'])
        return ret
  else:
    try:
      hand=Entrez.esummary(db="protein", id=acc)
    except SocketError:
#    if e.errno != errno.ECONNRESET:
#      raise # Not error we are looking for
      time.sleep(1.5)
      hand=Entrez.esummary(db="protein", id=acc)
    taxid=Entrez.read(hand)[0]['TaxId']
    hand.close()
    try:
      hand=Entrez.efetch(db='taxonomy',id=str(taxid))
    except SocketError:
#    if e.errno != errno.ECONNRESET:
#      raise # Not error we are looking for
      time.sleep(1.5)
      hand=Entrez.efetch(db='taxonomy',id=str(taxid))
    result=Entrez.read(hand)[0]
    hand.close()
    ret.append(result['ScientificName'])
    for i in result['LineageEx']:
      if i['Rank']=='superkingdom':
        ret.append(i['ScientificName'])
        return ret
    return ret
    
#def orgfetch(acc):
#    "Fetchs the organism from which acc accession number of a protein came from"
#    hand=Entrez.efetch(db="protein", id=acc ,rettype='gp')
#    record=''
#    while hand:
#        read=hand.readline()
#        record+=read.strip()+'  '
#        if read=='':
#            break
#    rec=record.split('  ')
#    recun=[]
#    for s in rec:
#        if s!='':
#            recun.append(s)
#    i=0
#    ret=[]
#    retu=['']
#    #print recun
#    mutisp=False
#    while i<len(recun):
#        j=0
#        if recun[i]=='DEFINITION':
#            multisp=recun[i+1].startswith("MULTISPECIES:")
#            if not multisp:
#              tab=recun[i+1].strip().split()
#              if not recun[i+2]=='ACCESSION':
#                  tab+=recun[i+2].strip().split()
#              for k in range(len(tab)):
#                  if tab[k].startswith('[') and len(tab[k])>5:
#                      retu=[tab[k][1:]]
#                      while k<len(tab):
#                          k+=1
#                          if tab[k].endswith('].'):
#                              retu[0]+=' '+tab[k][0:tab[k].index(']')]
#                              retu.append('Unknown.')
#                              break
#                          else:
#                              retu[0]+=' '+tab[k]
#            else:
#              return [recun[i+11]+' multispecies',recun[i+14].split(';')[0]]
##        elif recun[i+1].startswith("MULTISPECIES:"):
##          ret=
#        if recun[i]=='ORGANISM':
#            try:
#                ret=[recun[i+1].split()[0]+' '+recun[i+1].split()[1]]
#            except:
#                return retu
#            if ret[0].endswith('.'):
#                ret[0]+=' '+recun[i+1].split()[2]
#            ret.append(recun[i+2].split(';')[0])
#            return ret
#
#        i+=1

def toGI(acc):
  "Returns GI number for acc"
  hand=Entrez.esearch(db='protein', term=acc)
  res=Entrez.read(hand)
  return res['IdList'][0]

#def toGI(acc):
#    "Converts given accession no acc to a GI number"
#    hand=Entrez.efetch(db="protein", id=acc ,rettype='gp')
#    record=''
#    while hand:
#        read=hand.readline()
#        record+=read.strip()+'  '
#        if read=='':
#            break
#    rec=record.split('  ')
#    recun=[]
#    for s in rec:
#        if s!='':
#            recun.append(s)
#    string=''
#    for g in recun:
#        if g.startswith('GI:'):
#            for char in g[3:]:
#                if char.isdigit():
#                    string+=char
#            return string

def findlonglen(dicto):
    "finds the length of the longest key in dicto"
    keys=dicto.keys()
    pivot=''
    for k in keys:
        if len(k)>len(pivot):
            pivot=k
    return len(pivot)

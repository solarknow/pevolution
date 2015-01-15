from Bio import Entrez
import subprocess,os
Entrez.email=os.environ['ENTREZ_EMAIL']

def seqfetch(acc,local):
  "prints the sequence in fasta format to file"
  if local:
    org=orgfetch(acc,local)
    binom=org[0]
    bspl=binom.split()
    four=bspl[0][0]+bspl[1][:3]
    four=four.lower()
    ##DB params
    db='Proteomes/'+four
    out='Orthos/'+acc+'.fasta'
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
    lines=iter(hand.readlines())
    hand.close()
    if not os.path.exists('dicts'):
      os.mkdir('dicts')
   #print len(lines)
    org=orgfetch(acc,local)
    string+='>'+org[0]+': '+org[1]+' '+acc+'\n'
    for l in lines:
      inpout=l
      if inpout.startswith('>'):
        continue
      string+=inpout
    with open('Orthos/'+acc+'.fasta','w') as op:
      op.write(string)
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

def orgfetch(acc,local):
  "Fetches source organism for acc"
  ret=[]
  if local:
    dict_all={}
    with open('Proteomes/orgmap_all') as infile:
      for line in infile:
        spl=line.split()
        #print spl
        dict_all.update({spl[0]:spl[1]})
    taxid=dict_all[acc]
    taxall={}
    with open('Proteomes/all_tax') as infile:
      for line in infile:
        spl=line.split('\t')
        taxall.update({spl[0]:spl[1]})
      binom=taxall[taxid].strip()
      ret.append(binom)
      bspl=binom.split()
      four=bspl[0][0]+bspl[1][:3]
      four=four.lower()
      with open('Proteomes/'+four+'_tax_fetch') as taxf:
        read=Entrez.read(taxf)
        ret.append(read[0]['LineageEx'][1]['ScientificName'])
        return ret
  else:
    hand=Entrez.esummary(db="protein", id=acc)
    taxid=Entrez.read(hand)[0]['TaxId']
    hand.close()
    hand=Entrez.efetch(db='taxonomy',id=str(taxid))
    result=Entrez.read(hand)[0]
    hand.close()
    ret.append(result['ScientificName'])
    for i in result['LineageEx']:
      if i['Rank']=='superkingdom':
        ret.append(i['ScientificName'])
        return ret
    return ret
    
def toGI(acc):
  "Returns GI number for acc"
  hand=Entrez.esearch(db='protein', term=acc)
  res=Entrez.read(hand)
  return res['IdList'][0]

def findlonglen(dicto):
    "finds the length of the longest key in dicto"
    keys=dicto.keys()
    pivot=''
    for k in keys:
        if len(k)>len(pivot):
            pivot=k
    return len(pivot)

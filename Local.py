import sys, os, subprocess
from Bio import SeqIO
import Fetchutil

def generatescript(org):
  "generates a fetch script for all protein seqs of org organism from NCBI"
  DOL='$$'
  AMP='&&&'
  bi=org.split()
  four=bi[0][0]+bi[1][:3]
  four=four.lower()
  template=open('fetchtemp').readlines()
  if not os.path.exists('Proteomes'):
    subprocess.call(['mkdir','Proteomes'])
  genscript=open('Proteomes/fetch'+four,'w')

  for lin in template:
    spldol=lin.split(DOL)
    splamp=lin.split(AMP)
    if len(spldol)==3:
      genscript.write(spldol[0]+org+spldol[2])
      continue
    elif len(splamp)==3:
      genscript.write(splamp[0]+'Proteomes/'+four+splamp[2])
      continue
    else:
      genscript.write(lin)
  genscript.close()
  return four

def fetchfasta(org):
  "runs a fetch script of 4-letter org organism"
  subprocess.call(['perl', 'Proteomes/fetch'+org])
  return 'Proteomes/'+org+'.aa'

def taxidmap(list_binom):
  "Prints the taxid for binom organisms to file"
  org_map=open('orgmap','a')
  if type(list_binom) is str:
    list_binom=[list_binom]
  for org in list_binom:
    hand=Entrez.esearch(db='taxonomy', term=org+'[ORGN]')
    results=Entrez.read(hand)
    #print results
    hand.close()
    taxid=results['IdList']
    taxa=''
    for t in taxid:
      taxa+=str(t)+', '
    hand=Entrez.esummary(db='taxonomy',id=taxa)
    res=Entrez.read(hand)
    for i in range(len(res)):
      resi=res[i]
      org_map.write(resi['ScientificName']+'\t'+taxid[i]+'\n')
  org_map.close()

def seqmap(aa):
  "Prints a seq-taxid map to file for aa file"
  tax_map=open('orgmap').readlines()
  taxa_map={}
  print 'reading the orgmap'
  for i in tax_map:
    i=i.strip()
    spl=i.split('\t')
    taxa_map.update({spl[0]:spl[1]}) 
  seq_map=open(aa.split('.')[0]+'map','w')  
  hand=open(aa,'r')
  print 'parsing aa'
  for rec in SeqIO.parse(hand, 'fasta'):
    rid=rec.id.split('|')[1]
    #print rid
    org=Fetchutil.orgfetch(rid)
    try:
      seq_map.write(rid+'\t'+taxa_map[org[0]]+'\n')
    except KeyError:
      seq_map.write(rid+'\t'+org[0]+'\n')
  hand.close()
  seq_map.close()

if __name__=="__main__":
  org=sys.argv[1]
  fasta=fetchfasta(generatescript(org))
  seqmap(fasta)

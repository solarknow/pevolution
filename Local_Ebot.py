import sys, os, subprocess
from Bio import SeqIO, Entrez
import Fetchutil
####
## This version of Local makes use of NCBI's eBot and the scripts it generates.
## Thus it is not a pure pythonic implementation. You have been warned.
####

def generatescripts(org):
  "generates a fetch script for all protein seqs of org organism from NCBI"
  TEMP_DIR='templates'
  DOL='$$'
  AMP='&&&'
  AT='@@@@'
  bi=org.split()
  four=bi[0][0]+bi[1][:3]
  four=four.lower()
  for fil in os.listdir(TEMP_DIR):
    print fil
    template=open(TEMP_DIR+'/'+fil).readlines()
    if not os.path.exists('Proteomes'):
      subprocess.call(['mkdir','Proteomes'])
    genscript=open('Proteomes/'+four+'_'+fil.split('_')[1],'w')

    for lin in template:
      splat=lin.split(AT)
      spldol=lin.split(DOL)
      splamp=lin.split(AMP)
      if len(spldol)==3:
        #print spldol[0]+'+'.join(org.split())+spldol[2]
        genscript.write(spldol[0]+'+'.join(org.split())+spldol[2])
        continue
      elif len(splamp)==3:
        #print splamp[0]+'Proteomes/'+four+splamp[2]
        genscript.write(splamp[0]+'Proteomes/'+four+splamp[2])
        continue
      elif len(splat)==3:
       #print splat[0]+os.environ['ENTREZ_EMAIL']+splat[2]
        genscript.write(splat[0]+os.environ['ENTREZ_EMAIL']+splat[2])
        continue
      else:
        genscript.write(lin)
    genscript.close()
  return four

def fetchfasta(org):
  "runs a fetch script of 4-letter org organism"
  subprocess.call(['perl', 'Proteomes/'+org+'_fetch.pl'])
  return 'Proteomes/'+org+'.aa'

def taxidmap(org):
  "Prints the taxid for binom organisms to file"
  org_map=open('orgmap','a')
  subprocess.call(['perl', 'Proteomes/'+org+'_tax.pl'])
  results=Entrez.read(open('Proteomes/'+org+'_tax'))
  for i in results:
    taxid=i['TaxId']
    subprocess.call(['perl', 'Proteomes/'+org+'_summ.pl'])
    res=Entrez.read(open('Proteomes/'+org+'_summary'))
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
  print "Generating scripts"
  four=generatescripts(org)
  print "Scripts are generates"
  if not os.path.exists('Proteomes/'+four+'.aa'):
    print "Fetching seqs"
    fasta=fetchfasta(four)
  else:
    fasta='Proteomes/'+four+'.aa'
  print "Making seqmap"
  taxidmap(four)
  seqmap(fasta)

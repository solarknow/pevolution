import sys, os, subprocess, shutil, glob
from Bio import SeqIO, Entrez
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
  org_map=open('Proteomes/orgmap_'+four,'w')
  subprocess.call(['perl', 'Proteomes/'+org+'_tax.pl'])
  results=Entrez.read(open('Proteomes/'+org+'_tax'))
  #print results
  orgs={}
  for o in results:
    orgs.update({o['TaxId']:o['ScientificName']})
  ids=''
  with open('Proteomes/all_tax','a') as orgout:
    for ori in orgs:
      orgout.write(repr(ori)+'\t'+orgs[ori]+'\n')
      ids+=repr(ori)+','
  ids=ids[:-1]
  with Entrez.efetch(db='taxonomy',id=ids) as hand:
    for line in hand:
      with open(four+'_tax_fetch','w') as taxf:
        taxf.write(line)
  subprocess.call(['perl', 'Proteomes/'+org+'_summ.pl'])
  try:
    res=Entrez.read(open('Proteomes/'+org+'_summary'))
  except Entrez.Parser.CorruptedXMLError:
    fil=open('Proteomes/'+org+'_summary')
    fil_read=fil.read().split('<?xml version="1.0" encoding="UTF-8"?>\n')
    #print len(fil_read)
      #first=['<?xml version="1.0" encoding="UTF-8"?>'].append(fil_read[0].split('\n')[:-1])
    first=fil_read[1].split('\n')[:-2]
    first.insert(0,'<?xml version="1.0" encoding="UTF-8"?>')
    second=fil_read[2].split('\n')[2:]
    first.extend(second)
    fil.close()
    fil=open('Proteomes/'+org+'_summary','w')
    #print len(first)
    fil.write('\n'.join(first))
    fil.close()
    res=Entrez.read(open('Proteomes/'+org+'_summary'))
      
  for i in range(len(res)):
    resi=res[i]
    #print resi
    org_map.write(repr(resi['Gi'])+'\t'+repr(resi['TaxId'])+'\n')
  org_map.close()

#def seqmap(aa):
#  "Prints a seq-taxid map to file for aa file"
#  tax_map=open('orgmap').readlines()
#  taxa_map={}
#  print 'reading the orgmap'
#  for i in tax_map:
#    i=i.strip()
#    spl=i.split('\t')
#    taxa_map.update({spl[0]:spl[1]}) 
#  seq_map=open(aa.split('.')[0]+'map','w')  
#  hand=open(aa,'r')
#  print 'parsing aa'
#  for rec in SeqIO.parse(hand, 'fasta'):
#    rid=rec.id.split('|')[1]
#    #print rid
#    org=Fetchutil.orgfetch(rid)
#    try:
#      seq_map.write(rid+'\t'+taxa_map[org[0]]+'\n')
#    except KeyError:
#      seq_map.write(rid+'\t'+org[0]+'\n')
#  hand.close()
#  seq_map.close()

def makedb(binom,four):
  "runs makeblastdb with appropriate params"
  infile='Proteomes/'+four+'.aa'
  intype='fasta'
  dbtype='prot'
  title='Protein Database for all variants of '+binom
  parse='-parse_seqids'
  out=four
  taxmap='Proteomes/orgmap_'+four
  log='Proteomes/'+four+'.log'
  subprocess.call(['makeblastdb','-in',infile,'-input_type',intype,'-title',
                   title,parse,'-out',out,'-taxid_map',taxmap,'-logfile',log,
                   '-dbtype',dbtype])


if __name__=="__main__":
  orgs=sys.argv[1:]
  fours=[]
  taxa=open('Proteomes/taxa','w')
  for org in orgs:
    print "Fetching localdb for "+org
    print "Generating scripts"
    four=generatescripts(org)
    fours.append(four)
  #  taxa.write(four+'\n')
    print "Scripts are generated"
    if not os.path.exists('Proteomes/'+four+'.aa'):
      print "Fetching seqs"
      fasta=fetchfasta(four)
    else:
      fasta='Proteomes/'+four+'.aa'
    print "Making taxmap"
    taxidmap(four)
    if not os.path.exists('Proteomes/'+four+'.pin'):
      print 'making blastDB'
      makedb(org,four)
      print 'moving to database folder'
      dbs=glob.glob(four+'.p*')
      for g in dbs:
        shutil.move(g,'Proteomes/')
    else:
      print "DB "+four+" already made"  
 # taxa.close()
  print "merging orgmaps"
  with open('Proteomes/orgmap_all', 'w') as outfile:
    for f in fours:
      with open('Proteomes/orgmap_'+f) as infile:
        for line in infile:
          outfile.write(line)

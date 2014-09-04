import os
from Bio import Entrez
Entrez.email=os.environ['ENTREZ_EMAIL']
def post(binom):
  "Esearch for binom in taxology. Returns WebEnv and QueryKey"
  hand=Entrez.esearch(db='taxonomy', term=binom+'[ORGN]', usehistory="y")
  read=Entrez.read(hand)
  hand.close()
  spl=binom.split()
  four = spl[0][0]+spl[1][:3]
  four=four.lower()
  return (read["WebEnv"],read["QueryKey"],four)

def fetchtaxa(webenv,qkey,four):
  "Writes taxid, superkingdom and scientific names from specified post"
  hand=Entrez.efetch(db='taxonomy', webenv=webenv, query_key=qkey)
  read=Entrez.read(hand)
  hand.close()
  map=open('Proteomes/'+four+'.orgmap','w')
  for i in read:
    map.write(str(i['TaxId'])+'\t'+i['ScientificName']+'\t')
    for j in i['LineageEx']:
      if j['Rank']=='superkingdom':
        map.write(j['ScientificName']+'\n')
        break
  map.close()

def link(webenv,qkey):
  "links the taxonomy search designated, to the protein db. Returns the webenv and query_key of the new search"
  fromdb='taxonomy'
  todb='protein'
  hand=Entrez.elink(dbfrom=fromdb, db=todb, webenv=webenv, query_key=qkey, 
                    usehistory='y')
  read=Entrez.parse(hand)
#  print hand.readlines()
  hand.close()
  return (read['WebEnv'],read['QueryKey'])

def download(webenv,qkey,four):
  "Prints linked sequences to file: returns a handle to the esummary"
  batch_size=500
  out_hand=open(four+'.aa','w')
  for start in 
  hand=Entrez.efetch(db='protein', rettype='fasta', webenv=webenv,query_key=qkey)
  return Entrez.esummary(db='protein',webenv=webenv, query_key=qkey)

if __name__=="__main__":
  org=sys.argv[1]
  (we1,qk1,four)=post(org)
  open(org+'_keys','w').write(we1+'\t'+qk1+'\t'+four+'\n')
  fetchtaxa(we1,qk1,four)
  (we2,qk2)=link(we1,qk1)

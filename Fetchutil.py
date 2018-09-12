from Bio import Entrez
import subprocess,os

def set_email(email):
    Entrez.email=email

def fetch_protein(acc, fmt='fasta'):
    "returns a file-like handle of a Entrez search for acc accession number and returns in fmt format"
    fetched = Entrez.efetch(db='protein',id=str(acc),rettype=fmt)
    return fetched.readlines()
def fetch_fasta(acc):
    "prints the sequence in fasta format to file"
    hand = fetch_protein(acc)
    string=''
    org = fetch_organism(acc)
    string+='>'+org[0]+': '+org[1]+' '+acc+'\n'
    for line in hand[1:]:
        string+=line
        if line == '\n':
            if not os.path.exists('Orthos'):
                subprocess.call(['mkdir','Orthos'])
            with open('Orthos/'+acc+'.fasta','w') as writ_file:
                writ_file.write(string)
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

def fetch_definition(acc):
    "Fetchs definition of specified gene product"
    hand=fetch_protein(acc,'gp')
    record=''
    for line in hand:
        read=line
        record+=read.strip()+'  '
    rec=record.split('  ')
    recun=[]
    for s in rec:
        if s!='':
            recun.append(s)
    i=0
    name=''
    while i<len(recun):
        if recun[i].endswith('DEFINITION'):
            while i:
                i+=1
                if recun[i].startswith('ACCESSION'):
                    return name.strip()
                name+=recun[i]+' '
        i+=1


def fetch_organism(acc):
    "Fetchs the organism from which acc accession number of a protein came from"
    hand=fetch_protein(acc,'gp')
    record=''
    for line in hand:
        read=line
        record+=read.strip()+'  '
    rec=record.split('  ')
    recun=[]
    for s in rec:
        if s!='':
            recun.append(s)
    i=0
    ret=[]
    retu=['']
    multisp=False
    while i<len(recun):
        if recun[i]=='DEFINITION':
            multisp=recun[i+1].startswith("MULTISPECIES:")
            if not multisp:
              tab=recun[i+1].strip().split()
              if not recun[i+2]=='ACCESSION':
                  tab+=recun[i+2].strip().split()
              for k in range(len(tab)):
                  if tab[k].startswith('[') and len(tab[k])>5:
                      retu=[tab[k][1:]]
                      while k<len(tab):
                          k+=1
                          if tab[k].endswith('].'):
                              retu[0]+=' '+tab[k][0:tab[k].index(']')]
                              retu.append('Unknown.')
                              break
                          else:
                              retu[0]+=' '+tab[k]
            else:
              return [recun[i+11]+' multispecies',recun[i+14].split(';')[0]]

        if recun[i]=='ORGANISM':
            try:
                ret=[recun[i+1].split()[0]+' '+recun[i+1].split()[1]]
            except:
                return retu
            if ret[0].endswith('.'):
                ret[0]+=' '+recun[i+1].split()[2]
            ret.append(recun[i+2].split(';')[0])
            return ret

        i+=1

def findlonglen(dicto):
    "finds the length of the longest key in dicto"
    keys=dicto.keys()
    pivot=''
    for k in keys:
        if len(k)>len(pivot):
            pivot=k
    return len(pivot)

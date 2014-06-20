from Bio import Entrez
Entrez.email='mihir.sarwade@effem.com'

def seqfetch(acc):
    "prints the sequence in fasta format to file"
    hand=Entrez.efetch(db='protein',id=str(acc),rettype='fasta')
    string=''
    hand.readline()
    org=orgfetch(acc)
    string+='>'+org[0]+': '+org[1]+' '+acc+'\n'
    while hand:
        inpout=hand.readline()
        string+=inpout
        if inpout=='':
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
    "Fetchs definition of specified gene product"
    hand=Entrez.efetch(db="protein", id=acc ,rettype='gp')
    record=''
    while hand:
        read=hand.readline()
        record+=read.strip()+'  '
        if read=='':
            break
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
                    return name
                name+=recun[i]+' '
        i+=1


def orgfetch(acc):
    "Fetchs the organism from which acc accession number of a protein came from"
    hand=Entrez.efetch(db="protein", id=acc ,rettype='gp')
    record=''
    while hand:
        read=hand.readline()
        record+=read.strip()+'  '
        if read=='':
            break
    rec=record.split('  ')
    recun=[]
    for s in rec:
        if s!='':
            recun.append(s)
    i=0
    ret=[]
    retu=['']
    while i<len(recun):
        j=0
        if recun[i]=='DEFINITION':
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

def toGI(acc):
    "Converts given accession no acc to a GI number"
    hand=Entrez.efetch(db="protein", id=acc ,rettype='gp')
    record=''
    while hand:
        read=hand.readline()
        record+=read.strip()+'  '
        if read=='':
            break
    rec=record.split('  ')
    recun=[]
    for s in rec:
        if s!='':
            recun.append(s)
    string=''
    for g in recun:
        if g.startswith('GI:'):
            for char in g[3:]:
                if char.isdigit():
                    string+=char
            return string

def findlonglen(dicto):
    "finds the length of the longest key in dicto"
    keys=dicto.keys()
    pivot=''
    for k in keys:
        if len(k)>len(pivot):
            pivot=k
    return len(pivot)

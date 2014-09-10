from Bio.Blast import NCBIXML
import Fetchutil, random
import os, subprocess
from multiprocessing import Queue
#import time

if not os.path.exists('XML'):
  os.mkdir('XML')
if not os.path.exists('dicts'):
  os.mkdir('dicts')

def bestrecipblast(org, seed, thresh, queue):
    "Returns the best pairwise reciprocal BLAST using seed accession no. from seedorg organism against orgs list of organisms"
    #start=time.time()
    seedorg=Fetchutil.orgfetch(seed)
    binspl=seedorg[0].split()
    four=binspl[0][0]+binspl[1][:3]
    four.lower()
    orgspl=org.split()
    orgfour=orgspl[0][0]+orgspl[1][:3]
    orgfour.lower()
    acclist={}
    seed=Fetchutil.toGI(seed)
    #for i in orgs:
    print org
    ac=[]
    Fetchutil.seqfetch(seed)
    dum=str(int(int(seed)*random.random()))
        
    subprocess.call(['blastp', '-db', 'Proteomes/'+orgfour, '-query', 'Orthos/'+seed+'.fasta', '-evalue',str(thresh),
              '-out', 'XML/'+dum+'.xml', '-outfmt', '5', '-use_sw_tback'])
    qoutput=open('XML/'+dum+'.xml')
        
    parser=NCBIXML.parse(qoutput)
    for lin in parser:
        for align in lin.alignments:
            for hsp in align.hsps:
                if (hsp.positives/float(hsp.align_length))>=.4 and (float(hsp.align_length)/len(hsp.query))>=.25:
                    ac.append(align.title.split('|')[1])
    print "Done. Number of sequences found: "+repr(len(ac))

    for o in ac:
        print o
        Fetchutil.seqfetch(o)
        os.system('blastp -db Proteomes/'+four+' -query Orthos/'+o+'.fasta -evalue '+str(thresh)+
              ' -out XML/'+dum+'.xml -outfmt 5 -use_sw_tback')
        q1output=open('XML/'+dum+'.xml')
        parse=NCBIXML.parse(q1output)
        acc=[]
        print 'blasted'
        for lin in parse:
            for align in lin.alignments:
                for hsp in align.hsps:
                    if (hsp.positives/float(hsp.align_length))>=.4 and (float(hsp.align_length)/len(hsp.query))>.25:
                        acc.append(align.title.split('|')[1])
                    else:
                        continue

        print "Done. Number of sequences found: "+repr(len(acc))
            
        if seed in acc:
            print 'it\'s twue!'
            name=Fetchutil.orgfetch(o)[0]
            try:
                acclist[name]=[o,str(ac.index(o)+1)+'/'+str(len(ac)),str(acc.index(seed)+1)+'/'+str(len(acc))]
            except KeyError:
                acclist.update({name:[o,str(ac.index(o)+1)+'/'+str(len(ac)),str(acc.index(seed)+1)+'/'+str(len(acc))]})
                
            open('dicts/'+seed,'a').write(str(acclist)+'\n')
            break
	#elapsed=time.time()-start
	#print "Time elapsed: "+time.strftime('%M:%S',[elapsed])
    if queue is not None:
      queue.put(acclist)
    else:
      return acclist

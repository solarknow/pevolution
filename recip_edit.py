from Bio.Blast import NCBIXML
import Fetchutil, random
import os

def bestrecipblast(orgs, seed, thresh):
    "Returns the best pairwise reciprocal BLAST using seed accession no. from seedorg organism against orgs list of organisms"
    seedorg=Fetchutil.orgfetch(seed)
    acclist={}
    seed=Fetchutil.toGI(seed)
    for i in orgs:
        print i
        ac=[]
        Fetchutil.seqfetch(seed)
        dum=str(int(int(seed)*random.random()))
        os.system('blastp -db nr -query Orthos/'+seed+'.fasta -evalue '+str(thresh)+
                  ' -out XML/'+dum+'.xml -outfmt 5 -entrez_query \"'+i+'[ORGN]\" -use_sw_tback'+
                  ' -soft_masking True -remote')
        qoutput=open('XML/'+dum+'.xml')
        parser=NCBIXML.parse(qoutput)
        for lin in parser:
            for align in lin.alignments:
                for hsp in align.hsps:
                    if (hsp.positives/float(hsp.align_length))>=.3 and (float(hsp.align_length)/len(hsp.query))>=.25:
                        ac.append(align.title.split('|')[1])
        print ac

        for o in ac:
            print o
            Fetchutil.seqfetch(o)
            os.system('blastp -db nr -query Orthos/'+o+'.fasta -evalue '+str(thresh)+
                  ' -out XML/'+dum+'.xml -outfmt 5 -entrez_query \"'+seedorg[0]+'[ORGN]\" -use_sw_tback'+
                  ' -soft_masking True -remote')
            q1output=open('XML/'+dum+'.xml')
            parse=NCBIXML.parse(q1output)
            acc=[]
            print 'blasted'
            for lin in parse:
                for align in lin.alignments:
                    for hsp in align.hsps:
                        if (hsp.positives/float(hsp.align_length))>=.3 and (float(hsp.align_length)/len(hsp.query))>.25:
                            acc.append(align.title.split('|')[1])
                        else:
                            continue

            print acc

            if seed in acc:
                print 'it\'s twue!'
                name=Fetchutil.orgfetch(o)[0]
                try:
                    acclist[name]=[o,str(ac.index(o)+1)+'/'+str(len(ac)),str(acc.index(seed)+1)+'/'+str(len(acc))]
                except KeyError:
                    acclist.update({name:[o,str(ac.index(o)+1)+'/'+str(len(ac)),str(acc.index(seed)+1)+'/'+str(len(acc))]})
                open('dicts/'+seed,'a').write(str(acclist)+'\n')
                break
    return acclist

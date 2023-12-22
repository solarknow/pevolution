import os

from Bio.Blast import NCBIXML

import FetchUtil
import random

if not os.path.exists('XML'):
    os.mkdir('XML')
if not os.path.exists('dicts'):
    os.mkdir('dicts')


def best_reciprocal_blast(org, seed, thresh=5, queue=None):
    """Returns the best pairwise reciprocal BLAST using seed accession no. from
    against org organism"""
    seedorg = FetchUtil.fetch_organism(seed)[0]
    acclist = {}
    ac = []
    FetchUtil.fetch_fasta(seed)
    dum = str(int(int(seed.split('.')[0][-5:]) * random.random()))

    os.system('blastp -db nr -query Orthos/' + seed + '.fasta -evalue ' + str(thresh) +
              ' -out XML/' + dum + '.xml -outfmt 5 -entrez_query \"' + org + "[ORGN]\"" +
              ' -remote')
    with open('XML/' + dum + '.xml') as qoutput:
        parser = NCBIXML.parse(qoutput)
        for lin in parser:
            for align in lin.alignments:
                for hsp in align.hsps:
                    if (hsp.positives / float(hsp.align_length)) >= .4 and (
                            float(hsp.align_length) / len(hsp.query)) >= .25:
                        ac.append(align.title.split('|')[1])
        print("Done. Number of sequences found: " + repr(len(ac)))

    for o in ac:
        print(o)
        FetchUtil.fetch_fasta(o)
        os.system('blastp -db nr -query Orthos/' + o + '.fasta -evalue ' + str(thresh) +
                  ' -out XML/' + dum + '.xml -outfmt 5 -entrez_query \"' + seedorg[0] + '[ORGN]\"' +
                  ' -remote')
        with open('XML/' + dum + '.xml') as q1output:
            parse = NCBIXML.parse(q1output)
            acc = []
            print('blasted')
            for lin in parse:
                for align in lin.alignments:
                    for hsp in align.hsps:
                        if (hsp.positives / float(hsp.align_length)) >= .4 and (
                                float(hsp.align_length) / len(hsp.query)) > .25:
                            acc.append(align.title.split('|')[1])
                        else:
                            continue

            print("Done. Number of sequences found: " + repr(len(acc)))

        if seed in acc:
            print('it exists!')
            name = FetchUtil.fetch_organism(o)[0]
            try:
                acclist[name] = [o, str(ac.index(o) + 1) + '/' + str(len(ac)),
                                 str(acc.index(seed) + 1) + '/' + str(len(acc))]
            except KeyError:
                acclist.update({name: [o, str(ac.index(o) + 1) + '/' + str(len(ac)),
                                       str(acc.index(seed) + 1) + '/' + str(len(acc))]})

            open('dicts/' + seed, 'a').write(str(acclist) + '\n')
            break
    if queue is not None:
        queue.put(acclist)
    else:
        return acclist

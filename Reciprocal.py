import random

from Bio.Blast import NCBIXML

import FetchUtil
from helpers import commands
from helpers.constants import XMLPath, DictsPath


def best_reciprocal_blast(org, seed, thresh=5):
    """Returns the best pairwise reciprocal BLAST using seed accession no. from against org organism"""
    seedorg = FetchUtil.fetch_organism(seed)[0]
    acclist = {}
    ac = []
    FetchUtil.fetch_fasta(seed)
    dum = str(int(int(seed.split('.')[0][-5:]) * random.random()))

    commands.run_blast(seed, thresh, dum, org)
    with open(str(XMLPath(dum + '.xml'))) as qoutput:
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
        commands.run_blast(o, thresh, dum, seedorg[0])
        with open(str(XMLPath(dum + '.xml'))) as q1output:
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
            print("it's twue!")
            name = FetchUtil.fetch_organism(o)[0]
            acclist[name] = [o, str(ac.index(o) + 1) + '/' + str(len(ac)),
                             str(acc.index(seed) + 1) + '/' + str(len(acc))]
            with open(str(DictsPath(seed)), 'a') as dicts:
                dicts.write(str(acclist) + '\n')
            break
    else:
        return acclist

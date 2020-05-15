import os

from Bio.Blast import NCBIXML

import FetchUtil
import random

os.makedirs('XML', exist_ok=True)
os.makedirs('dicts', exist_ok=True)


def bestrecipblast(org, source, thresh=5, queue=None):
    """Returns the best pairwise reciprocal BLAST using source accession no. from
    against org organism"""

    sourceorg = FetchUtil.fetch_organism(source)[0]
    acclist = {}
    ac = []
    print("Source Organism: " + str(sourceorg))
    FetchUtil.fetch_fasta(source)
    file_name = str(int(int(source.split('.')[0][-5:]) * random.random()))
    print("Using {} as filename".format(file_name))
    print('blastp -db nr -query Orthos' + os.sep + source + '.fasta -evalue ' + str(thresh) +
              ' -out XML' + os.sep + file_name + '.xml -outfmt 5 -entrez_query \"' + org + '[ORGN]\"' +
              ' -remote')
    os.system('blastp -db nr -query Orthos' + os.sep + source + '.fasta -evalue ' + str(thresh) +
              ' -out XML' + os.sep + file_name + '.xml -outfmt 5 -entrez_query \"' + org + '[ORGN]\"' +
              ' -remote')
    with open('XML' + os.sep + file_name + '.xml') as qoutput:
        parser = NCBIXML.parse(qoutput)
        for lin in parser:
            for align in lin.alignments:
                for hsp in align.hsps:
                    if (hsp.positives / float(hsp.align_length)) >= .4 and (
                            float(hsp.align_length) / len(hsp.query)) >= .25:
                        ac.append(align.title.split('|')[1])
    print("First BLAST Done. Number of sequences found: " + repr(ac))

    for o in ac:
        print("BLASTING back to " + sourceorg[0])
        print(o)
        acc = []
        FetchUtil.fetch_fasta(o)
        os.system('blastp -db nr -query Orthos' + os.sep + o + '.fasta -evalue ' + str(thresh) +
                  ' -out XML' + os.sep + file_name + '.xml -outfmt 5 -entrez_query \"' + sourceorg[0] + '[ORGN]\" -use_sw_tback' +
                  ' -remote')
        with open('XML' + os.sep + file_name + '.xml') as q1output:
            parse = NCBIXML.parse(q1output)
            print('blasted')
            for lin in parse:
                for align in lin.alignments:
                    for hsp in align.hsps:
                        if (hsp.positives / float(hsp.align_length)) >= .4 and (
                                float(hsp.align_length) / len(hsp.query)) > .25:
                            acc.append(align.title.split('|')[1])
        print("Done. Number of sequences found: " + repr(len(acc)))

        if source in acc:
            print('it\'s twue!')
            name = FetchUtil.fetch_organism(o)[0]
            acclist[name] = [
                o,
                str(ac.index(o) + 1) + '/' + str(len(ac)),
                str(acc.index(source) + 1) + '/' + str(len(acc))
            ]
            with open('dicts' + os.sep + source, 'a') as seed_file:
                seed_file.write(str(acclist) + '\n')
            break
    if queue is not None:
        queue.put(acclist)
    else:
        return acclist

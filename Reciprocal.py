import os

from Bio.Blast import NCBIXML

import FetchUtil
from constants import ORTHOS_PATH, XML_PATH, DICTS_PATH

os.makedirs(XML_PATH, exist_ok=True)
os.makedirs(DICTS_PATH, exist_ok=True)


def bestrecipblast(org, source, thresh=5):
    """
    Returns the best pairwise reciprocal BLAST using source accession no. from
    against org organism
    :param org: Target binomial organism
    :param source: source accession Id
    :param thresh: E-value threshold
    :return: Dict[Binomial Organism: List[Best Accession Number, place in target search, place in source search]]
    """
    source_protein = FetchUtil.fetch_protein(source)
    source_binomial = source_protein.binomial
    acclist = {}
    ac = []
    print("Source Organism: " + str(source_binomial))
    FetchUtil.write_fasta(source_protein)
    file_name = source_protein.accession + '_' + ''.join(w[:2] for w in org.split())
    print("Using {} as filename".format(file_name))
    outfile = XML_PATH + file_name + '.xml'
    if not os.path.exists(outfile):
        FetchUtil.remote_blast(ORTHOS_PATH + source + '.fasta', thresh, outfile, org)
    else:
        print('Outfile exists')
    with open(outfile) as qoutput:
        parser = NCBIXML.parse(qoutput)
        for lin in parser:
            for align in lin.alignments:
                for hsp in align.hsps:
                    if (hsp.positives / float(hsp.align_length)) >= .4 and (
                            float(hsp.align_length) / len(hsp.query)) >= .25:
                        ac.append(align.accession)
    print("First BLAST Done. Number of sequences found: " + str(ac))

    for o in ac:
        print("BLASTING back to " + source_binomial)
        print(o)
        o_prot = FetchUtil.fetch_protein(o)
        acc = []
        FetchUtil.write_fasta(o_prot)
        file_name = o_prot.accession + '_' + ''.join(w[:2] for w in source_binomial.split())
        print("Using " + file_name + " as a new filename")
        outfile = XML_PATH + file_name + '.xml'
        if not os.path.exists(outfile):
            FetchUtil.remote_blast(ORTHOS_PATH + o_prot.accession + '.fasta', thresh, outfile, source_binomial)
        else:
            print('Outfile exists')
        with open(outfile) as q1output:
            parse = NCBIXML.parse(q1output)
            print('blasted')
            for lin in parse:
                for align in lin.alignments:
                    for hsp in align.hsps:
                        if (hsp.positives / float(hsp.align_length)) >= .4 and (
                                float(hsp.align_length) / len(hsp.query)) > .25:
                            acc.append(align.accession)
        print("Done. Number of sequences found: " + repr(len(acc)))

        if source in acc:
            print("It's twue!")
            name = o_prot.binomial
            acclist[name] = [
                o,
                str(ac.index(o) + 1) + '/' + str(len(ac)),
                str(acc.index(source) + 1) + '/' + str(len(acc))
            ]
            with open(DICTS_PATH + source, 'a') as seed_file:
                seed_file.write(str(acclist) + '\n')
            break
    return acclist

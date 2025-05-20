import os.path

from Utils import FetchUtil
from Utils.FileUtil import xml_parse_and_extract_accession_numbers
from helpers import commands
from helpers.constants import XMLPath, DictsPath


def best_reciprocal_blast(org, seed, thresh=5):
    """Returns the best pairwise reciprocal BLAST using seed accession no. from against org organism
    @returns {
    Organism binomial name:
    [Accession number, rank of search in target organism, rank of search in source organism]
    }
    """
    seedorg = FetchUtil.fetch_organism(seed)[0]
    FetchUtil.fetch_fasta(seed)
    dum = '_'.join((seed+seedorg+org).split())
    print("Run: " + dum)
    if not os.path.isfile(str(XMLPath(dum + '.xml'))) or not os.path.getsize(str(XMLPath(dum + '.xml'))):
        print('blasting')
        commands.run_blast(seed, thresh, dum, org)
    ac = xml_parse_and_extract_accession_numbers(str(XMLPath(dum + '.xml')))
    print("Done. Number of sequences found: " + repr(len(ac)))
    acclist = {}
    for o in ac:
        print(o)
        if len(o) <= 4:
            print('Skipping')
            continue
        FetchUtil.fetch_fasta(o)
        dum2 = '_'.join((o+org+seedorg).split())
        if not os.path.isfile(str(XMLPath(dum2 + '.xml'))) or not os.path.getsize(str(XMLPath(dum2 + '.xml'))):
            print('blasting back')
            commands.run_blast(o, thresh, dum2, seedorg)
        acc = xml_parse_and_extract_accession_numbers(str(XMLPath(dum2 + '.xml')))
        print("Done. Number of sequences found: " + repr(len(acc)))

        if seed in acc:
            print("it's twue!")
            name = FetchUtil.fetch_organism(o)[0]
            acclist[name] = [o, str(ac.index(o) + 1) + '/' + str(len(ac)),
                             str(acc.index(seed) + 1) + '/' + str(len(acc))]
            with open(str(DictsPath(seed)), 'a') as dicts:
                dicts.write(str(acclist) + '\n')
            break
    return acclist

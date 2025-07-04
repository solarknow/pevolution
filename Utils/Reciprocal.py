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
    run_id = '_'.join((seed+seedorg+org).split())
    print("Run: " + run_id)
    if not os.path.isfile(str(XMLPath(run_id + '.xml'))) or not os.path.getsize(str(XMLPath(run_id + '.xml'))):
        print('blasting')
        commands.run_blast(seed, thresh, run_id, org)
    ac = xml_parse_and_extract_accession_numbers(str(XMLPath(run_id + '.xml')))
    print("Done. Number of sequences found: " + repr(len(ac)))
    acclist = {}
    for o in ac:
        print(o)
        if len(o) <= 4:
            print('Skipping')
            continue
        FetchUtil.fetch_fasta(o)
        rev_run_id = '_'.join((o+org+seedorg).split())
        if not os.path.isfile(str(XMLPath(rev_run_id + '.xml'))) or not os.path.getsize(str(XMLPath(rev_run_id + '.xml'))):
            print('blasting back')
            commands.run_blast(o, thresh, rev_run_id, seedorg)
        acc = xml_parse_and_extract_accession_numbers(str(XMLPath(rev_run_id + '.xml')))
        print("Done. Number of sequences found: " + repr(len(acc)))

        if seed in acc:
            print("The original query sequence was found!")
            name = FetchUtil.fetch_organism(o)[0]
            acclist[name] = [o, str(ac.index(o) + 1) + '/' + str(len(ac)),
                             str(acc.index(seed) + 1) + '/' + str(len(acc))]
            with open(str(DictsPath(seed)), 'a') as dicts:
                dicts.write(str(acclist) + '\n')
            break
    return acclist

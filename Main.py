import getopt
import os
import sys

from Utils import SeqUtil, FetchUtil, FileUtil, Reciprocal
from Utils.dom import DomainRun
from helpers.constants import DataPath, domain_thresholds, organism_list, Domains


def main(argv):
    # Expected: 1=query accession no.; 2=out prefix; 3=domain (euk,bac,arch,all)
    query = ''
    out = ''
    dom = ''
    phy = False

    try:
        opts, args = getopt.getopt(argv, 'q:o:d:e:yh', ['query=', 'output=', 'domain=', 'email=', 'with_phyml'])
    except getopt.GetoptError:
        print('Main.py -q|--query <query accession number> -o|--output <output name> -d|--domain <domain of life to '
              'query> -e|--email <email of user> [-y | --with_phyml]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Main.py -q|--query <query accession number> -o|--output <output name> -d|--domain <domain of life '
                  'to query> -e|--email <email of user> [-y | --with_phyml]')
            sys.exit()
        elif opt in ('-q', '--query'):
            query = arg
        elif opt in ('-o', '--output'):
            out = arg
        elif opt in ('-d', '--domain'):
            dom = arg
        elif opt in ('-y', '--with_phyml'):
            phy = True
        elif opt in ('-e', '--email'):
            FetchUtil.set_email(arg)
    if not os.path.exists(str(DataPath(dom + '-' + out + '.fas'))):
        dom_query = FetchUtil.fetch_organism(query)[1]
        thresh = domain_thresholds(dom_query)

        arch_accs = {}
        bac_accs = {}
        euk_accs = {}

        print("Blasting")
        if dom == 'arch' or dom == 'all':
            for a in organism_list(Domains.ARCHAEA):
                arch_accs.update(Reciprocal.best_reciprocal_blast(a, query, thresh.arch))
        if dom == 'bac' or dom == 'all':
            for b in organism_list(Domains.BACTERIA):
                bac_accs.update(Reciprocal.best_reciprocal_blast(b, query, thresh.bac))
        if dom == 'euk' or dom == 'all':
            for e in organism_list(Domains.EUKARYOTA):
                euk_accs.update(Reciprocal.best_reciprocal_blast(e, query, thresh.euk))

        all_accs = dict(list(arch_accs.items()) + list(bac_accs.items()) + list(euk_accs.items()))
        if dom == 'all':
            num_seqs = sum(len(v) for v in all_accs.values())
            print(f"Dictionary generated with {len(all_accs)} keys and {num_seqs} sequences.")
        # Fetching the sequences and writing them to file
        print("Writing seqs to file.")
        if arch_accs:
            FileUtil.merge_domain_fastas('arch-' + out + '.fas', arch_accs)
            SeqUtil.append_sequences(DataPath('all-' + out + '.fas'), DataPath('arch-' + out + '.fas'))
        if bac_accs:
            FileUtil.merge_domain_fastas('bac-' + out + '.fas', bac_accs)
            SeqUtil.append_sequences(DataPath('all-' + out + '.fas'), DataPath('bac-' + out + '.fas'))
        if euk_accs:
            FileUtil.merge_domain_fastas('euk-' + out + '.fas', euk_accs)
            SeqUtil.append_sequences(DataPath('all-' + out + '.fas'), DataPath('euk-' + out + '.fas'))

    DomainRun(out, query, dom, phy).run_report()


if __name__ == "__main__":
    main(sys.argv[1:])

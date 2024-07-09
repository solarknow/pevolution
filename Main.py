import getopt
import os
import sys
from multiprocessing import Process, Queue

import psutil

import FetchUtil
import Reciprocal
import SeqUtil
from helpers.commands import run_dom_file
from helpers.constants import DataPath, OrthosPath

PROCS = psutil.cpu_count()


def main(argv):
    # Expected: 1=query accession no.; 2=out prefix; 3=domain (euk,bac,arch,all)
    query = ''
    out = ''
    dom = ''
    phy = False

    try:
        opts, args = getopt.getopt(argv, 'q:o:d:yh', ['query=', 'output=', 'domain=', 'with_phyml'])
    except getopt.GetoptError:
        print('Main.py -q|--query <query accession number> -o|--output <output name> -d|--domain <domain of life to '
              'query> [-y | --with_phyml]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Main.py -q|--query <query accession number> -o|--output <output name> -d|--domain <domain of life '
                  'to query> [-y | --with_phyml]')
            sys.exit()
        elif opt in ('-q', '--query'):
            query = arg
        elif opt in ('-o', '--output'):
            out = arg
        elif opt in ('-d', '--domain'):
            dom = arg
        elif opt in ('-y', '--with_phyml'):
            phy = True
    if not os.path.exists(str(DataPath(dom + '-' + out + '.fas'))):
        arch_list = ['Haloferax volcanii', 'Sulfolobus tokodaii', 'Methanococcus aeolicus',
                     'Methanobrevibacter smithii', 'Thermococcus sibiricus', 'Archaeoglobus fulgidus',
                     'Nanoarchaeum equitans', 'Thermoplasma acidophilum']
        bac_list = ['Gemmata obscuriglobus', 'Prosthecobacter dejongeii', 'Verrucomicrobium spinosum',
                    'Rickettsia prowazekii', 'Agrobacterium tumefaciens', 'Escherichia coli', 'Bacillus subtilis',
                    'Anabaena variabilis', 'Thermotoga maritima']
        euk_list = ['Drosophila melanogaster', 'Homo sapiens', 'Oryza sativa', 'Trypanosoma brucei',
                    'Plasmodium falciparum', 'Saccharomyces cerevisiae', 'Neurospora crassa', 'Arabidopsis thaliana']
        # subject to change
        # setting threshold values: thresh1-w/ arch ;thresh2-w/ bac;
        dom_query = FetchUtil.fetch_organism(query)[1]
        if dom_query == 'Archaea':
            thresh1 = 1e-10
            thresh2 = 1e-5
            thresh3 = 5
        elif dom_query == 'Eukaryota':
            thresh1 = 5
            thresh2 = 5
            thresh3 = 1e-10
        else:
            thresh1 = 1e-5
            thresh2 = 1e-10
            thresh3 = 5

        arch_accs = {}
        bac_accs = {}
        euk_accs = {}
        print("Blasting")
        if dom == 'arch' or dom == 'all':
            queue_arch = Queue()
            for a in arch_list:
                p = Process(target=Reciprocal.best_reciprocal_blast, args=(a, query, thresh1,
                                                                           queue_arch))
                p.start()
                p.join()
            while not queue_arch.empty():
                arch_accs.update(queue_arch.get())
        if dom == 'bac' or dom == 'all':
            queue_bac = Queue()
            for b in bac_list:
                p = Process(target=Reciprocal.best_reciprocal_blast, args=(b, query, thresh2,
                                                                           queue_bac))
                p.start()
                p.join()
            while not queue_bac.empty():
                bac_accs.update(queue_bac.get())
        if dom == 'euk' or dom == 'all':
            queue_euk = Queue()
            for e in euk_list:
                p = Process(target=Reciprocal.best_reciprocal_blast, args=(e, query, thresh3,
                                                                           queue_euk))
                p.start()
                p.join()
            while not queue_euk.empty():
                euk_accs.update(queue_euk.get())

        all_accs = {}
        if dom == all:
            all_accs.update(arch_accs)
            all_accs.update(bac_accs)
            all_accs.update(euk_accs)
            num_seqs = 0
            for j in all_accs:
                num_seqs += len(all_accs[j])
            print(f"Dictionary generated with {len(all_accs)} keys and {num_seqs} sequences.")
        # Fetching the sequences and writing them to file
        print("Writing seqs to file.")
        for a in arch_accs:
            FetchUtil.fetch_fasta(arch_accs[a][0])
            with open(str(OrthosPath(arch_accs[a][0] + '.fasta'))) as fil:
                fil_arr = fil.readlines()

            with open(str(OrthosPath(arch_accs[a][0] + '.fasta')), 'w') as fil:
                for i in range(len(fil_arr)):
                    if i == 0:
                        fil.write(fil_arr[i].strip() + ' ' + arch_accs[a][1] + '  ' + arch_accs[a][2] + '\n')
                    else:
                        fil.write(fil_arr[i])

            SeqUtil.addseq(DataPath('arch-' + out + '.fas'), OrthosPath(arch_accs[a][0] + '.fasta'))
            os.remove(str(OrthosPath(arch_accs[a][0] + '.fasta')))
        for b in bac_accs:
            FetchUtil.fetch_fasta(bac_accs[b][0])
            with open(str(OrthosPath(bac_accs[b][0] + '.fasta'))) as fil:
                fil_arr = fil.readlines()

            with open(str(OrthosPath(bac_accs[b][0] + '.fasta')), 'w') as fil:
                for i in range(len(fil_arr)):
                    if i == 0:
                        fil.write(fil_arr[i].strip() + ' ' + bac_accs[b][1] + '  ' + bac_accs[b][2] + '\n')
                    else:
                        fil.write(fil_arr[i])

            SeqUtil.addseq(DataPath('bac-' + out + '.fas'), OrthosPath(bac_accs[b][0] + '.fasta'))
            os.remove(str(OrthosPath(bac_accs[b][0] + '.fasta')))
        for e in euk_accs:
            FetchUtil.fetch_fasta(euk_accs[e][0])
            with open(str(OrthosPath(euk_accs[e][0] + '.fasta'))) as fil:
                fil_arr = fil.readlines()

            with open(str(OrthosPath(euk_accs[e][0] + '.fasta')), 'w') as fil:
                for i in range(len(fil_arr)):
                    if i == 0:
                        fil.write(fil_arr[i].strip() + ' ' + euk_accs[e][1] + '  ' + euk_accs[e][2] + '\n')
                    else:
                        fil.write(fil_arr[i])

            SeqUtil.addseq(DataPath('euk-' + out + '.fas'), OrthosPath(euk_accs[e][0] + '.fasta'))
            os.remove(str(OrthosPath(euk_accs[e][0] + '.fasta')))
        SeqUtil.addseq(DataPath('all-' + out + '.fas'), DataPath('arch-' + out + '.fas'))
        SeqUtil.addseq(DataPath('all-' + out + '.fas'), DataPath('bac-' + out + '.fas'))
        SeqUtil.addseq(DataPath('all-' + out + '.fas'), DataPath('euk-' + out + '.fas'))
        for c in all_accs:
            FetchUtil.fetch_fasta(all_accs[c][0])
            with open(str(OrthosPath(all_accs[c][0] + '.fasta'))) as fil:
                fil_arr = fil.readlines()

            with open(str(OrthosPath(all_accs[c][0] + '.fasta')), 'w') as fil:
                for i in range(len(fil_arr)):
                    if i == 0:
                        fil.write(fil_arr[i].strip() + ' ' + all_accs[c][1] + '  ' + all_accs[c][2] + '\n')
                    else:
                        fil.write(fil_arr[i])

            SeqUtil.addseq(DataPath('all-' + out + '.fas'), OrthosPath(all_accs[c][0] + '.fasta'))
            os.remove(str(OrthosPath(all_accs[c][0] + '.fasta')))

    run_dom_file(out, query, dom, phy)


if __name__ == "__main__":
    main(sys.argv[1:])

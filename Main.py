import os
import sys

import psutil
from Bio import Entrez

import AlignUtil
import FetchUtil
import Reciprocal
import SeqUtil
from FetchUtil import ORTHOS_PATH
from constants import DATA_PATH, ALIGNS_PATH

PROCS = psutil.cpu_count()

os.makedirs('Data', exist_ok=True)
os.makedirs('Orthos', exist_ok=True)

# Expected: 1=query accession no.; 2=out prefix; 3=domain (euk,bac,arch,all)
phyml = ''
if '-y' in sys.argv:
    phyml = '-y'
try:
    query = sys.argv[1]
    out = sys.argv[2]
    dom = sys.argv[3]
    if Entrez.email is None:
        FetchUtil.set_email(input('Email: '))
except IndexError:
    query = input('Query accession no: ')
    out = input('Output file: ')
    dom = input('Organism Domain to be explored: euk,bac,arch, or all ').lower()
    phy = input("Should we use PhyML to find trees, in addition to MrBayes? [y/N]")
    FetchUtil.set_email(input('Email: '))
    if phy.lower() == 'y':
        phyml = '-y'
if os.path.exists(DATA_PATH + dom + '-' + out + '.fas'):
    AlignUtil.align_sequences_prank(out, dom)
    models = SeqUtil.prottest_best_models(ALIGNS_PATH + dom + '-' + out + '.best.nex')
    if phyml == '-y':
        AlignUtil.run_phyml(out, dom, models)
    AlignUtil.prepare_bayes_files(out, dom, models)
    AlignUtil.run_bayes_and_report(out, dom, FetchUtil.fetch_protein(query), models)
else:
    arch_list = ['Haloferax volcanii', 'Sulfolobus tokodaii', 'Methanococcus aeolicus', 'Methanobrevibacter smithii',
                 'Thermococcus sibiricus', 'Archaeoglobus fulgidus', 'Nanoarchaeum equitans',
                 'Thermoplasma acidophilum']
    bac_list = ['Gemmata obscuriglobus', 'Prosthecobacter dejongeii', 'Verrucomicrobium spinosum',
                'Rickettsia prowazekii', 'Agrobacterium tumefaciens', 'Escherichia coli', 'Bacillus subtilis',
                'Anabaena variabilis', 'Thermotoga maritima']
    euk_list = ['Drosophila melanogaster', 'Homo sapiens', 'Oryza sativa', 'Trypanosoma brucei',
                'Plasmodium falciparum', 'Saccharomyces cerevisiae', 'Neurospora crassa',
                'Arabidopsis thaliana']  # subject to change
    # setting threshold values: thresh1-w/ arch ;thresh2-w/ bac;
    dom_query = FetchUtil.fetch_protein(query).domain
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

    accs = {}
    print("Blasting")
    if dom == 'arch' or dom == 'all':
        print("Archea")
        for a in arch_list:
            p = Reciprocal.bestrecipblast(a, query, thresh1)
            print(p)
            if 'arch' in accs:
                accs['arch'].update(p)
            else:
                accs['arch'] = p
    if dom == 'bac' or dom == 'all':
        print("Bacteria")
        for b in bac_list:
            p = Reciprocal.bestrecipblast(b, query, thresh2)
            if 'bac' in accs:
                accs['bac'].update(p)
            else:
                accs['bac'] = p
    if dom == 'euk' or dom == 'all':
        print("Eukarya")
        for e in euk_list:
            p = Reciprocal.bestrecipblast(e, query, thresh3)
            if 'euk' in accs:
                accs['euk'].update(p)
            else:
                accs['euk'] = p

    if dom == 'all':
        print("All")
        all_accs = {}
        for v in accs.values():
            all_accs.update(v)
        num_seqs = sum(len(all_accs[j]) for j in all_accs)
        print("Dictionary generated with " + repr(len(all_accs)) + " keys and " + repr(num_seqs) + " sequences.")
    # Fetching the sequences and writing them to file
    print("Writing seqs to file.")
    print(accs)
    for d, daccs in accs.items():
        print(daccs)
        for seq in daccs.values():
            FetchUtil.write_fasta(seq[0])
            with open(ORTHOS_PATH + seq[0] + '.fasta') as fil_arr:
                with open(ORTHOS_PATH + seq[0] + '.fasta', 'w') as fil:
                    for i in fil_arr:
                        if i.startswith('>'):
                            fil.write(i.strip() + ' ' + seq[1] + '  ' + seq[2] + '\n')
                        else:
                            fil.write(i)
            SeqUtil.fasta_add_sequence(DATA_PATH + d + '-' + out + '.fas', ORTHOS_PATH + seq[0] + '.fasta')
        SeqUtil.fasta_add_sequence(DATA_PATH + 'all-' + out + '.fas', DATA_PATH + d + '-' + out + '.fas')

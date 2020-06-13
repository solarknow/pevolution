import os
import sys


from Bio import Entrez

import AlignUtil
import FetchUtil
import Reciprocal
import SeqUtil
from FetchUtil import ORTHOS_PATH
from constants import DATA_PATH, ALIGNS_PATH


os.makedirs('Data', exist_ok=True)
os.makedirs('Orthos', exist_ok=True)

# Expected: 1=query accession no.; 2=out prefix; 3=domain (euk,bac,arch,all)
local = '-y' in sys.argv
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
    loc = input('Run locally? [y/N]')
    FetchUtil.set_email(input('Email: '))
    local = loc.lower() == 'y'
if os.path.exists(DATA_PATH + dom + '-' + out + '.fas'):
    AlignUtil.align_sequences_prank(out, dom)
    models = SeqUtil.prottest_best_models(ALIGNS_PATH + dom + '-' + out + '.best.nex')
    AlignUtil.prepare_bayes_files(out, dom, models)
    AlignUtil.run_bayes_and_report(out, dom, FetchUtil.fetch_protein(query), models)
else:
    arch_list = {'309800': 'Haloferax volcanii DS2', '273063': 'Sulfurisphaera tokodaii str. 7',
                 '419665': 'Methanococcus aeolicus Nankai-3', '2173': 'Methanobrevibacter smithii',
                 '604354': 'Thermococcus sibiricus MM 739', '2234': 'Archaeoglobus fulgidus',
                 '228908': 'Nanoarchaeum equitans Kin4-M', '273075': 'Thermoplasma acidophilum DSM 1728'}
    bac_list = {'214688': 'Gemmata obscuriglobus UQM 2246', '48465': 'Prosthecobacter dejongeii',
                '2736': 'Verrucomicrobium spinosum', '782': 'Rickettsia prowazekii', '358': 'Agrobacterium tumefaciens',
                '562': 'Escherichia coli', '1423': 'Bacillus subtilis', '264691': 'Anabaena variabilis',
                '243274': 'Thermotoga maritima MSB8'}
    euk_list = {'7227': 'Drosophila melanogaster', '9606': 'Homo sapiens', '4530': 'Oryza sativa',
                '5691': 'Trypanosoma brucei', '5833': 'Plasmodium falciparum', '4932': 'Saccharomyces cerevisiae',
                '5141': 'Neurospora crassa', '3702': 'Arabidopsis thaliana'}  # subject to change
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
        for a, binom in arch_list.items():
            p = Reciprocal.bestrecipblast([a, binom], query, thresh1, local)
            if 'arch' in accs:
                accs['arch'].update(p)
            else:
                accs['arch'] = p
    if dom == 'bac' or dom == 'all':
        print("Bacteria")
        for b, binom in bac_list.items():
            p = Reciprocal.bestrecipblast([b, binom], query, thresh2, local)
            if 'bac' in accs:
                accs['bac'].update(p)
            else:
                accs['bac'] = p
    if dom == 'euk' or dom == 'all':
        print("Eukarya")
        for e, binom in euk_list.items():
            p = Reciprocal.bestrecipblast([e, binom], query, thresh3, local)
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
            FetchUtil.write_fasta(FetchUtil.fetch_protein(seq[0]))
            with open(ORTHOS_PATH + seq[0] + '.fasta') as fil_arr:
                with open(ORTHOS_PATH + seq[0] + '.fas', 'w') as fil:
                    for i in fil_arr:
                        if i.startswith('>'):
                            fil.write(i.strip() + ' ' + seq[1] + '  ' + seq[2] + '\n')
                        else:
                            fil.write(i)
            SeqUtil.fasta_add_sequence(DATA_PATH + d + '-' + out + '.fas', ORTHOS_PATH + seq[0] + '.fas')
        SeqUtil.fasta_add_sequence(DATA_PATH + 'all-' + out + '.fas', DATA_PATH + d + '-' + out + '.fas')

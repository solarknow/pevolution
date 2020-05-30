import os
import subprocess
import sys

import psutil

import FetchUtil
import Reciprocal
import SeqUtil
from FetchUtil import ORTHOS_PATH
from constants import DATA_PATH

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
    subprocess.call(['python', dom + '.py', out, query, phyml])
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

    arch_accs = {}
    bac_accs = {}
    euk_accs = {}
    print("Blasting")
    if dom == 'arch' or dom == 'all':
        for a in arch_list:
            p = Reciprocal.bestrecipblast(a, query, thresh1)
            arch_accs.update(p)
    if dom == 'bac' or dom == 'all':
        # bac_accs=Reciprocal.bestrecipblast(bac_list,query,thresh2)
        for b in bac_list:
            p = Reciprocal.bestrecipblast(b, query, thresh2)
            bac_accs.update(p)
    if dom == 'euk' or dom == 'all':
        # euk_accs=Reciprocal.bestrecipblast(euk_list,query,thresh3)
        for e in euk_list:
            p = Reciprocal.bestrecipblast(e, query, thresh3)
            euk_accs.update(p)

    if dom == 'all':
        all_accs = {}
        all_accs.update(arch_accs)
        all_accs.update(bac_accs)
        all_accs.update(euk_accs)
        num_seqs = sum(len(all_accs[j]) for j in all_accs)
        print("Dictionary generated with " + repr(len(all_accs)) + " keys and " + repr(num_seqs) + " sequences.")
    # Fetching the sequences and writing them to file
    print("Writing seqs to file.")
    for a in arch_accs:
        FetchUtil.write_fasta(arch_accs[a][0])
        with open(ORTHOS_PATH + arch_accs[a][0] + '.fasta') as fil:
            fil_arr = fil.readlines()
        with open(ORTHOS_PATH + arch_accs[a][0] + '.fasta', 'w') as fil:
            for i in range(len(fil_arr)):
                if i == 0:
                    fil.write(fil_arr[i].strip() + ' ' + arch_accs[a][1] + '  ' + arch_accs[a][2] + '\n')
                else:
                    fil.write(fil_arr[i])
        SeqUtil.fasta_add_sequence(DATA_PATH + 'arch-' + out + '.fas', ORTHOS_PATH + arch_accs[a][0] + '.fasta')
        # os.remove(ORTHOS_PATH + arch_accs[a][0] + '.fasta')
    for b in bac_accs.keys():
        FetchUtil.write_fasta(bac_accs[b][0])
        fil = open(ORTHOS_PATH + bac_accs[b][0] + '.fasta')
        fil_arr = fil.readlines()
        fil.close()
        fil = open(ORTHOS_PATH + bac_accs[b][0] + '.fasta', 'w')
        for i in range(len(fil_arr)):
            if i == 0:
                fil.write(fil_arr[i].strip() + ' ' + bac_accs[b][1] + '  ' + bac_accs[b][2] + '\n')
            else:
                fil.write(fil_arr[i])
        fil.close()
        SeqUtil.fasta_add_sequence(DATA_PATH + 'bac-' + out + '.fas', ORTHOS_PATH + bac_accs[b][0] + '.fasta')
        # os.remove(ORTHOS_PATH + bac_accs[b][0] + '.fasta')
    for e in euk_accs:
        FetchUtil.write_fasta(euk_accs[e][0])
        with open(ORTHOS_PATH + euk_accs[e][0] + '.fasta') as fil:
            fil_arr = fil.readlines()
        with open(ORTHOS_PATH + euk_accs[e][0] + '.fasta', 'w') as fil:
            for i in range(len(fil_arr)):
                if i == 0:
                    fil.write(fil_arr[i].strip() + ' ' + euk_accs[e][1] + '  ' + euk_accs[e][2] + '\n')
                else:
                    fil.write(fil_arr[i])
        SeqUtil.fasta_add_sequence(DATA_PATH + 'euk-' + out + '.fas', ORTHOS_PATH + euk_accs[e][0] + '.fasta')
        # os.remove(ORTHOS_PATH + euk_accs[e][0] + '.fasta')
    SeqUtil.fasta_add_sequence(DATA_PATH + 'all-' + out + '.fas', ORTHOS_PATH + 'arch-' + out + '.fas')
    SeqUtil.fasta_add_sequence(DATA_PATH + 'all-' + out + '.fas', ORTHOS_PATH + 'bac-' + out + '.fas')
    SeqUtil.fasta_add_sequence(DATA_PATH + 'all-' + out + '.fas', ORTHOS_PATH + 'euk-' + out + '.fas')



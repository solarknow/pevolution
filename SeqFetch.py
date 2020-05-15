import FetchUtil
import Reciprocal
import SeqUtil
import os
import sys

out = sys.argv[1]
query = sys.argv[2]

try:
    open('arch-' + out + '.fas')
    open('bac-' + out + '.fas')
    open('euk-' + out + '.fas')
    open('all-' + out + '.fas')
except IOError:
    arch_list = ['Haloferax volcanii', 'Sulfolobus tokodaii', 'Methanococcus aeolicus', 'Methanobrevibacter smithii',
                 'Thermococcus sibiricus', 'Archaeoglobus fulgidus', 'Nanoarchaeum equitans',
                 'Thermoplasma acidophilum']
    bac_list = ['Gemmata obscuriglobus', 'Prosthecobacter dejongeii', 'Verrucomicrobium spinosum',
                'Rickettsia prowazekii', 'Agrobacterium tumefaciens', 'Escherichia coli',
                'Bacillus subtilis', 'Anabaena variabilis', 'Thermotoga maritima']
    euk_list = ['Drosophila melanogaster', 'Homo sapiens', 'Oryza sativa',
                'Trypanosoma brucei', 'Plasmodium falciparum', 'Saccharomyces cerevisiae',
                'Neurospora crassa', 'Arabidopsis thaliana']  # subject to change
    # setting threshold values; thresh1-
    dom = FetchUtil.fetch_organism(query)[1]
    if dom == 'Archaea':
        thresh1 = 1e-10
        thresh2 = 1e-5
        thresh3 = 5
    elif dom == 'Eukaryota':
        thresh1 = 5
        thresh2 = 5
        thresh3 = 1e-10
    else:
        thresh1 = 1e-5
        thresh2 = 1e-10
        thresh3 = 5
    os.makedirs('XML', exist_ok=True)
    os.makedirs('Orthos', exist_ok=True)
    arch_accs = Reciprocal.bestrecipblast(arch_list, query, thresh1)
    bac_accs = Reciprocal.bestrecipblast(bac_list, query, thresh2)
    euk_accs = Reciprocal.bestrecipblast(euk_list, query, thresh3)
    all_accs = {}
    all_accs.update(arch_accs)
    all_accs.update(bac_accs)
    all_accs.update(euk_accs)

    print(all_accs)
    # Fetching the sequences and writing them to file
    for a in arch_accs.keys():
        FetchUtil.fetch_fasta(arch_accs[a][0])
        with open('Orthos' + os.sep + arch_accs[a][0] + '.fasta') as fil:
            fil_arr = fil.readlines()
        with open('Orthos' + os.sep + arch_accs[a][0] + '.fasta', 'w') as fil:
            for i in range(len(fil_arr)):
                if i == 0:
                    fil.write(fil_arr[i].strip() + ' ' + arch_accs[a][1] + '  ' + arch_accs[a][2] + '\n')
                else:
                    fil.write(fil_arr[i])
        SeqUtil.addseq('arch-' + out + '.fas', 'Orthos' + os.sep + arch_accs[a][0] + '.fasta')

    for accs in bac_accs.values():
        FetchUtil.fetch_fasta(accs[0])
        with open('Orthos' + os.sep + accs[0] + '.fasta') as fil:
            fil_arr = fil.readlines()
        with open('Orthos' + os.sep + accs[0] + '.fasta', 'w') as fil:
            for i in range(len(fil_arr)):
                if i == 0:
                    fil.write(fil_arr[i].strip() + ' ' + accs[1] + '  ' + accs[2] + '\n')
                else:
                    fil.write(fil_arr[i])
        SeqUtil.addseq('bac-' + out + '.fas', 'Orthos' + os.sep + accs[0] + '.fasta')

    for e in euk_accs.values():
        FetchUtil.fetch_fasta(e[0])
        fil = open('Orthos' + os.sep + e[0] + '.fasta')
        fil_arr = fil.readlines()
        fil.close()
        fil = open('Orthos' + os.sep + e[0] + '.fasta', 'w')
        for i in range(len(fil_arr)):
            if i == 0:
                fil.write(fil_arr[i].strip() + ' ' + e[1] + '  ' + e[2] + '\n')
            else:
                fil.write(fil_arr[i])
        fil.close()
        SeqUtil.addseq('euk-' + out + '.fas', 'Orthos\\' + e[0] + '.fasta')

    for c in all_accs.keys():
        FetchUtil.fetch_fasta(all_accs[c][0])
        fil = open('Orthos\\' + all_accs[c][0] + '.fasta')
        fil_arr = fil.readlines()
        fil.close()
        fil = open('Orthos\\' + all_accs[c][0] + '.fasta', 'w')
        for i in range(len(fil_arr)):
            if i == 0:
                fil.write(fil_arr[i].strip() + ' ' + all_accs[c][1] + '  ' + all_accs[c][2] + '\n')
            else:
                fil.write(fil_arr[i])
        fil.close()
        SeqUtil.addseq('all-' + out + '.fas', 'Orthos\\' + all_accs[c][0] + '.fasta')

import os
import psutil
import sys
from multiprocessing import Process, Queue

import FetchUtil
import Reciprocal
import SeqUtil

PROCS = psutil.cpu_count()
# Check for directories' existance, if not create them.
if not os.path.exists('Data'):
    os.mkdir('Data')
if not os.path.exists('Orthos'):
    os.mkdir('Orthos')

# Expected: 1=query accession no.; 2=out prefix; 3=domain (euk,bac,arch,all)
phyml = ''
if '-y' in sys.argv:
    phyml = '-y'
try:
    query = sys.argv[1]
    out = sys.argv[2]
    dom = sys.argv[3]
except IndexError:
    query = input('Query accession no: ')
    out = input('Output file: ')
    dom = input('Organism Domain to be explored: euk,bac,arch, or all ').lower()
    phy = input("Should we use PhyML to find trees, in addition to MrBayes? [y/N]")
    if phy.lower() == 'y':
        phyml = '-y'
try:
    open('Data/' + dom + '-' + out + '.fas').read()
except IOError:
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
        #       p.join()
        while not queue_arch.empty():
            arch_accs.update(queue_arch.get())
    if dom == 'bac' or dom == 'all':
        # bac_accs=Reciprocal.best_reciprocal_blast(bac_list,query,thresh2)
        queue_bac = Queue()
        for b in bac_list:
            p = Process(target=Reciprocal.best_reciprocal_blast, args=(b, query, thresh2,
                                                                       queue_bac))
            p.start()
            p.join()
        while not queue_bac.empty():
            bac_accs.update(queue_bac.get())
    if dom == 'euk' or dom == 'all':
        # euk_accs=Reciprocal.best_reciprocal_blast(euk_list,query,thresh3)
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
        print("Dictionary generated with " + repr(len(all_accs.keys())) + " keys and " + repr(num_seqs) + " sequences.")
    # Fetching the sequences and writing them to file
    print("Writing seqs to file.")
    for a in arch_accs.keys():
        FetchUtil.fetch_fasta(arch_accs[a][0])
        fil = open('Orthos/' + arch_accs[a][0] + '.fasta')
        fil_arr = fil.readlines()
        fil.close()
        fil = open('Orthos/' + arch_accs[a][0] + '.fasta', 'w')
        for i in range(len(fil_arr)):
            if i == 0:
                fil.write(fil_arr[i].strip() + ' ' + arch_accs[a][1] + '  ' + arch_accs[a][2] + '\n')
            else:
                fil.write(fil_arr[i])
        fil.close()
        SeqUtil.addseq('Data/arch-' + out + '.fas', 'Orthos/' + arch_accs[a][0] + '.fasta')
        os.remove('Orthos/' + arch_accs[a][0] + '.fasta')
    for b in bac_accs.keys():
        FetchUtil.fetch_fasta(bac_accs[b][0])
        fil = open('Orthos/' + bac_accs[b][0] + '.fasta')
        fil_arr = fil.readlines()
        fil.close()
        fil = open('Orthos/' + bac_accs[b][0] + '.fasta', 'w')
        for i in range(len(fil_arr)):
            if i == 0:
                fil.write(fil_arr[i].strip() + ' ' + bac_accs[b][1] + '  ' + bac_accs[b][2] + '\n')
            else:
                fil.write(fil_arr[i])
        fil.close()
        SeqUtil.addseq('Data/bac-' + out + '.fas', 'Orthos/' + bac_accs[b][0] + '.fasta')
        os.remove('Orthos/' + bac_accs[b][0] + '.fasta')
    for e in euk_accs.keys():
        FetchUtil.fetch_fasta(euk_accs[e][0])
        fil = open('Orthos/' + euk_accs[e][0] + '.fasta')
        fil_arr = fil.readlines()
        fil.close()
        fil = open('Orthos/' + euk_accs[e][0] + '.fasta', 'w')
        for i in range(len(fil_arr)):
            if i == 0:
                fil.write(fil_arr[i].strip() + ' ' + euk_accs[e][1] + '  ' + euk_accs[e][2] + '\n')
            else:
                fil.write(fil_arr[i])
        fil.close()
        SeqUtil.addseq('Data/euk-' + out + '.fas', 'Orthos/' + euk_accs[e][0] + '.fasta')
        os.remove('Orthos/' + euk_accs[e][0] + '.fasta')
    SeqUtil.addseq('Data/all-' + out + '.fas', 'Data/arch-' + out + '.fas')
    SeqUtil.addseq('Data/all-' + out + '.fas', 'Data/bac-' + out + '.fas')
    SeqUtil.addseq('Data/all-' + out + '.fas', 'Data/euk-' + out + '.fas')
#    for c in all_accs.keys():
#        Fetchutil.seqfetch(all_accs[c][0])
#        fil=open('Orthos/'+all_accs[c][0]+'.fasta')
#        fil_arr=fil.readlines()
#        fil.close()
#        fil=open('Orthos/'+all_accs[c][0]+'.fasta','w')
#        for i in range(len(fil_arr)):
#            if i==0:
#                fil.write(fil_arr[i].strip()+' '+all_accs[c][1]+'  '+all_accs[c][2]+'\n')
#            else:
#                fil.write(fil_arr[i])
#        fil.close()
#        Fetchutil.addseq('Data/all-'+out+'.fas','Orthos/'+all_accs[c][0]+'.fasta')
#        os.remove('Orthos/'+all_accs[c][0]+'.fasta')

os.system('python ' + dom + '.py ' + out + ' ' + query + ' ' + phyml)

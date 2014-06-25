import sys,os
import Fetchutil, recip_edit, Report

#Expected: 1=query accession no.; 2=out prefix; 3=domain (euk,bac,arch,all)
phyml=''
if '-y' in sys.argv:
	phyml='-y'
try:
    query= sys.argv[1]
    out=sys.argv[2]
    dom=sys.argv[3]
except IndexError:
    query=raw_input('Query accession no: ')
    out=raw_input('Output file: ')
    dom=raw_input('Organism Domain to be explored: euk,bac,arch, or all ').lower()
    phy=raw_input("Should we use PhyML to find trees, in addition to MrBayes? [y/N]'
if phy.lower()=='y':
	phyml='-y'
try:
    open(dom+'-'+out+'.fasta').read()
except IOError:
    arch_list=['Haloferax volcanii','Sulfolobus tokodaii','Methanococcus aeolicus','Methanobrevibacter smithii', 'Thermococcus sibiricus','Archaeoglobus fulgidus','Nanoarchaeum equitans','Thermoplasma acidophilum']
    bac_list= ['Gemmata obscuriglobus', 'Prosthecobacter dejongeii', 'Verrucomicrobium spinosum','Rickettsia prowazekii', 'Agrobacterium tumefaciens','Escherichia coli', 'Bacillus subtilis','Anabaena variabilis', 'Thermotoga maritima']
    euk_list= ['Drosophila melanogaster','Homo sapiens','Oryza sativa', 'Trypanosoma brucei','Plasmodium falciparum','Saccharomyces cerevisiae', 'Neurospora crassa','Arabidopsis thaliana']#subject to change
    #setting threshold values: thresh1-w/ arch ;thresh2-w/ bac; 
    dom_query=Fetchutil.orgfetch(query)[1]
    if dom_query=='Archaea':
        thresh1=1e-10
        thresh2=1e-5
        thresh3=5
    elif dom_query=='Eukaryota':
        thresh1=5
        thresh2=5
        thresh3=1e-10
    else:
        thresh1=1e-5
        thresh2=1e-10
        thresh3=5
    try:
        os.mkdir('XML')
        os.mkdir('Orthos')
    except:
        pass
    arch_accs={}
    bac_accs={}
    euk_accs={}
    print "Blasting"
    if dom=='arch' or dom=='all':
      arch_accs=recip_edit.bestrecipblast(arch_list,query,thresh1)
    if dom=='bac' or dom=='all':
      bac_accs=recip_edit.bestrecipblast(bac_list,query,thresh2)
    if dom=='euk' or dom=='all':
      euk_accs=recip_edit.bestrecipblast(euk_list,query,thresh3)
    all_accs={}
    all_accs.update(arch_accs)
    all_accs.update(bac_accs)
    all_accs.update(euk_accs)
    num_seqs=0
    for j in all_accs:
    	num_seqs+=len(all_accs[j])
    print "Dictionary generated with "+repr(len(all_accs.keys()))+" keys and "+repr(num_seqs)+" sequences."
##Fetching the sequences and writing them to file
    print "Writing seqs to file."
    for a in arch_accs.keys():
        Fetchutil.seqfetch(arch_accs[a][0])
        fil=open('Orthos/'+arch_accs[a][0]+'.fasta')
        fil_arr=fil.readlines()
        fil.close()
        fil=open('Orthos/'+arch_accs[a][0]+'.fasta','w')
        for i in range(len(fil_arr)):
            if i==0:
                fil.write(fil_arr[i].strip()+' '+arch_accs[a][1]+'  '+arch_accs[a][2]+'\n')
            else:
                fil.write(fil_arr[i])
        fil.close()
        Fetchutil.addseq('arch-'+out+'.fas','Orthos/'+arch_accs[a][0]+'.fasta')
	os.remove('Orthos/'+arch_accs[a][0]+'.fasta')
    for b in bac_accs.keys():
        Fetchutil.seqfetch(bac_accs[b][0])
        fil=open('Orthos/'+bac_accs[b][0]+'.fasta')
        fil_arr=fil.readlines()
        fil.close()
        fil=open('Orthos/'+bac_accs[b][0]+'.fasta','w')
        for i in range(len(fil_arr)):
            if i==0:
                fil.write(fil_arr[i].strip()+' '+bac_accs[b][1]+'  '+bac_accs[b][2]+'\n')
            else:
                fil.write(fil_arr[i])
        fil.close()
        Fetchutil.addseq('bac-'+out+'.fas','Orthos/'+bac_accs[b][0]+'.fasta')
        os.remove('Orthos/'+bac_accs[b][0]+'.fasta')
    for e in euk_accs.keys():
        Fetchutil.seqfetch(euk_accs[e][0])
        fil=open('Orthos/'+euk_accs[e][0]+'.fasta')
        fil_arr=fil.readlines()
        fil.close()
        fil=open('Orthos/'+euk_accs[e][0]+'.fasta','w')
        for i in range(len(fil_arr)):
            if i==0:
                fil.write(fil_arr[i].strip()+' '+euk_accs[e][1]+'  '+euk_accs[e][2]+'\n')
            else:
                fil.write(fil_arr[i])
        fil.close()
        Fetchutil.addseq('euk-'+out+'.fas','Orthos/'+euk_accs[e][0]+'.fasta')
        os.remove('Orthos/'+euk_accs[e][0]+'.fasta')

    for c in all_accs.keys():
        Fetchutil.seqfetch(all_accs[c][0])
        fil=open('Orthos/'+all_accs[c][0]+'.fasta')
        fil_arr=fil.readlines()
        fil.close()
        fil=open('Orthos/'+all_accs[c][0]+'.fasta','w')
        for i in range(len(fil_arr)):
            if i==0:
                fil.write(fil_arr[i].strip()+' '+all_accs[c][1]+'  '+all_accs[c][2]+'\n')
            else:
                fil.write(fil_arr[i])
        fil.close()
        Fetchutil.addseq('all-'+out+'.fas','Orthos/'+all_accs[c][0]+'.fasta')
        os.remove('Orthos/'+all_accs[c][0]+'.fasta')

try:
    os.mkdir('aligns')
    os.mkdir('Bayes')
    os.mkdir('ML')
    os.mkdir('Prot')
except:
    pass
os.system('python '+dom+'.py '+out+' '+phyml)
Report.generateReport(out,query,models,dom)

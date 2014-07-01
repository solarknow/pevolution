import Fetchutil, SeqUtil, os, string
if not os.path.exists('Reports'):
  os.mkdir('Reports')
def generateReport(name,quer,models,dom):
    "Generates a report summarizing the analysis done"
    ret=open('Reports/Report-'+dom+'-'+name+'.txt','w')
    ret.write('Orthologous sequence Search and Alignment\n'+
              'Python Scripts Written by Mihir Sarwade\n\n')
    ##List accession numbers, names of genes and their respective organisms
    ret.write('-----------------------------\n')
    ret.write('Reciprocal Best BLAST Results\n')
    ret.write('-----------------------------\n\n')
    ret.write('Query sequence Accension number: '+quer+' Query sequence organism: '+
                Fetchutil.orgfetch(quer)[0]+'\n')
    ret.write('BLASTs are performed using expected value (E-value) thresholds based on the kingdom (Bacteria, Archaea, and Eukaryota).\n'+
    ' These lists are then refined by picking only those sequences that are at least 50% similar and\n'+
    ' whose aligned portion is at least 25% that of the query, as recommended by Moreno-Hagelsieb and Latimer (2008).\n')
    import_file=open('Data/'+dom+'-'+name+'.fas')
    info={}
    while import_file:
        lin=import_file.readline()
        if lin=='':
            break
        spl=lin.split(':')
        
        try:
            spl_data=spl[1].split()
            #print spl_data
            org=spl[0][1:]
            info.update({org:spl_data})
        except IndexError:
            continue
        
    #print info
    for i in info:
        info[i].append(Fetchutil.namefetch(info[i][1]))
    ret.write('Organism                       Kingdom     Accession No.(GI)      No. on list/length of list     No. on list/length of list         ')

    ret.write('Seq Definition\n')
    #print info
    keys=list(info.keys())
    keys.sort()
    length=Fetchutil.findlonglen(info)+6
    for i in keys:
        #print info[i]
    #writing organism
        ret.write(i)
        for j in range(length-len(i)):
            ret.write(' ')
    #writing Kingdom
        ret.write(info[i][0])
        for j in range(20-len(info[i][0])):
            ret.write(' ')
    #writing Accession No.
        ret.write(info[i][1])
        for j in range(25-len(info[i][1])):
            ret.write(' ')
    #writing 1st list no
        ret.write(info[i][2])
        for j in range(20-len(info[i][2])):
            ret.write(' ')
    #writing 2nd list no
        ret.write(info[i][3])
        for j in range(20-len(info[i][3])):
            ret.write(' ')
    #writing Seq Def
        ret.write(info[i][-1]+'\n')
        

    
    ret.write('\n\tThe higher up on the list of accension number garned by the best BLAST protocol, i.e. the smaller the ratio,\n'+
        'the more likely the chosen sequence is an ortholog of the query sequence. Granted that several of the species may have \nseveral copies of the gene in question due to gene duplication events, this function chooses only the one BEST match of the several copies it may encounter.\n\n')
    ##Draw initial alignment +length, final alignment +length
    ret.write('\nAlignments:\n'+
              'Original alignment: Length: ')
    alig=open('aligns/'+dom+'-'+name+'.best.nex')
    length=''
    while alig:
        lin=alig.readline()
        if lin.startswith('dimensions'):
            length=lin.split()[2][6:-1]
            break
    ret.write(length+'\n'+
                'Final alignment: Length: ')
    alig=open('aligns/'+dom+'-'+name+'.best.nex')
    length=''
    print 'alignment printing'
    while alig:
        lin=alig.readline()
        if lin.startswith('dimensions'):
            length=lin.split()[2][6:-1]
            ret.write(length+'\n\n')
        elif lin.startswith('matrix'):
            alig.readline()
            while alig:
                lin=alig.readline()
                if lin==';\n':
                    break
                ret.write(lin)
            break
        elif lin=='end;\n':
            break
    print 'alignment printing done'
    ##How much was taken out of original seqs, i.e. that is left
#    ret.write('\nShowing the original alignment and the final alignment may not give enough information.\n'+
#              'How much of each sequence was conserved in the final alignment:\n')
#    oriseqs={}
#    fil=open(dom+'-'+name)
#    while fil:
#        lin=fil.readline().strip()
#        #print len(oriseqs)
#        if lin=='':
#            break
#        if lin.startswith('>'):
#            title=lin[1:]
#            seq=''
#            while fil:
#                #print 'blah...'
#                lan=fil.readline().strip()
#                if lan=='':
#                    oriseqs.update({title:seq})
#                    break
#                seq+=lan
#    print 'modding'
#    moddedseqs=SeqUtil.modseqs('Bayes/'+dom+'-'+name+'-mod.nxs')
#    print len(oriseqs),len(moddedseqs)
#    for m in moddedseqs:
#        lenmod=len(moddedseqs[m])
#        lenori=len(oriseqs[m])
#        ret.write(m+' : '+str(lenmod)+' residues from the original '+str(lenori)+' residues. About %.2f percent.\n' % ((lenmod+0.0)/lenori*100))
#    print 'all alignment stats written'
    ##Best model(s)+ respective parameters and BIC values(?)
    prot_hand=open('Prot/'+dom+'-'+name+'.pro')
    while prot_hand:
        prot=prot_hand.readline()
        #print prot
        if prot.startswith('Best model'):
            prot_hand.readline()
            prot_hand.readline()
            prot_hand.readline()
            prot_hand.readline()
            ret.write('\nModel          deltaBIC*    BIC          BICw       -lnL    \n'+
                    '------------------------------------------------------------\n')
            while prot_hand:
                prot=prot_hand.readline()
                #print prot
                if float(prot.split()[1])<=200:
                    ret.write(prot)
                else:
                    break
            break
    print 'models printed'              
    
    trees=''
## print trees: ProtTest best tree
#    ret.write('\nTree found by ProtTest using the best model:\n')
#    tree=open('Prot/'+dom+'-'+name+'-mod.pro.tre').readline()
#    trees+=tree
#    ret.write(tree+'\n')
    ##              PhyML1 + (PhyML2)
    try:
      for i in models:
        ret.write('\nTree found by PhyML using the '+i.split('+')[0]+' model:\n')
        tree=consense('ML/'+dom+'-'+name+i.split('+')[0]+'_phyml_boot_trees.txt')
        trees+=tree+'\n'
        ret.write(tree+'\n')
    except:
      ret.write('Tree not found')
    ##              Bayesian selected tree
    ret.write('\nTree found by MrBayes using the best model:\n')
    try:
      read=open('Bayes/'+dom+'-'+name+'-bayes.nxs.con')
    except:
      read=open('Bayes/'+dom+'-'+name+'-bayes.nxs.con.tre')
    taxa={}
    while read:
        lin=read.readline()
        #print lin
        spl=lin.split()
        if lin=='':
            break
        elif len(spl)==0:
            continue
        if spl[0]=='translate':
          spl=read.readline().split()
          #print spl
          while not spl[0]==';':
            taxa.update({spl[0]:spl[1].split(',')[0]})
            spl=read.readline().split()
        print taxa
        if spl[0]=='tree':
          tree_temp=spl[4]
          for i in taxa.keys():
            print i, taxa[i]
            tree_temp=string.replace(tree_temp,i+'[&prob',taxa[i]+'[&prob')
          trees+=tree_temp+'\n'
          ret.write(trees+'\n')
          break
    print 'trees printed'
    ret.close()
    open('Reports/'+dom+'-'+name+'-trees.tre','w').write(trees)

def consense(fil):
    "Returns a newick consensus tree of the trees in fil"
    dum=open('inputer','w')
    dum.write(fil+'\nf\n'+fil.split('.')[0]+'_cons\ny\nf\n'+fil.split('.')[0]+'.tre')
    dum.close()
    os.system('./consense < inputer')
    tree=''
    treefil=open('Reports/'+fil.split('.')[0]+'.tre')
    for i in treefil:
        tree+=i.strip()
    try:
      os.remove('inputer')
    except OSError:
      pass
    return tree

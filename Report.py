import os

import FetchUtil
import SeqUtil

if not os.path.exists('Reports'):
    os.mkdir('Reports')


def generate_report(name, quer, models, dom):
    """Generates a report summarizing the analysis done"""
    with open('Reports' + os.sep + 'Report-' + dom + '-' + name + '.txt', 'w') as ret:
        ret.write('Orthologous sequence Search and Alignment\n' +
                  'Python Scripts Written by Mihir Sarwade\n\n')
        # List accession numbers, names of genes and their respective organisms
        ret.write('-----------------------------\n')
        ret.write('Reciprocal Best BLAST Results\n')
        ret.write('-----------------------------\n\n')
        ret.write('Query sequence Accension number: ' + quer + ' Query sequence organism: ' +
                  FetchUtil.fetch_organism(quer)[0] + '\n')
        ret.write(
            'BLASTs are performed using expected value (E-value) thresholds based on the kingdom (Bacteria, Archaea, '
            'and Eukaryota).\n These lists are then refined by picking only those sequences that are at least '
            '50% similar and\n whose aligned portion is at least 25% that of the query, as recommended by Moreno-Hagelsieb '
            'and Latimer (2008).\n')
        with open('Data' + os.sep + dom + '-' + name + '.fas') as import_file:
            info = {}
            while import_file:
                lin = import_file.readline()
                if lin == '':
                    break
                spl = lin.split(':')

                try:
                    spl_data = spl[1].split()
                    # print spl_data
                    org = spl[0][1:]
                    info.update({org: spl_data})
                except IndexError:
                    continue

        # print info
        for i in info:
            info[i].append(FetchUtil.fetch_definition(info[i][1]))
        ret.write(
            'Organism                       Kingdom     Accession No.(GI)      No. on list/length of list     No. on '
            'list/length of list         ')

        ret.write('Seq Definition\n')
        # print info
        keys = list(info.keys())
        keys.sort()
        length = SeqUtil.find_longest_key_length(info) + 6
        for i in keys:
            # print info[i]
            # writing organism
            ret.write(i)
            for j in range(length - len(i)):
                ret.write(' ')
            # writing Kingdom
            ret.write(info[i][0])
            for j in range(20 - len(info[i][0])):
                ret.write(' ')
            # writing Accession No.
            ret.write(info[i][1])
            for j in range(25 - len(info[i][1])):
                ret.write(' ')
            # writing 1st list no
            ret.write(info[i][2])
            for j in range(20 - len(info[i][2])):
                ret.write(' ')
            # writing 2nd list no
            ret.write(info[i][3])
            for j in range(20 - len(info[i][3])):
                ret.write(' ')
            # writing Seq Def
            ret.write(info[i][-1] + '\n')

        ret.write(
            '\n\tThe higher up on the list of accension number garned by the best BLAST protocol, '
            'i.e. the smaller the ratio,\nthe more likely the chosen sequence is an ortholog of the query sequence. '
            'Granted that several of the species may have \nseveral copies of the gene in question due to gene duplication '
            'events, this function chooses only the one BEST match of the several copies it may encounter.\n\n')
        # Draw initial alignment +length, final alignment +length
        ret.write('\nAlignments:\n' + 'Original alignment: Length: ')
        with open('aligns' + os.sep + dom + '-' + name + '.best.nex') as alig:
            length = ''
            while alig:
                lin = alig.readline()
                if lin.startswith('dimensions'):
                    length = lin.split()[2][6:-1]
                    break
        ret.write(length + '\n' + 'Final alignment: Length: ')
        with open('aligns' + os.sep + dom + '-' + name + '.best.nex') as alig:
            print('alignment printing')
            while alig:
                lin = alig.readline()
                if lin.startswith('dimensions'):
                    length = lin.split()[2][6:-1]
                    ret.write(length + '\n\n')
                elif lin.startswith('matrix'):
                    alig.readline()
                    while alig:
                        lin = alig.readline()
                        if lin == ';\n':
                            break
                        ret.write(lin)
                    break
                elif lin == 'end;\n':
                    break
            print('alignment printing done')

        with open('Prot' + os.sep + dom + '-' + name + '.pro') as prot_hand:
            while prot_hand:
                prot = prot_hand.readline()
                if prot.startswith('Best model'):
                    prot_hand.readline()
                    prot_hand.readline()
                    prot_hand.readline()
                    prot_hand.readline()
                    ret.write('\nModel          deltaBIC*    BIC          BICw       -lnL    \n' +
                              '------------------------------------------------------------\n')
                    while prot_hand:
                        prot = prot_hand.readline()
                        if float(prot.split()[1]) <= 200:
                            ret.write(prot)
                        else:
                            break
                    break
        print('models printed')

        trees = ''
        # print trees: ProtTest best tree
        #    ret.write('\nTree found by ProtTest using the best model:\n')
        #    tree=open('Prot' + os.sep +dom+'-'+name+'-mod.pro.tre').readline()
        #    trees+=tree
        #    ret.write(tree+'\n')
        #              PhyML1 + (PhyML2)
        try:
            for i in models:
                ret.write('\nTree found by PhyML using the ' + i.split('+')[0] + ' model:\n')
                tree = consense('ML' + os.sep + dom + '-' + name + i.split('+')[0] + '_phyml_boot_trees.txt')
                trees += tree + '\n'
                ret.write(tree + '\n')
        except:
            ret.write('Tree not found')
        #              Bayesian selected tree
        ret.write('\nTree found by MrBayes using the best model:\n')
        if os.path.exists('Bayes' + os.sep + dom + '-' + name + '-bayes.nxs.con'):
            read = open('Bayes' + os.sep + dom + '-' + name + '-bayes.nxs.con')
        else:
            read = open('Bayes' + os.sep + dom + '-' + name + '-bayes.nxs.con.tre')
        taxa = {}
        while read:
            lin = read.readline()
            spl = lin.split()
            if lin == '':
                break
            elif len(spl) == 0:
                continue
            if spl[0] == 'translate':
                spl = read.readline().split()
                while not spl[0] == ';':
                    taxa.update({spl[0]: spl[1].split(',')[0]})
                    spl = read.readline().split()
            if spl[0] == 'tree':
                tree_temp = spl[4]
                for i in taxa.keys():
                    tree_temp = tree_temp.replace(i + '[&prob', taxa[i] + '[&prob')
                trees += tree_temp + '\n'
                ret.write(trees + '\n')
                break
        read.close()
        print('trees printed')
    with open('Reports' + os.sep + dom + '-' + name + '-trees.tre', 'w') as tree_file:
        tree_file.write(trees)


def consense(fil):
    """Returns a newick consensus tree of the trees in fil"""
    with open('inputer', 'w') as dum:
        dum.write(fil + '\nf\n' + fil.split('.')[0] + '_cons\ny\nf\n' + fil.split('.')[0] + '.tre')
    os.system('./consense < inputer')
    tree = ''
    with open('Reports' + os.sep + fil.split('.')[0] + '.tre') as treefil:
        for i in treefil:
            tree += i.strip()
    try:
        os.remove('inputer')
    except OSError:
        pass
    return tree

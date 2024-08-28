import os

from Utils import SeqUtil, FetchUtil
from helpers.constants import ReportsPath, DataPath, AlignsPath, ProtPath, MLPath, BayesPath


def generate_report(name, quer, models, dom):
    """Generates a report summarizing the analysis done"""
    with open(str(ReportsPath(f'Report-{dom}-{name}.txt')), 'w') as ret:
        ret.write('Orthologous sequence Search and Alignment\n' +
                  'Python Scripts Written by Mihir Sarwade\n\n')
        # List accession numbers, names of genes and their respective organisms
        ret.write('-----------------------------\n')
        ret.write('Reciprocal Best BLAST Results\n')
        ret.write('-----------------------------\n\n')
        ret.write(f'Query sequence\nAccension number: {quer}\nQuery sequence organism: ' +
                 f'{FetchUtil.fetch_organism(quer)[0]}\n\n')
        ret.write(
            'BLASTs are performed using expected value (E-value) thresholds based on the domain (Bacteria, Archaea, '
            'and Eukaryota).\nThese lists are then refined by picking only those sequences that are at least '
            '50% similar and\nwhose aligned portion is at least 25% that of the query, as recommended by '
            'Moreno-Hagelsieb and Latimer (2008).\n\n')
        with open(str(DataPath(dom + '-' + name + '.fas'))) as import_file:
            info = {}
            while import_file:
                lin = import_file.readline()
                if lin == '':
                    break
                spl = lin.split(':')

                if len(spl) >= 2:
                    spl_data = spl[1].strip().split()
                    # print spl_data
                    org = spl[0]
                    info.update({org: spl_data})

        for i in info:
            info[i].append(FetchUtil.fetch_definition(info[i][1]))
        orgs = list(info.keys())
        orgs.sort()
        max_org_len = SeqUtil.find_longest_key_length(info)[1] + 3
        columns = ['Organism', 'Domain', 'Accession No.(GI)',
                   'No. on list/length of list', 'No. on list/length of list', 'Seq Definition']
        rows = []
        max_lens = [max(len(c), max_org_len) if c == 'Organism' else len(c) for c in columns]
        for o in orgs:
            row = []
            for idx in range(len(columns)):
                if columns[idx] == 'Organism':
                    row.append(o)
                else:
                    max_lens[idx] = max(max_lens[idx], len(str(info[o][idx-1])))
                    row.append(str(info[o][idx-1]))
            rows.append(row)

        padding = 3*' '
        header = padding.join([k+' '*(v-len(k)) for k,v in zip(columns, max_lens)]) + '\n'
        header_spacer = padding.join(['-'*len(k)+' '*(v-len(k)) for k,v in zip(columns, max_lens)]) + '\n'
        ret.write(header)
        ret.write(header_spacer)
        # print info
        for r in rows:
            row = padding.join([k+' '*(v-len(k)) for k,v in zip(r, max_lens)]) + '\n'
            ret.write(row)

        ret.write(
            '\nThe higher up on the list of accession number garnered by the best BLAST protocol, \n'
            'i.e. the smaller the ratio, the more likely the chosen sequence is an ortholog of the query sequence. \n'
            'Granted that several of the species may have several copies of the gene in question due to gene '
            'duplication events, \nthis function chooses only the one BEST match of the several copies it may '
            'encounter.\n\n')
        # Draw initial alignment +length, final alignment +length
        ret.write('\nAlignments:\nOriginal alignment length: ')
        with open(str(AlignsPath(dom + '-' + name + '.best.nex'))) as alig:
            length = ''
            while alig:
                lin = alig.readline()
                if lin.startswith('dimensions'):
                    length = lin.split()[2][6:-1]
                    break
        ret.write(length + '\n' + 'Final alignment length: ')
        with open(str(AlignsPath(dom + '-' + name + '.best.nex'))) as alig:
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

        with open(str(ProtPath(dom + '-' + name + '.pro'))) as prot_hand:
            while prot_hand:
                prot = prot_hand.readline()
                if prot.startswith('Best model'):
                    prot_hand.readline()
                    prot_hand.readline()
                    prot_hand.readline()
                    prot_hand.readline()
                    ret.write('\nModel          deltaBIC*      BIC            BICw           -lnL\n' +
                              '-------------------------------------------------------------------\n')
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
        if os.path.exists(str(MLPath(dom + '-' + name + i.split('+')[0] + '_phyml_boot_trees.txt'))):
            for i in models:
                ret.write('\nTree found by PhyML using the ' + i.split('+')[0] + ' model:\n')
                tree = consense(str(MLPath(dom + '-' + name + i.split('+')[0] + '_phyml_boot_trees.txt')))
                trees += tree + '\n'
                ret.write(tree + '\n')
        #              Bayesian selected tree
        ret.write('\nTree found by MrBayes using the best model:\n')
        if os.path.exists(str(BayesPath(dom + '-' + name + '-bayes.nxs.con'))):
            path = BayesPath(dom + '-' + name + '-bayes.nxs.con')
        else:
            path = BayesPath(dom + '-' + name + '-bayes.nxs.con.tre')
        taxa = {}
        with open(str(path)) as read:
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
                    for i in taxa:
                        tree_temp = tree_temp.replace(i + '[&prob', taxa[i] + '[&prob')
                    trees += tree_temp + '\n'
                    ret.write(trees + '\n')
                    with open(str(ReportsPath(dom + '-' + name + '-trees.tre')), 'w') as tree_file:
                        tree_file.write(trees)
                    break
        print('trees printed')


def consense(fil):
    """Returns a newick consensus tree of the trees in fil"""
    filename = str(fil).split('.')[0]
    with open('inputer', 'w') as dum:
        dum.write(str(fil) + '\nf\n' + filename + '_cons\ny\nf\n' + filename + '.tre')
    os.system('./consense < inputer')
    tree = ''
    with open(str(ReportsPath(filename + '.tre'))) as treefil:
        for i in treefil:
            tree += i.strip()
    return tree

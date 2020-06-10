import AlignUtil
import FetchUtil
import SeqUtil
import os

from constants import ML_PATH, BAYES_PATH, REPORTS_PATH, DATA_PATH, ALIGNS_PATH, PROT_PATH

os.makedirs('Reports', exist_ok=True)


def generateReport(name, protein, models, dom):
    """Generates a report summarizing the analysis done"""
    ret = open(REPORTS_PATH + 'Report-' + dom + '-' + name + '.txt', 'w')
    ret.write('Orthologous sequence Search and Alignment\n' +
              'Python Scripts Written by Mihir Sarwade\n\n')
    # List accession numbers, names of genes and their respective organisms
    ret.write('-----------------------------\n')
    ret.write('Reciprocal Best BLAST Results\n')
    ret.write('-----------------------------\n\n')
    ret.write('Query sequence Accension number: ' + protein.accession + ' Query sequence organism: ' +
              protein.binomial + '\n')
    ret.write(
        'BLASTs are performed using expected value (E-value) thresholds based on the kingdom (Bacteria, Archaea, '
        'and Eukaryota).\n' +
        ' These lists are then refined by picking only those sequences that are at least 50% similar and\n' +
        'whose aligned portion is at least 25% that of the query, as recommended by Moreno-Hagelsieb and Latimer ('
        '2008).\n')
    import_file = open(DATA_PATH + dom + '-' + name + '.fas')
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
        info[i].append(FetchUtil.fetch_protein(info[i][1]).definition)
    ret.write(
        'Organism                       Kingdom     Accession No.(GI)      No. on list/length of list     No. on '
        'list/length of list         ')

    ret.write('Seq Definition\n')
    # print info
    keys = list(info.keys())
    keys.sort()
    length = SeqUtil.longest_key_length(info) + 6
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
        '\n\tThe higher up on the list of accension number garned by the best BLAST protocol, i.e. the smaller the '
        'ratio,\nthe more likely the chosen sequence is an ortholog of the query sequence. Granted that several of '
        'the species may have \nseveral copies of the gene in question due to gene duplication events, this function '
        'chooses only the one BEST match of the several copies it may encounter.\n\n')
    # Draw initial alignment +length, final alignment +length
    ret.write('\nAlignments:\n' +
              'Original alignment: Length: ')
    alig = open(ALIGNS_PATH + dom + '-' + name + '.best.nex')
    length = ''
    while alig:
        lin = alig.readline()
        if lin.startswith('dimensions'):
            length = lin.split()[2][6:-1]
            break
    ret.write(length + '\n' +
              'Final alignment: Length: ')
    alig = open(ALIGNS_PATH + dom + '-' + name + '.best.nex')
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

    # Best model(s)+ respective parameters and BIC values(?)
    prot_hand = open(PROT_PATH + dom + '-' + name + '.pro')
    while prot_hand:
        prot = prot_hand.readline()
        # print prot
        if prot.startswith('Best model'):
            prot_hand.readline()
            prot_hand.readline()
            prot_hand.readline()
            prot_hand.readline()
            ret.write('\nModel          deltaBIC*    BIC          BICw       -lnL    \n' +
                      '------------------------------------------------------------\n')
            while prot_hand:
                prot = prot_hand.readline()
                # print prot
                if float(prot.split()[1]) <= 200:
                    ret.write(prot)
                else:
                    break
            break
    print('models printed')

    trees = ''

    for i in models:
        ret.write('\nTree found by PhyML using the ' + i.split('+')[0] + ' model:\n')
        tree = AlignUtil.consense(ML_PATH + dom + '-' + name + i.split('+')[0] + '_phyml_boot_trees.txt')
        trees += tree + '\n'
        ret.write(tree + '\n')
    #              Bayesian selected tree
    ret.write('\nTree found by MrBayes using the best model:\n')
    if os.path.exists(BAYES_PATH + dom + '-' + name + '-bayes.nxs.con'):
        read = open(BAYES_PATH + dom + '-' + name + '-bayes.nxs.con')
    else:
        read = open(BAYES_PATH + dom + '-' + name + '-bayes.nxs.con.tre')
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
                taxa[spl[0]] = spl[1].split(',')[0]
                spl = read.readline().split()
        if spl[0] == 'tree':
            tree_temp = spl[4]
            for i in taxa:
                tree_temp = tree_temp.replace(i + '[&prob', taxa[i] + '[&prob')
            trees += tree_temp + '\n'
            ret.write(trees + '\n')
            break
    print('trees printed')
    read.close()
    ret.close()
    with open(REPORTS_PATH + dom + '-' + name + '-trees.tre', 'w') as report:
        report.write(trees)


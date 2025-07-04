import os

from Utils import SeqUtil, FetchUtil
from helpers.constants import ReportsPath, DataPath, AlignsPath, ProtPath, MLPath, BayesPath

class Report:
    def __init__(self, name, quer, models, dom):
        self.info = None
        self.name=name
        self.query_accession=quer
        self.best_models=models
        self.domain=dom
        self.domain_name=f'{self.domain}-{self.name}'
        self.report_path=f'Report-{self.domain_name}.txt'

    def generate_header(self):
        return ('Orthologous sequence Search and Alignment\n' 
                'Python Scripts Written by Mihir Sarwade\n\n' 
                '-----------------------------\n' 
                'Reciprocal Best BLAST Results\n' 
                '-----------------------------\n\n' 
                'Query sequence\n'
                f'Accension number: {self.query_accession}: {FetchUtil.fetch_definition(self.query_accession)}\n'
                f'Query sequence organism: '
                f'{FetchUtil.fetch_organism(self.query_accession)[0]}\n\n' 
                'BLASTs are performed using expected value (E-value) thresholds based on the domain (Bacteria, Archaea, ' 
                'and Eukaryota).\nThese lists are then refined by picking only those sequences that are at least ' 
                '50% similar and\nwhose aligned portion is at least 25% that of the query, as recommended by ' 
                'Moreno-Hagelsieb and Latimer (2008).\n\n')

    def read_in_sequence_data(self):
        info = {}
        with open(str(DataPath(self.domain_name + '.fas'))) as import_file:
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
        self.info = info

    def generate_sequence_data(self):
        ret = ''
        orgs = list(self.info.keys())
        orgs.sort()
        max_org_len = SeqUtil.find_longest_key_length(self.info)[1] + 3
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
                    max_lens[idx] = max(max_lens[idx], len(str(self.info[o][idx - 1])))
                    row.append(str(self.info[o][idx - 1]))
            rows.append(row)

        padding = 3 * ' '
        header = padding.join([k + ' ' * (v - len(k)) for k, v in zip(columns, max_lens)]) + '\n'
        header_spacer = padding.join(['-' * len(k) + ' ' * (v - len(k)) for k, v in zip(columns, max_lens)]) + '\n'
        ret += header
        ret += header_spacer
        # print info
        for r in rows:
            row = padding.join([k + ' ' * (v - len(k)) for k, v in zip(r, max_lens)]) + '\n'
            ret += row
        ret += \
            '\nThe higher up on the list of accession number garnered by the best BLAST protocol, \n' +\
            'i.e. the smaller the ratio, the more likely the chosen sequence is an ortholog of the query sequence. \n' +\
            'Granted that several of the species may have several copies of the gene in question due to gene ' +\
            'duplication events, \nthis function chooses only the one BEST match of the several copies it may ' +\
            'encounter.\n\n'
        return ret

    def generate_alignment(self):
        # Draw initial alignment +length, final alignment +length
        ret = ''
        ret += '\nAlignments:\nOriginal alignment length: '
        with open(str(AlignsPath(self.domain_name + '.best.nex'))) as alig:
            length = ''
            while alig:
                lin = alig.readline()
                if lin.startswith('dimensions'):
                    length = lin.split()[2][6:-1]
                    break
        ret += length + '\n' + 'Final alignment length: '
        with open(str(AlignsPath(self.domain_name + '.best.nex'))) as alig:
            print('alignment printing')
            while alig:
                lin = alig.readline()
                if lin.startswith('dimensions'):
                    length = lin.split()[2][6:-1]
                    ret+=length + '\n\n'
                elif lin.startswith('matrix'):
                    alig.readline()
                    while alig:
                        lin = alig.readline()
                        if lin == ';\n':
                            break
                        ret+=lin
                    break
                elif lin == 'end;\n':
                    break
            print('alignment printing done')
        return ret
    def generate_protein_models(self):
        ret = ''
        with open(str(ProtPath(self.domain_name + '.pro'))) as prot_hand:
            while prot_hand:
                prot = prot_hand.readline()
                if prot.startswith('Best model'):
                    prot_hand.readline()
                    prot_hand.readline()
                    prot_hand.readline()
                    prot_hand.readline()
                    ret+= '\nModel          deltaBIC*      BIC            BICw           -lnL\n' +\
                          '-------------------------------------------------------------------\n'
                    while prot_hand:
                        prot = prot_hand.readline()
                        if float(prot.split()[1]) <= 200:
                            ret += prot
                        else:
                            break
                    break

        print('models printed')
        return ret

    def generate_trees(self):
        ret=''
        trees = ''
        for i in self.best_models:
            prefix = f"{self.domain_name}-{i.name.split('+')[0]}-ori"
            if os.path.exists(str(MLPath(prefix + '_phyml_boot_trees.txt'))):
                ret+='\nTree found by PhyML using the ' + i.name.split('+')[0] + ' model:\n'
                tree = self.consense(i)
                trees += tree + '\n'
                ret += tree + '\n'
        #              Bayesian selected tree
        ret+='\nTree found by MrBayes using the best model:\n'
        if os.path.exists(str(BayesPath(self.domain_name + '-bayes.nxs.con'))):
            path = BayesPath(self.domain_name + '-bayes.nxs.con')
        else:
            path = BayesPath(self.domain_name + '-bayes.nxs.con.tre')
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
                    ret+=trees + '\n'
                    with open(str(ReportsPath(self.domain_name + '-trees.tre')), 'w') as tree_file:
                        tree_file.write(trees)
                    break
        print('trees printed')
        return ret

    def consense(self, model):
        """Returns a newick consensus tree of the trees in fil"""
        fil = str(MLPath(f"{self.domain_name}-{model.name.split('+')[0]}-ori_phyml_boot_trees.txt"))
        filename = str(fil).split('.')[0]
        with open('inputer', 'w') as dum:
            dum.write(str(fil) + '\nf\n' + filename + '_cons\ny\nf\n' + filename + '.tre')
        os.system("./consense < inputer")
        tree = ''
        with open(filename + '.tre') as treefil:
            for i in treefil:
                tree += i.strip()
        return tree
    def generate_report(self):
        with open(str(ReportsPath(self.report_path)), 'w') as ret:
            print('Printing header to file')
            ret.write(self.generate_header())
            print('Printed header')
            self.read_in_sequence_data()
            print("Writing out Sequence data")
            ret.write(self.generate_sequence_data())
            print("Done printing out sequence data")
            ret.write(self.generate_alignment())
            ret.write(self.generate_protein_models())
            ret.write(self.generate_trees())


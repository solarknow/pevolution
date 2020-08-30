import os
import subprocess

import Report
import SeqUtil
from constants import DATA_PATH, ALIGNS_PATH, BAYES_PATH, ML_PATH, REPORTS_PATH


def align_sequences_prank(out, dom):
    """
    Aligns found sequences and prepares alignment for tree inference
    :param out: analysis title
    :param dom: domain prefix
    :return: None
    """
    print("Beginning alignment")
    data_prefix = DATA_PATH + dom + '-' + out
    aligns_prefix = ALIGNS_PATH + dom + '-' + out
    SeqUtil.shorten_sequence_names(data_prefix + '.fas')
    if not os.path.exists(aligns_prefix + '.best.nex'):
        subprocess.run(['prank', f'-d={data_prefix}', f'-o={aligns_prefix}',
                        '-f=nexus', '-quiet'])
        SeqUtil.prank_to_mrbayes(aligns_prefix + '.best.nex')
    print("Alignment Complete")


def run_phyml(out, dom, models_ori):
    """
    Runs phyml
    :param models_ori: dict with best protein models
    :param out: analysis title
    :param dom: domain prefix
    :return: None
    """
    for mod in models_ori:
        SeqUtil.nexus_to_proml(ALIGNS_PATH + dom + '-' + out + '.3.nex',
                               ML_PATH + dom + '-' + out + mod.split('+')[0] + '-ori')
        if models_ori[mod][0] == '0' and models_ori[mod][1] == '0':
            subprocess.run(['phyml', '-i', ML_PATH + dom + '-' + out + mod.split('+')[0] + '-ori', '-d', 'aa', '-b',
                            '100', '-m', mod.split('+')[0], '-f', 'e', '-s', 'BEST', '-u',
                            ALIGNS_PATH + dom + '-' + out + '.ed.2.dnd', '-o', 'tl'])
        elif models_ori[mod][0] == '0':
            subprocess.run(['phyml', '-i', ML_PATH + dom + '-' + out + mod.split('+')[0] + '-ori', '-d', 'aa', '-b',
                            '100', '-m', mod.split('+')[0], '-f', 'e', '-v', models_ori[mod][1], '-s', 'BEST', '-u',
                            ALIGNS_PATH + dom + '-' + out + '.2.dnd', '-o', 'tl'])
        elif models_ori[mod][1] == '0':
            subprocess.run(['phyml', '-i', ML_PATH + dom + '-' + out + mod.split('+')[0] + '-ori', '-d', 'aa', '-b',
                            '100', '-m', mod.split('+')[0], '-f', 'e', '-a', models_ori[mod][0], '-s', 'BEST', '-u',
                            ALIGNS_PATH + dom + '-' + out + '.2.dnd', '-o', 'tl'])
        else:
            subprocess.run(['phyml', '-i', ML_PATH + dom + '-' + out + mod.split('+')[0] + '-ori', '-d', 'aa', '-b',
                            '100', '-m', mod.split('+')[0], '-f', 'e', '-v', models_ori[mod][1], '-a',
                            models_ori[mod][0], '-s', 'BEST', '-u', ALIGNS_PATH + dom + '-' + out + '.2.dnd', '-o',
                            'tl'])


def prepare_bayes_files(out, dom, models_ori):
    """
    Finds best protein model and prepares input file for MrBayes tree inference
    :param models_ori: dict with best protein models
    :param out: analysis title
    :param dom: domain prefix
    :return: None
    """
    print('Writing mrbayes input file')
    bayes_prefix = BAYES_PATH + dom + '-' + out
    if not os.path.exists(bayes_prefix + '-bayes.nxs'):
        SeqUtil.write_to_mrbayes(ALIGNS_PATH + dom + '-' + out + '.best.nex', models_ori,
                                 bayes_prefix + '-bayes.nxs')


def run_bayes_and_report(out, dom, protein, models_ori):
    """
    Runs MrBayes tree inference and runs overall report
    :param models_ori: dict with best protein models
    :param out: analysis title
    :param dom: domain prefix
    :param protein: Protein object of original query
    :return: None
    """
    subprocess.run(['mb', BAYES_PATH + dom + '-' + out + '-bayes.nxs'])
    Report.generateReport(out, protein, models_ori, dom)


def consense(filename):
    """Returns a newick consensus tree of the trees in fil"""
    file_prefix = filename.split('.')[0]
    with open('inputer', 'w') as dum:
        dum.write(f'{filename}\nf\n{file_prefix}_cons\ny\nf\n{file_prefix}.tre')
    os.system('./consense < inputer')
    tree = ''
    with open(REPORTS_PATH + file_prefix + '.tre') as treefil:
        for i in treefil:
            tree += i.strip()
    return tree

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
    SeqUtil.shorten_sequence_names(DATA_PATH + dom + '-' + out + '.fas')
    if not os.path.exists(ALIGNS_PATH + dom + '-' + out + '.best.nex'):
        subprocess.run(['prank', '-d=' + DATA_PATH + dom + '-' + out, '-o=' + ALIGNS_PATH + dom + '-' + out,
                        '-f=nexus', '-quiet'])
        SeqUtil.prank_to_mrbayes(ALIGNS_PATH + dom + '-' + out + '.best.nex')
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
    print('Calculating best model for tree inference')
    print('Found, writing mrbayes input file')
    if not os.path.exists(BAYES_PATH + dom + '-' + out + '-bayes.nxs'):
        SeqUtil.write_to_mrbayes(ALIGNS_PATH + dom + '-' + out + '.best.nex', models_ori,
                                 BAYES_PATH + dom + '-' + out + '-bayes.nxs')


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
    Report.generateReport(out, protein, models_ori, 'arch')


def consense(fil):
    """Returns a newick consensus tree of the trees in fil"""
    with open('inputer', 'w') as dum:
        dum.write(fil + '\nf\n' + fil.split('.')[0] + '_cons\ny\nf\n' + fil.split('.')[0] + '.tre')
    os.system('./consense < inputer')
    tree = ''
    with open(REPORTS_PATH + fil.split('.')[0] + '.tre') as treefil:
        for i in treefil:
            tree += i.strip()
    return tree

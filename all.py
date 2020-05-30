import subprocess
import sys

import Report
import SeqUtil
from constants import DATA_PATH, ALIGNS_PATH, BAYES_PATH, ML_PATH

out = sys.argv[1]
query = sys.argv[2]
try:
    paml = sys.argv[3]
    paml = paml == '-y'
except IndexError:
    paml = False
print("Beginning alignment")
SeqUtil.shorten_seqence_names(DATA_PATH + 'all-' + out + '.fas')
subprocess.call(['prank', '-d=' + DATA_PATH + 'all-' + out, '-o=' + ALIGNS_PATH + 'all-' + out, '-f=nexus', '-quiet'])
SeqUtil.prank_to_mrbayes(ALIGNS_PATH + 'all-' + out + '.best.nex')
print("Alignment complete.\nCalculating best model for tree finding")
models_ori = SeqUtil.prottest_best_models(ALIGNS_PATH + 'all-' + out + '.best.nex')
if paml:
    # for mod in models.keys():
    #     SeqUtil.pamlseqnex('Bayes/all-' + out + '-mod.nxs', 'ML/all-' + out + mod.split('+')[0])
    #     if models[mod][0] == '0' and models[mod][1] == '0':
    #         os.system('phyml -i ML/all-' + out + mod.split('+')[0] + ' -d aa -b 100 -m ' + mod.split('+')[0] +
    #                   ' -f e -s BEST -u aligns/all-' + out + '.ed.2.dnd -o tl')
    #     elif models[mod][0] == '0':
    #         os.system('phyml -i ' + 'ML/all-' + out + mod.split('+')[0] + ' -d aa -b 100 -m ' + mod.split('+')[0] +
    #                   ' -f e -v ' + models[mod][1] + ' -s BEST ' +
    #                   '-u aligns/all-' + out + '.ed.2.dnd -o tl')
    #     elif models[mod][1] == '0':
    #         os.system('phyml -i ' + 'ML/all-' + out + mod.split('+')[0] + ' -d aa -b 100 -m ' + mod.split('+')[0] +
    #                   ' -f e -a ' + models[mod][0] + ' -s BEST ' +
    #                   '-u aligns/all-' + out + '.ed.2.dnd -o tl')
    #     else:
    #         os.system('phyml -i ' + 'ML/all-' + out + mod.split('+')[0] + ' -d aa -b 100 -m ' + mod.split('+')[0] +
    #                   ' -f e -v ' + models[mod][1] + ' -a ' + models[mod][0] + ' -s BEST ' +
    #                   '-u aligns/all-' + out + '.ed.2.dnd -o tl')

    for mod in models_ori.keys():
        SeqUtil.nexus_to_proml(ALIGNS_PATH + 'all-' + out + '.3.nex', 'ML/all-' + out + mod.split('+')[0] + '-ori')
        if models_ori[mod][0] == '0' and models_ori[mod][1] == '0':
            subprocess.call(['phyml', '-i', ML_PATH + 'all-' + out + mod.split('+')[0] + '-ori', '-d', 'aa', '-b',
                             '100', '-m', mod.split('+')[0], '-f', 'e', '-s', 'BEST', '-u',
                             ALIGNS_PATH + 'all-' + out + '.ed.2.dnd', '-o', 'tl'])
        elif models_ori[mod][0] == '0':
            subprocess.call(['phyml', '-i', ML_PATH + 'all-' + out + mod.split('+')[0] + '-ori', '-d', 'aa', '-b',
                             '100', '-m', mod.split('+')[0], '-f', 'e', '-v', models_ori[mod][1], '-s', 'BEST', '-u',
                             ALIGNS_PATH + 'all-' + out + '.2.dnd', '-o', 'tl'])
        elif models_ori[mod][1] == '0':
            subprocess.call(['phyml', '-i', ML_PATH + 'all-' + out + mod.split('+')[0] + '-ori', '-d', 'aa', '-b',
                             '100', '-m', mod.split('+')[0], '-f', 'e', '-a', models_ori[mod][0], '-s', 'BEST', '-u',
                             ALIGNS_PATH + 'all-' + out + '.2.dnd', '-o', 'tl'])
        else:
            subprocess.call(['phyml', '-i', ML_PATH + 'all-' + out + mod.split('+')[0] + '-ori', '-d', 'aa', '-b',
                             '100', '-m', mod.split('+')[0], '-f', 'e', '-v', models_ori[mod][1], '-a',
                             models_ori[mod][0], '-s', 'BEST', '-u', ALIGNS_PATH + 'all-' + out + '.2.dnd', '-o', 'tl'])

SeqUtil.write_to_mrbayes('aligns/all-' + out + '.best.nex', models_ori, 'Bayes/all-' + out + '-bayes.nxs')
subprocess.call(['mb', BAYES_PATH + 'all-' + out + '-bayes.nxs'])
Report.generateReport(out, query, models_ori, 'all')

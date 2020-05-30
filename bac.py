import os
import subprocess
import sys

import Report
import SeqUtil
from constants import DATA_PATH, ALIGNS_PATH, BAYES_PATH

out = sys.argv[1]
query = sys.argv[2]
SeqUtil.shorten_seqence_names(DATA_PATH + 'bac-' + out + '.fas')
if not os.path.exists(ALIGNS_PATH + 'bac-' + out + '.best.nex'):
    subprocess.call(['prank', '-d=' + DATA_PATH + 'bac-' + out, '-o=' + ALIGNS_PATH + 'bac-' + out, '-f=nexus',
                     '-quiet'])
    SeqUtil.prank_to_mrbayes(ALIGNS_PATH + 'bac-' + out + '.best.nex')
models_ori = SeqUtil.prottest_best_models(ALIGNS_PATH + 'bac-' + out + '.best.nex')
if not os.path.exists(BAYES_PATH + 'bac-' + out + '-bayes.nxs'):
    SeqUtil.write_to_mrbayes(ALIGNS_PATH + 'bac-' + out + '.best.nex', models_ori,
                             BAYES_PATH + 'bac-' + out + '-bayes.nxs')
subprocess.call(['mb', BAYES_PATH + 'bac-' + out + '-bayes.nxs'])

Report.generateReport(out, query, models_ori, 'bac')

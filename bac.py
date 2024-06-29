import os
import sys

import Report
import SeqUtil

if not os.path.exists('aligns'):
    os.mkdir('aligns')
if not os.path.exists('Bayes'):
    os.mkdir('Bayes')

out = sys.argv[1]
query = sys.argv[2]
SeqUtil.rename_seqs('Data' + os.sep + 'bac-' + out + '.fas')
if not os.path.exists('aligns' + os.sep + 'bac-' + out + '.best.nex'):
    os.system('prank -d=Data' + os.sep + 'bac-' + out + ' -o=aligns' + os.sep + 'bac-' + out + ' -f=nexus -quiet')
    SeqUtil.bayes_in_nex('aligns' + os.sep + 'bac-' + out + '.best.nex')
models_ori = SeqUtil.bestmod('aligns' + os.sep + 'bac-' + out + '.best.nex')
if not os.path.exists('Bayes' + os.sep + 'bac-' + out + '-bayes.nxs'):
    SeqUtil.bayesfile('aligns' + os.sep + 'bac-' + out + '.best.nex', models_ori,
                      'Bayes' + os.sep + 'bac-' + out + '-bayes.nxs')
os.system('mb Bayes' + os.sep + 'bac-' + out + '-bayes.nxs')
Report.generate_report(out, query, models_ori, 'bac')

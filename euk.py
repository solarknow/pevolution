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
SeqUtil.rename_seqs('Data/euk-' + out + '.fas')
if not os.path.exists('aligns/euk-' + out + '.best.nex'):
    os.system('prank -d=Data/euk-' + out + ' -o=aligns/euk-' + out + ' -f=nexus -quiet')
    SeqUtil.bayes_in_nex('aligns/euk-' + out + '.best.nex')
models_ori = SeqUtil.bestmod('aligns/euk-' + out + '.best.nex')
if not os.path.exists('Bayes/euk-' + out + '-bayes.nxs'):
    SeqUtil.bayesfile('aligns/euk-' + out + '.best.nex', models_ori, 'Bayes/euk-' + out + '-bayes.nxs')
os.system('mb Bayes/euk-' + out + '-bayes.nxs')
Report.generate_report(out, query, models_ori, 'euk')

import os, sys
import SeqUtil, Report

if not os.path.exists('aligns'):
    os.mkdir('aligns')
if not os.path.exists('Bayes'):
    os.mkdir('Bayes')

out = sys.argv[1]
query = sys.argv[2]
SeqUtil.rename_seqs('Data/bac-' + out + '.fas')
if not os.path.exists('aligns/bac-' + out + '.best.nex'):
    os.system('prank -d=Data/bac-' + out + ' -o=aligns/bac-' + out + ' -f=nexus -quiet')
    SeqUtil.bayes_in_nex('aligns/bac-' + out + '.best.nex')
models_ori = SeqUtil.bestmod('aligns/bac-' + out + '.best.nex')
if not os.path.exists('Bayes/bac-' + out + '-bayes.nxs'):
    SeqUtil.bayesfile('aligns/bac-' + out + '.best.nex', models_ori, 'Bayes/bac-' + out + '-bayes.nxs')
os.system('mb Bayes/bac-' + out + '-bayes.nxs')
Report.generateReport(out, query, models_ori, 'bac')

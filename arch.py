import os, sys
import SeqUtil, Report

if not os.path.exists('aligns'):
    os.mkdir('aligns')
if not os.path.exists('Bayes'):
    os.mkdir('Bayes')

out = sys.argv[1]
query = sys.argv[2]
SeqUtil.rename_seqs('Data/arch-' + out + '.fas')
if not os.path.exists('aligns/arch-' + out + '.best.nex'):
    os.system('prank -d=Data/arch-' + out + ' -o=aligns/arch-' + out + ' -f=nexus -quiet')
    SeqUtil.bayes_in_nex('aligns/arch-' + out + '.best.nex')
models_ori = SeqUtil.bestmod('aligns/arch-' + out + '.best.nex')
if not os.path.exists('Bayes/arch-' + out + '-bayes.nxs'):
    SeqUtil.bayesfile('aligns/arch-' + out + '.best.nex', models_ori, 'Bayes/arch-' + out + '-bayes.nxs')
os.system('mb Bayes/arch-' + out + '-bayes.nxs')
Report.generateReport(out, query, models_ori, 'arch')

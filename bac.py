import os
import sys

import Report
import SeqUtil

if not os.path.exists('aligns'):
    os.mkdir('aligns')
if not os.path.exists('Bayes'):
    os.mkdir('Bayes')
# if not os.path.exists('ML'):
#  os.mkdir('ML')

out = sys.argv[1]
query = sys.argv[2]
SeqUtil.rename('Data/bac-' + out + '.fas')
if not os.path.exists('aligns/bac-' + out + '.best.nex'):
    os.system('prank -d=Data/bac-' + out + ' -o=aligns/bac-' + out + ' -f=nexus -quiet')
    SeqUtil.bayesinNex('aligns/bac-' + out + '.best.nex')
# SeqUtil.splicealign('aligns/bac-'+out+'.best.nex','Bayes/bac-'+out+'-mod.nxs')
# models=SeqUtil.bestmod('Bayes/bac-'+out+'-mod.nxs')
models_ori = SeqUtil.bestmod('aligns/bac-' + out + '.best.nex')
if not os.path.exists('Bayes/bac-' + out + '-bayes.nxs'):
    SeqUtil.bayesfile('aligns/bac-' + out + '.best.nex', models_ori, 'Bayes/bac-' + out + '-bayes.nxs')
# SeqUtil.bayesfile('Bayes/bac-'+out+'-mod.nxs',models,'Bayes/bac-'+out+'-bayes.nxs')
os.system('mb Bayes/bac-' + out + '-bayes.nxs')
# SeqUtil.pamlseqnex('Bayes/bac-'+out+'-mod.nxs','ML/bac-'+out)
# for mod in models.keys():
#    SeqUtil.pamlinput('ML/bac-'+out,'ML/bac-'+out+'.out','ML/bac-'+out+'.ctl',{models.keys()[mod].split('+')[0]:models[models.keys()[mod]][1]})
#    os.system('codeml ML/bac-'+out+'.ctl')
#    SeqUtil.extractMLtree('ML/bac-'+out+'.out')
Report.generateReport(out, query, models_ori, 'bac')

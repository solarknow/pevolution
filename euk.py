import os, sys
import SeqUtil, Report

if not os.path.exists('aligns'):
  os.mkdir('aligns')
if not os.path.exists('Bayes'):
  os.mkdir('Bayes')
#if not os.path.exists('ML'):
#  os.mkdir('ML')
  
out= sys.argv[1]
query = sys.argv[2]
SeqUtil.rename('Data/euk-'+out+'.fas')
if not os.path.exists('aligns/euk-'+out+'.best.nex'):
  os.system('prank/bin/prank -d=Data/euk-'+out+' -o=aligns/euk-'+out+' -f=nexus -quiet')
  SeqUtil.bayesinNex('aligns/euk-'+out+'.best.nex')
#SeqUtil.splicealign('aligns/euk-'+out+'.best.nex','Bayes/euk-'+out+'-mod.nxs')
#models=SeqUtil.bestmod('Bayes/euk-'+out+'-mod.nxs')
models_ori=SeqUtil.bestmod('aligns/euk-'+out+'.best.nex')
if not os.path.exists('Bayes/euk-'+out+'-bayes.nxs'):
  SeqUtil.bayesfile('aligns/euk-'+out+'.best.nex',models_ori,'Bayes/euk-'+out+'-bayes.nxs')
#SeqUtil.bayesfile('Bayes/euk-'+out+'-mod.nxs',models,'Bayes/euk-'+out+'-bayes.nxs')
os.system('mb Bayes/euk-'+out+'-bayes.nxs')
#SeqUtil.pamlseqnex('Bayes/euk-'+out+'-mod.nxs','ML/euk-'+out)
#for mod in models.keys():
#    SeqUtil.pamlinput('ML/euk-'+out,'ML/euk-'+out+'.out','ML/euk-'+out+'.ctl',{models.keys()[mod].split('+')[0]:models[models.keys()[mod]][1]})
#    os.system('codeml ML/euk-'+out+'.ctl')
#    SeqUtil.extractMLtree('ML/euk-'+out+'.out')
Report.generateReport(out,query,models_ori,'euk')

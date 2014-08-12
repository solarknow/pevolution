import os, sys
import SeqUtil, Report

if not os.path.exists('aligns'):
  os.mkdir('aligns')
if not os.path.exists('Bayes'):
  os.mkdir('Bayes')
  
out= sys.argv[1]
query = sys.argv[2]
try:
  paml=sys.argv[3]
  paml= paml=='-y'
except IndexError:
  paml=False

SeqUtil.rename('Data/arch-'+out+'.fas')
if not os.path.exists('aligns/arch-'+out+'.best.nex'):
  os.system('prank -d=Data/arch-'+out+' -o=aligns/arch-'+out+' -f=nexus -quiet')
  SeqUtil.bayesinNex('aligns/arch-'+out+'.best.nex')
#SeqUtil.splicealign('aligns/arch-'+out+'.best.nex','Bayes/arch-'+out+'-mod.nxs')
#models=SeqUtil.bestmod('Bayes/arch-'+out+'-mod.nxs')
models_ori=SeqUtil.bestmod('aligns/arch-'+out+'.best.nex')
if not os.path.exists('Bayes/arch-'+out+'-bayes.nxs'):
  SeqUtil.bayesfile('aligns/arch-'+out+'.best.nex',models_ori,'Bayes/arch-'+out+'-bayes.nxs')
#SeqUtil.bayesfile('Bayes/arch-'+out+'-mod.nxs',models,'Bayes/arch-'+out+'-bayes.nxs')
os.system('mb Bayes/arch-'+out+'-bayes.nxs')
if paml:
  if not os.path.exists('ML'):
    os.mkdir('ML')
  SeqUtil.pamlseqnex('Bayes/arch-'+out+'-mod.nxs','ML/arch-'+out)
  for mod in models.keys():
    SeqUtil.pamlinput('ML/arch-'+out,'ML/arch-'+out+'.out','ML/arch-'+out+'.ctl',{models.keys()[mod].split('+')[0]:models[models.keys()[mod]][1]})
    os.system('codeml ML/arch-'+out+'.ctl')
    SeqUtil.extractMLtree('ML/arch-'+out+'.out')
Report.generateReport(out,query,models_ori,'arch',paml)

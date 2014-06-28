import os, sys
import SeqUtil, Report

if not os.path.exists('aligns'):
  os.mkdir('aligns')
if not os.path.exists('Bayes'):
  os.mkdir('Bayes')
#if not os.path.exists('ML'):
#  os.mkdir('ML')

out= sys.argv[1]
query=sys.argv[2]
paml=sys.argv[3]
paml= paml=='-y'
SeqUtil.rename('Data/arch-'+out+'.fas')
os.system('prank/bin/prank -d=Data/arch-'+out+' -o=aligns/arch-'+out+' -f=nexus -quiet')
SeqUtil.splicealign('aligns/arch-'+out+'.2.nex','Bayes/arch-'+out+'-mod.nxs')
models=SeqUtil.bestmod('Bayes/arch-'+out+'-mod.nxs')
models_ori=SeqUtil.bestmod('aligns/arch-'+out+'.2.nex')

if paml:
  for mod in models.keys():
     SeqUtil.pamlseqnex('Bayes/arch-'+out+'-mod.nxs','ML/arch-'+out)
     SeqUtil.pamlinput('ML/arch-'+out,'ML/arch-'+out+'.out','ML/arch-'+out+'.ctl',{models.keys()[mod].split('+')[0]:models[models.keys()[mod]][1]})
     os.system('codeml ML/arch-'+out+'.ctl')
     SeqUtil.extractMLtree('ML/arch-'+out+'.out')
  SeqUtil.pamlseqnex('Bayes/arch-'+out+'-mod.nxs','ML/arch-'+out)
  for mod in models.keys():
     SeqUtil.pamlinput('ML/arch-'+out,'ML/arch-'+out+'.out','ML/arch-'+out+'.ctl',{models.keys()[mod].split('+')[0]:models[models.keys()[mod]][1]})
     os.system('codeml ML/arch-'+out+'.ctl')
     SeqUtil.extractMLtree('ML/arch-'+out+'.out')

SeqUtil.bayesfile('Bayes/arch-'+out+'-mod.nxs',models,'Bayes/arch-'+out+'-bayes.nxs')
os.system('mb Bayes/arch-'+out+'-bayes.nxs')
Report.generateReport(out,query,models_ori,'arch')

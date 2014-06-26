import os, sys
import SeqUtil, Report
out= sys.argv[1]
query = sys.argv[2]
SeqUtil.rename('euk-'+out+'.fas')
os.system('prank/bin/prank -d=euk-'+out+' -o=aligns/euk-'+out+' -f=nexus -quiet')
SeqUtil.splicealign('aligns/euk-'+out+'.best.nex','Bayes/euk-'+out+'-mod.nxs')
#models=SeqUtil.bestmod('Bayes/euk-'+out+'-mod.nxs')
models_ori=SeqUtil.bestmod('aligns/euk-'+out+'.best.nex')

SeqUtil.bayesfile('Bayes/euk-'+out+'-mod.nxs',models,'Bayes/euk-'+out+'-bayes.nxs')
os.system('mb Bayes/euk-'+out+'-bayes.nxs')
#SeqUtil.pamlseqnex('Bayes/euk-'+out+'-mod.nxs','ML/euk-'+out)
#for mod in models.keys():
#    SeqUtil.pamlinput('ML/euk-'+out,'ML/euk-'+out+'.out','ML/euk-'+out+'.ctl',{models.keys()[mod].split('+')[0]:models[models.keys()[mod]][1]})
#    os.system('codeml ML/euk-'+out+'.ctl')
#    SeqUtil.extractMLtree('ML/euk-'+out+'.out')
Report.generateReport(out,query,models_ori,'euk')

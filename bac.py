import os, sys
import SeqUtil
out= sys.argv[1]
SeqUtil.rename('bac-'+out+'.fas')
os.system('clustalw -align -infile=bac-'+out+' -outfile=aligns/bac-'+out+' -output=nexus -quiet')
SeqUtil.splicealign('aligns/bac-'+out,'Bayes/bac-'+out+'-mod.nxs')
models=SeqUtil.bestmod('Bayes/bac-'+out+'-mod.nxs')
models_ori=SeqUtil.bestmod('aligns/bac-'+out)

SeqUtil.bayesfile('Bayes/bac-'+out+'-mod.nxs',models,'Bayes/bac-'+out+'-bayes.nxs')
os.system('mb Bayes/bac-'+out+'-bayes.nxs')
#SeqUtil.pamlseqnex('Bayes/bac-'+out+'-mod.nxs','ML/bac-'+out)
#for mod in models.keys():
#    SeqUtil.pamlinput('ML/bac-'+out,'ML/bac-'+out+'.out','ML/bac-'+out+'.ctl',{models.keys()[mod].split('+')[0]:models[models.keys()[mod]][1]})
#    os.system('codeml ML/bac-'+out+'.ctl')
#    SeqUtil.extractMLtree('ML/bac-'+out+'.out')

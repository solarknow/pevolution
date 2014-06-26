import os, sys
import SeqUtil, Report

out=sys.argv[1]
query=sys.argv[2]
paml=sys.argv[3]
paml= paml=='-y'
print "Beginning alignment"
SeqUtil.rename('all-'+out+'.fas')
os.system('prank/bin/prank -d=all-'+out+' -o=aligns/all-'+out+' -f=nexus -quiet')
SeqUtil.splicealign('aligns/all-'+out+'.2.nex','Bayes/all-'+out+'-mod.nxs')
print "Alignment complete.\nCalculating best model for tree finding"
models_ori=SeqUtil.bestmod('aligns/all-'+out+'.2.nex')
models=SeqUtil.bestmod('Bayes/all-'+out+'-mod.nxs')
print models_ori, models
if paml:
  for mod in models.keys():
     SeqUtil.pamlseqnex('Bayes/all-'+out+'-mod.nxs','ML/all-'+out+mod.split('+')[0])
     if models[mod][0]=='0' and models[mod][1]=='0':
         os.system('phyml -i ml/all-'+out+mod.split('+')[0]+' -d aa -b 100 -m '+mod.split('+')[0]+
                   ' -f e -s BEST -u aligns/all-'+out+'.ed.2.dnd -o tl')
     elif models[mod][0]=='0':
         os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+' -d aa -b 100 -m '+mod.split('+')[0]+
                   ' -f e -v '+models[mod][1]+' -s BEST '+
                   '-u aligns/all-'+out+'.ed.2.dnd -o tl')
     elif models[mod][1]=='0':
         os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+' -d aa -b 100 -m '+mod.split('+')[0]+
                   ' -f e -a '+models[mod][0]+' -s BEST '+
                   '-u aligns/all-'+out+'.ed.2.dnd -o tl')
     else:
         os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+' -d aa -b 100 -m '+mod.split('+')[0]+
                   ' -f e -v '+models[mod][1]+' -a '+models[mod][0]+' -s BEST '+
                   '-u aligns/all-'+out+'.ed.2.dnd -o tl')

  for mod in models_ori.keys():
    SeqUtil.pamlseqnex('aligns/all-'+out+'.3.nex','ML/all-'+out+mod.split('+')[0]+'-ori')
     if models_ori[mod][0]=='0' and models_ori[mod][1]=='0':
         os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+'-ori -d aa -b 100 -m '+mod.split('+')[0]+
                   ' -f e -s BEST -u aligns/all-'+out+'.ed.2.dnd -o tl')
     elif models_ori[mod][0]=='0':
         os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+'-ori -d aa -b 100 -m '+mod.split('+')[0]+
                   ' -f e -v '+models_ori[mod][1]+' -s BEST '+
                   '-u aligns/all-'+out+'.2.dnd -o tl')
     elif models_ori[mod][1]=='0':
         os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+'-ori -d aa -b 100 -m '+mod.split('+')[0]+
                   ' -f e -a '+models_ori[mod][0]+' -s BEST '+
                   '-u aligns/all-'+out+'.2.dnd -o tl')
     else:
         os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+'-ori -d aa -b 100 -m '+mod.split('+')[0]+
                   ' -f e -v '+models_ori[mod][1]+' -a '+models_ori[mod][0]+' -s BEST '+
                   '-u aligns/all-'+out+'.2.dnd -o tl')
        
SeqUtil.bayesfile('Bayes/all-'+out+'-mod.nxs',models,'Bayes/all-'+out+'-bayes.nxs')
os.system('mb Bayes/all-'+out+'-bayes.nxs')
SeqUtil.bayesfile('aligns/all-'+out+'.2.nex',models_ori,'Bayes/all-'+out+'-ori-bayes.nxs')
os.system('mb Bayes/all-'+out+'-ori-bayes.nxs')
Report.generateReport(out,query,models_ori,'all')


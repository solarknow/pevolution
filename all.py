import os, sys
import SeqUtil

out=sys.argv[1]

SeqUtil.rename('all-'+out+'.fas')
os.system('./prank -d=all-'+out+' -o=aligns/all-'+out+' -f=nexus -quiet')
SeqUtil.splicealign('aligns/all-'+out+'.2.nex','Bayes/all-'+out+'-mod.nxs')
os.system('java -jar prottest3.3.jar -i aligns/all-'+out+'.3.nex -o Prot/all-'+out+
            '-ori.pro -all-distributions -all T -S 1 -threads 6 -BIC')
os.system('java -jar prottest3.3.jar -i Bayes/all-'+out+'-mod.nxs -o Prot/all-'+out+
            '-mod.pro -all-distributions -all T -S 1 -threads 6 -BIC')
models_ori=SeqUtil.bestmod('aligns/all-'+out+'.2.nex')
models=SeqUtil.bestmod('Bayes/all-'+out+'-mod.nxs')
print models_ori, models

#for mod in models.keys():
#    SeqUtil.pamlseqnex('Bayes/all-'+out+'-mod.nxs','ML/all-'+out+mod.split('+')[0])
#    if models[mod][0]=='0' and models[mod][1]=='0':
#        os.system('phyml -i ml/all-'+out+mod.split('+')[0]+' -d aa -b 100 -m '+mod.split('+')[0]+
#                  ' -f e -s BEST -u aligns/all-'+out+'.ed.2.dnd -o tl')
#    elif models[mod][0]=='0':
#        os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+' -d aa -b 100 -m '+mod.split('+')[0]+
#                  ' -f e -v '+models[mod][1]+' -s BEST '+
#                  '-u aligns/all-'+out+'.ed.2.dnd -o tl')
#    elif models[mod][1]=='0':
#        os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+' -d aa -b 100 -m '+mod.split('+')[0]+
#                  ' -f e -a '+models[mod][0]+' -s BEST '+
#                  '-u aligns/all-'+out+'.ed.2.dnd -o tl')
#    else:
#        os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+' -d aa -b 100 -m '+mod.split('+')[0]+
#                  ' -f e -v '+models[mod][1]+' -a '+models[mod][0]+' -s BEST '+
#                  '-u aligns/all-'+out+'.ed.2.dnd -o tl')

#for mod in models_ori.keys():
#   SeqUtil.pamlseqnex('aligns/all-'+out+'.3.nex','ML/all-'+out+mod.split('+')[0]+'-ori')
#    if models_ori[mod][0]=='0' and models_ori[mod][1]=='0':
#        os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+'-ori -d aa -b 100 -m '+mod.split('+')[0]+
#                  ' -f e -s BEST -u aligns/all-'+out+'.ed.2.dnd -o tl')
#    elif models_ori[mod][0]=='0':
#        os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+'-ori -d aa -b 100 -m '+mod.split('+')[0]+
#                  ' -f e -v '+models_ori[mod][1]+' -s BEST '+
#                  '-u aligns/all-'+out+'.2.dnd -o tl')
#    elif models_ori[mod][1]=='0':
#        os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+'-ori -d aa -b 100 -m '+mod.split('+')[0]+
#                  ' -f e -a '+models_ori[mod][0]+' -s BEST '+
#                  '-u aligns/all-'+out+'.2.dnd -o tl')
#    else:
#        os.system('phyml -i '+'ml/all-'+out+mod.split('+')[0]+'-ori -d aa -b 100 -m '+mod.split('+')[0]+
#                  ' -f e -v '+models_ori[mod][1]+' -a '+models_ori[mod][0]+' -s BEST '+
#                  '-u aligns/all-'+out+'.2.dnd -o tl')
        
SeqUtil.bayesfile('Bayes/all-'+out+'-mod.nxs',models,'Bayes/all-'+out+'-bayes.nxs')
os.system('mb Bayes/all-'+out+'-bayes.nxs')
SeqUtil.bayesfile('aligns/all-'+out+'.2.nex',models_ori,'Bayes/all-'+out+'-ori-bayes.nxs')
os.system('mb Bayes/all-'+out+'-ori-bayes.nxs')


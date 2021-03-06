import os
import sys

import Report
import SeqUtil

if not os.path.exists('aligns'):
    os.mkdir('aligns')
if not os.path.exists('Bayes'):
    os.mkdir('Bayes')

out = sys.argv[1]
query = sys.argv[2]
try:
    paml = sys.argv[3]
    paml = paml == '-y'
except IndexError:
    paml = False
print("Beginning alignment")
SeqUtil.rename_seqs('Data/all-' + out + '.fas')
os.system('prank -d=Data/all-' + out + ' -o=aligns/all-' + out + ' -f=nexus -quiet')
SeqUtil.bayes_in_nex('aligns/all-' + out + '.best.nex')
print("Alignment complete.\nCalculating best model for tree finding")
models_ori = SeqUtil.bestmod('aligns/all-' + out + '.best.nex')
models = SeqUtil.bestmod('Bayes/all-'+out+'-mod.nxs')
if paml:
    for mod in models.keys():
        SeqUtil.nexus_to_proml('Bayes/all-' + out + '-mod.nxs', 'ML/all-' + out + mod.split('+')[0])
        if models[mod][0] == '0' and models[mod][1] == '0':
            os.system('phyml -i ML/all-' + out + mod.split('+')[0] + ' -d aa -b 100 -m ' + mod.split('+')[0] +
                      ' -f e -s BEST -u aligns/all-' + out + '.ed.2.dnd -o tl')
        elif models[mod][0] == '0':
            os.system('phyml -i ' + 'ML/all-' + out + mod.split('+')[0] + ' -d aa -b 100 -m ' + mod.split('+')[0] +
                      ' -f e -v ' + models[mod][1] + ' -s BEST ' +
                      '-u aligns/all-' + out + '.ed.2.dnd -o tl')
        elif models[mod][1] == '0':
            os.system('phyml -i ' + 'ML/all-' + out + mod.split('+')[0] + ' -d aa -b 100 -m ' + mod.split('+')[0] +
                      ' -f e -a ' + models[mod][0] + ' -s BEST ' +
                      '-u aligns/all-' + out + '.ed.2.dnd -o tl')
        else:
            os.system('phyml -i ' + 'ML/all-' + out + mod.split('+')[0] + ' -d aa -b 100 -m ' + mod.split('+')[0] +
                      ' -f e -v ' + models[mod][1] + ' -a ' + models[mod][0] + ' -s BEST ' +
                      '-u aligns/all-' + out + '.ed.2.dnd -o tl')

    for mod in models_ori.keys():
        SeqUtil.nexus_to_proml('aligns/all-' + out + '.3.nex', 'ML/all-' + out + mod.split('+')[0] + '-ori')
        if models_ori[mod][0] == '0' and models_ori[mod][1] == '0':
            os.system('phyml -i ' + 'ML/all-' + out + mod.split('+')[0] + '-ori -d aa -b 100 -m ' + mod.split('+')[0] +
                      ' -f e -s BEST -u aligns/all-' + out + '.ed.2.dnd -o tl')
        elif models_ori[mod][0] == '0':
            os.system('phyml -i ' + 'ML/all-' + out + mod.split('+')[0] + '-ori -d aa -b 100 -m ' + mod.split('+')[0] +
                      ' -f e -v ' + models_ori[mod][1] + ' -s BEST ' +
                      '-u aligns/all-' + out + '.2.dnd -o tl')
        elif models_ori[mod][1] == '0':
            os.system('phyml -i ' + 'ML/all-' + out + mod.split('+')[0] + '-ori -d aa -b 100 -m ' + mod.split('+')[0] +
                      ' -f e -a ' + models_ori[mod][0] + ' -s BEST ' +
                      '-u aligns/all-' + out + '.2.dnd -o tl')
        else:
            os.system('phyml -i ' + 'ML/all-' + out + mod.split('+')[0] + '-ori -d aa -b 100 -m ' + mod.split('+')[0] +
                      ' -f e -v ' + models_ori[mod][1] + ' -a ' + models_ori[mod][0] + ' -s BEST ' +
                      '-u aligns/all-' + out + '.2.dnd -o tl')

SeqUtil.bayesfile('aligns/all-' + out + '.best.nex', models_ori, 'Bayes/all-' + out + '-bayes.nxs')
os.system('mb Bayes/all-' + out + '-bayes.nxs')
Report.generate_report(out, query, models_ori, 'all')

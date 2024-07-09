import os
import sys

import Report
import SeqUtil
from helpers.commands import run_prank, run_bayes, run_phyml
from helpers.constants import DataPath, AlignsPath, BayesPath, MLPath

out = sys.argv[1]
query = sys.argv[2]
dom = sys.argv[3]
try:
    paml = sys.argv[4]
    paml = paml == '-y'
except IndexError:
    paml = False
print("Beginning alignment")

SeqUtil.rename_seqs(DataPath(f'{dom}-{out}.fas'))
if not os.path.exists(str(AlignsPath(f'{dom}-{out}.best.nex'))):
    run_prank(str(DataPath(f'{dom}-{out}')), str(AlignsPath(f'{dom}-{out}')))
    SeqUtil.bayes_in_nex(AlignsPath(f'{dom}-{out}.best.nex'))
print("Alignment complete.\nCalculating best model for tree finding")
models_ori = SeqUtil.best_model(AlignsPath(f'{dom}-{out}.best.nex'))
models = SeqUtil.best_model(BayesPath(f'{dom}-{out}-mod.nxs'))
if paml:
    for mod, params in models.items():
        SeqUtil.nexus_to_proml(BayesPath(f'{dom}-{out}-mod.nxs'), MLPath(f'{dom}-' + out + mod.split('+')[0]))
        if params[0] == '0' and params[1] == '0':
            run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
                      AlignsPath(f'{dom}-{out}.ed.2.dnd'))
        elif params[0] == '0':
            run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
                      AlignsPath(f'{dom}-{out}.ed.2.dnd'), v=params[1])
        elif params[1] == '0':
            run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
                      AlignsPath(f'{dom}-{out}.ed.2.dnd'), a=params[0])
        else:
            run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
                      AlignsPath(f'{dom}-{out}.ed.2.dnd'), a=params[0], v=params[1])

    for mod, params in models_ori.items():
        SeqUtil.nexus_to_proml(AlignsPath(f'{dom}-{out}.3.nex'), MLPath(f"{dom}-{out}{mod.split('+')[0]}-ori"))
        if params[0] == '0' and params[1] == '0':
            run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
                      AlignsPath(f'{dom}-{out}.2.dnd'), ori=True)
        elif params[0] == '0':
            run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
                      AlignsPath(f'{dom}-{out}.2.dnd'), v=params[1], ori=True)
        elif params[1] == '0':
            run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
                      AlignsPath(f'{dom}-{out}.2.dnd'), a=params[0], ori=True)
        else:
            run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
                      AlignsPath(f'{dom}-{out}.2.dnd'), a=params[0], v=params[1], ori=True)

if not os.path.exists(str(BayesPath(f'{dom}-{out}-bayes.nxs'))):
    SeqUtil.bayesfile(AlignsPath(f'{dom}-{out}.best.nex'), models_ori,
                      BayesPath(f'{dom}-{out}-bayes.nxs'))
run_bayes(BayesPath(f'{dom}-{out}-bayes.nxs'))
Report.generate_report(out, query, models_ori, dom)

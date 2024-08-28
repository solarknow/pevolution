import os
import sys

import Report
from Utils import SeqUtil
from helpers.commands import run_prank, run_bayes, run_phyml, run_prottest
from helpers.constants import DataPath, AlignsPath, BayesPath, MLPath, ProtPath

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
    print('Aligning sequences')
    run_prank(str(DataPath(f'{dom}-{out}')), str(AlignsPath(f'{dom}-{out}')))
    SeqUtil.bayes_in_nex(AlignsPath(f'{dom}-{out}.best.nex'))
print("Alignment complete.\nCalculating best model for tree finding")
if not os.path.exists(str(ProtPath(f'{dom}-{out}.pro'))):
    run_prottest(str(AlignsPath(f'{dom}-{out}.best.nex')), str(ProtPath(f'{dom}-{out}.pro')))
original_models = SeqUtil.best_model(str(ProtPath(f'{dom}-{out}.pro')))
print(original_models)
# if paml:
#     print('pamling')
#     for mod, params in original_models.items():
#         print(mod, params)
#         SeqUtil.nexus_to_proml(AlignsPath(f'{dom}-{out}.best.nex'), MLPath(f"{dom}-{out}{mod.split('+')[0]}-ori"))
#         if params[0] == '0' and params[1] == '0':
#             run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
#                       AlignsPath(f'{dom}-{out}.best.dnd'), ori=True)
#         elif params[0] == '0':
#             run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
#                       AlignsPath(f'{dom}-{out}.2.dnd'), v=params[1], ori=True)
#         elif params[1] == '0':
#             run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
#                       AlignsPath(f'{dom}-{out}.2.dnd'), a=params[0], ori=True)
#         else:
#             run_phyml(MLPath(f'{dom}-{out}' + mod.split('+')[0]), mod.split('+')[0],
#                       AlignsPath(f'{dom}-{out}.2.dnd'), a=params[0], v=params[1], ori=True)
print('Getting ready for building a tree')
if not os.path.exists(str(BayesPath(f'{dom}-{out}-bayes.nxs'))):
    SeqUtil.bayesfile(AlignsPath(f'{dom}-{out}.best.nex'), original_models,
                      BayesPath(f'{dom}-{out}-bayes.nxs'))
if not os.path.exists(str(BayesPath(f'{dom}-{out}-bayes.nxs.con.tre'))):
    run_bayes(str(BayesPath(f'{dom}-{out}-bayes.nxs')))
print('Tree built. Generating report.')
Report.generate_report(out, query, original_models, dom)

import os
import sys

import Report
from helpers.commands import run_bayes, run_phyml, run_prank, run_prottest
from helpers.constants import AlignsPath, BayesPath, DataPath, MLPath, ProtPath
from Utils import SeqUtil

out = sys.argv[1]
query = sys.argv[2]
dom = sys.argv[3]
paml = sys.argv[4].strip().lower() in {"1", "true", "yes", "y"}

print("Beginning alignment")
dom_out = f"{dom}-{out}"
SeqUtil.rename_seqs(DataPath(f"{dom_out}.fas"))
if not os.path.exists(str(AlignsPath(f"{dom_out}.best.nex"))):
    print("Aligning sequences")
    run_prank(str(DataPath(dom_out)), str(AlignsPath(dom_out)))
    SeqUtil.bayes_in_nex(AlignsPath(f"{dom_out}.best.nex"))
print("Alignment complete.\nCalculating best model for tree finding")
if not os.path.exists(str(ProtPath(f"{dom_out}.pro"))):
    run_prottest(str(AlignsPath(f"{dom_out}.best.nex")), str(ProtPath(f"{dom_out}.pro")))
original_models = SeqUtil.best_model(str(ProtPath(f"{dom_out}.pro")))
print(original_models)
if paml:
    print("pamling")
    for mod, params in original_models.items():
        print(mod, params)
        prefix = f"{dom_out}-{mod.split('+')[0]}-ori"
        if not os.path.exists(str(MLPath(f"{prefix}_phyml_boot_trees.txt"))):
            print("Path doesn't exist", str(MLPath(f"{prefix}_phyml_boot_trees.txt")))
            SeqUtil.nexus_to_proml(AlignsPath(f"{dom_out}.best.nex"), MLPath(prefix))
            if params[0] == "0" and params[1] == "0":
                print("No extra parameters")
                run_phyml(MLPath(prefix), mod.split("+")[0], AlignsPath(f"{dom_out}.best.dnd"))
            elif params[0] == "0":
                print("adding proportion of invariable sites")
                run_phyml(MLPath(prefix), mod.split("+")[0], AlignsPath(f"{dom_out}.best.dnd"), v=params[1])
            elif params[1] == "0":
                print("Adding gamma shape")
                run_phyml(MLPath(prefix), mod.split("+")[0], AlignsPath(f"{dom_out}.best.dnd"), a=params[0])
            else:
                print("Adding both proportion of invariable sites and gamma shape")
                run_phyml(
                    MLPath(prefix), mod.split("+")[0], AlignsPath(f"{dom_out}.best.dnd"), a=params[0], v=params[1]
                )
    print("Done Pamling")
print("Getting ready for building a tree")
if not os.path.exists(str(BayesPath(f"{dom_out}-bayes.nxs"))):
    SeqUtil.bayesfile(AlignsPath(f"{dom_out}.best.nex"), original_models, BayesPath(f"{dom_out}-bayes.nxs"))
if not os.path.exists(str(BayesPath(f"{dom_out}-bayes.nxs.con.tre"))):
    run_bayes(str(BayesPath(f"{dom_out}-bayes.nxs")))
print("Tree built. Generating report.")
Report.generate_report(out, query, original_models, dom)

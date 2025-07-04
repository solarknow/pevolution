import os

from Report import Report
from Utils import SeqUtil
from helpers.commands import run_prank, run_bayes, run_phyml, run_prottest
from helpers.constants import DataPath, AlignsPath, BayesPath, MLPath, ProtPath


class DomainRun:
    def __init__(self, out, query, dom, paml):
        self.out_prefix = out
        self.query_acc_number = query
        self.domain = dom
        self.is_paml = paml
        self.domain_out = f'{dom}-{out}'

    def run_alignment(self):
        print("Beginning alignment")
        SeqUtil.rename_seqs(DataPath(f'{self.domain_out}.fas'))
        if not os.path.exists(str(AlignsPath(f'{self.domain_out}.best.nex'))):
            print('Aligning sequences')
            run_prank(str(DataPath(self.domain_out)), str(AlignsPath(self.domain_out)))
            SeqUtil.make_bayes_compatible(AlignsPath(f'{self.domain_out}.best.nex'))
        print("Alignment complete.")

    def find_best_model(self):
        print('Calculating best model for tree finding')
        if not os.path.exists(str(ProtPath(f'{self.domain_out}.pro'))):
            run_prottest(str(AlignsPath(f'{self.domain_out}.best.nex')), str(ProtPath(f'{self.domain_out}.pro')))
        original_models = SeqUtil.best_protein_model(str(ProtPath(f'{self.domain_out}.pro')))
        print(original_models)
        return original_models

    def run_paml(self, models):
        if self.is_paml:
            print('pamling')
            for model in models:
                print(model)
                prefix = f"{self.domain_out}-{model.name.split('+')[0]}-ori"
                if not os.path.exists(str(MLPath(f"{prefix}_phyml_boot_trees.txt"))):
                    print("Path doesn't exist", str(MLPath(f"{prefix}_phyml_boot_trees.txt")))
                    SeqUtil.nexus_to_proml(AlignsPath(f'{self.domain_out}.best.nex'), MLPath(prefix))
                    if model.gamma == '0' and model.proportion == '0':
                        print('No extra parameters')
                        run_phyml(MLPath(prefix), model.name.split('+')[0], AlignsPath(f'{self.domain_out}.best.dnd'))
                    elif model.gamma == '0':
                        print('adding proportion of invariable sites')
                        run_phyml(MLPath(prefix), model.name.split('+')[0],
                                  AlignsPath(f'{self.domain_out}.best.dnd'), v=model.proportion)
                    elif model.proportion == '0':
                        print('Adding gamma shape')
                        run_phyml(MLPath(prefix), model.name.split('+')[0],
                                  AlignsPath(f'{self.domain_out}.best.dnd'), a=model.gamma)
                    else:
                        print("Adding both proportion of invariable sites and gamma shape")
                        run_phyml(MLPath(prefix), model.name.split('+')[0], AlignsPath(f'{self.domain_out}.best.dnd'),
                                  a=model.gamma, v=model.proportion)
            print("Done Pamling")
        else:
            print('Paml was not selected')

    def prep_and_run_bayes(self, models):
        print('Getting ready for building a tree')
        if not os.path.exists(str(BayesPath(f'{self.domain_out}-bayes.nxs'))):
            SeqUtil.create_mrbayes_cmd_file(AlignsPath(f'{self.domain_out}.best.nex'),
                                            models,
                                            BayesPath(f'{self.domain_out}-bayes.nxs'),
                                            )
        if not os.path.exists(str(BayesPath(f'{self.domain_out}-bayes.nxs.con.tre'))):
            run_bayes(str(BayesPath(f'{self.domain_out}-bayes.nxs')))
        print('Tree built. Generating report.')

    def run_report(self):
        self.run_alignment()
        models = self.find_best_model()
        self.run_paml(models)
        self.prep_and_run_bayes(models)
        Report(self.out_prefix, self.query_acc_number, models, self.domain).generate_report()

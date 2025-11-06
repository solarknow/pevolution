import getopt
import os
import sys

import Reciprocal
from Utils import SeqUtil, FetchUtil, FileUtil
from helpers.commands import run_domain_file
from helpers.constants import DataPath, BlastThresholds


def main(argv):
    # Expected: 1=query accession no.; 2=out prefix; 3=domain (euk,bac,arch,all)
    query = ""
    out = ""
    dom = ""
    phy = False

    try:
        opts, args = getopt.getopt(argv, "q:o:d:e:yh", ["query=", "output=", "domain=", "email=", "with_phyml"])
    except getopt.GetoptError:
        print(
            "Main.py -q|--query <query accession number> -o|--output <output name> -d|--domain <domain of life to "
            "query> -e|--email <email of user> [-y | --with_phyml]"
        )
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print(
                "Main.py -q|--query <query accession number> -o|--output <output name> -d|--domain <domain of life "
                "to query> -e|--email <email of user> [-y | --with_phyml]"
            )
            sys.exit()
        elif opt in ("-q", "--query"):
            query = arg
        elif opt in ("-o", "--output"):
            out = arg
        elif opt in ("-d", "--domain"):
            dom = arg
        elif opt in ("-y", "--with_phyml"):
            phy = True
        elif opt in ("-e", "--email"):
            FetchUtil.set_email(arg)
    if not os.path.exists(str(DataPath(dom + "-" + out + ".fas"))):
        arch_list = [
            "Haloferax volcanii",
            "Sulfolobus tokodaii",
            "Methanococcus aeolicus",
            "Methanobrevibacter smithii",
            "Thermococcus sibiricus",
            "Archaeoglobus fulgidus",
            "Nanoarchaeum equitans",
            "Thermoplasma acidophilum",
        ]
        bac_list = [
            "Gemmata obscuriglobus",
            "Prosthecobacter dejongeii",
            "Verrucomicrobium spinosum",
            "Rickettsia prowazekii",
            "Agrobacterium tumefaciens",
            "Escherichia coli",
            "Bacillus subtilis",
            "Anabaena variabilis",
            "Thermotoga maritima",
        ]
        euk_list = [
            "Drosophila melanogaster",
            "Homo sapiens",
            "Oryza sativa",
            "Trypanosoma brucei",
            "Plasmodium falciparum",
            "Saccharomyces cerevisiae",
            "Neurospora crassa",
            "Arabidopsis thaliana",
        ]
        # subject to change
        # setting threshold values: arch_thresh-w/ arch ;bac_thresh-w/ bac;
        dom_query = FetchUtil.fetch_organism(query)[1]
        if dom_query == "Archaea":
            thresh = BlastThresholds(arch=1e-10, bac=1e-5, euk=5)
        elif dom_query == "Eukaryota":
            thresh = BlastThresholds(arch=5, bac=5, euk=1e-10)
        else:
            thresh = BlastThresholds(arch=1e-5, bac=1e-10, euk=5)

        arch_accs = {}
        bac_accs = {}
        euk_accs = {}
        print("Blasting")
        if dom == "arch" or dom == "all":
            for a in arch_list:
                arch_accs.update(Reciprocal.best_reciprocal_blast(a, query, thresh.arch))
        if dom == "bac" or dom == "all":
            for b in bac_list:
                bac_accs.update(Reciprocal.best_reciprocal_blast(b, query, thresh.bac))
        if dom == "euk" or dom == "all":
            for e in euk_list:
                euk_accs.update(Reciprocal.best_reciprocal_blast(e, query, thresh.euk))

        all_accs = dict(list(arch_accs.items()) + list(bac_accs.items()) + list(euk_accs.items()))
        if dom == "all":
            num_seqs = sum(len(v) for v in all_accs.values())
            print(f"Dictionary generated with {len(all_accs)} keys and {num_seqs} sequences.")
        # Fetching the sequences and writing them to file
        print("Writing seqs to file.")
        if arch_accs:
            FileUtil.merge_domain_fastas("arch-" + out + ".fas", arch_accs)
            SeqUtil.addseq(DataPath("all-" + out + ".fas"), DataPath("arch-" + out + ".fas"))
        if bac_accs:
            FileUtil.merge_domain_fastas("bac-" + out + ".fas", bac_accs)
            SeqUtil.addseq(DataPath("all-" + out + ".fas"), DataPath("bac-" + out + ".fas"))
        if euk_accs:
            FileUtil.merge_domain_fastas("euk-" + out + ".fas", euk_accs)
            SeqUtil.addseq(DataPath("all-" + out + ".fas"), DataPath("euk-" + out + ".fas"))

    run_domain_file(out, query, dom, phy)


if __name__ == "__main__":
    main(sys.argv[1:])

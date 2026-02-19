import sys

import Reciprocal
from helpers.constants import BlastThresholds
from Utils import FetchUtil

if len(sys.argv) > 1:
    query = sys.argv[1]
    FetchUtil.set_email(sys.argv[2])
else:
    query = input("Query: ")
    FetchUtil.set_email(input("Email: "))
dom = FetchUtil.fetch_organism(query)
if dom[1] == "Archaea":
    thresh = BlastThresholds(arch=1e-10, bac=1e-5, euk=5)
elif dom[1] == "Eukaryota":
    thresh = BlastThresholds(arch=5, bac=5, euk=1e-10)
else:
    thresh = BlastThresholds(arch=1e-5, bac=1e-10, euk=5)

orgs = ["Homo sapiens", "Bacteroidota bacterium", "Haloferax volcanii"]
org_dict = zip(orgs, [thresh.euk, thresh.bac, thresh.arch], strict=True)
init_acc = [Reciprocal.best_reciprocal_blast(k, query, v) for k, v in org_dict]
print(init_acc)

runs = []
for idx in range(len(init_acc)):
    e = init_acc[idx]
    if not e:
        continue
    print("Pass " + repr(idx))
    for o in orgs:
        runs.append(Reciprocal.best_reciprocal_blast(o, list(e.values())[0][0], 5))

print(runs)

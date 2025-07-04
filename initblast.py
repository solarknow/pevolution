import sys

from Utils import FetchUtil, Reciprocal
from helpers.constants import domain_thresholds, organism_list, Domains

if len(sys.argv) > 1:
    query = sys.argv[1]
    FetchUtil.set_email(sys.argv[2])
else:
    query = input('Query: ')
    FetchUtil.set_email(input('Email: '))
dom = FetchUtil.fetch_organism(query)
thresh = domain_thresholds(dom[0])

orgs = {
    organism_list(Domains.EUKARYOTA)[0],
    organism_list(Domains.BACTERIA)[0],
    organism_list(Domains.ARCHAEA)[0]
}
org_dict= zip(orgs, [thresh.euk, thresh.bac, thresh.arch])
init_acc = [Reciprocal.best_reciprocal_blast(k, query, v) for k,v in org_dict]
print(init_acc)

runs = []
for idx in range(len(init_acc)):
    e = init_acc[idx]
    if not e:
        continue
    print("Pass " + str(idx))
    for o in orgs:
        runs.append(Reciprocal.best_reciprocal_blast(o, list(e.values())[0][0], 5))

print(runs)

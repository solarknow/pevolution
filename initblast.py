import sys

import FetchUtil
import Reciprocal

if len(sys.argv) > 1:
    query = sys.argv[1]
    FetchUtil.set_email(input("Email: "))
else:
    query = input('Query: ')
    FetchUtil.set_email(input("Email: "))
dom = FetchUtil.fetch_protein(query).domain
print("Source domain: " + str(dom))
if dom == 'Archaea':
    arch_e_thresh = 1e-10
    bac_e_thresh = 1e-5
    euk_e_thresh = 5
elif dom == 'Eukaryota':
    arch_e_thresh = 5
    bac_e_thresh = 5
    euk_e_thresh = 1e-10
else:
    arch_e_thresh = 1e-5
    bac_e_thresh = 1e-10
    euk_e_thresh = 5

init_acc = [Reciprocal.bestrecipblast('Homo sapiens', query, euk_e_thresh),
            Reciprocal.bestrecipblast('Escherichia coli', query, bac_e_thresh),
            Reciprocal.bestrecipblast('Haloferax volcanii', query, arch_e_thresh)]
runs = []
count = 0
orgs = ['Homo sapiens', 'Escherichia coli', 'Haloferax volcanii']
for acc in init_acc:
    if acc == {}:
        continue
    count += 1
    print("Pass " + repr(count))
    acs = [Reciprocal.bestrecipblast(o, list(acc.values())[0][0], 5) for o in orgs]
    runs.append(acs)

print(runs)

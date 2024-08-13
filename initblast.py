import sys
import FetchUtil
import Reciprocal
from helpers.constants import BlastThresholds

if len(sys.argv) > 1:
    query = sys.argv[1]
    FetchUtil.set_email(sys.argv[2])
else:
    query = input('Query: ')
    FetchUtil.set_email(input('Email: '))
dom = FetchUtil.fetch_organism(query)
if dom[1] == 'Archaea':
    thresh = BlastThresholds(arch=1e-10, bac=1e-5, euk=5)
elif dom[1] == 'Eukaryota':
    thresh = BlastThresholds(arch=5, bac=5, euk=1e-10)
else:
    thresh = BlastThresholds(arch=1e-5, bac=1e-10, euk=5)

init_acc = [Reciprocal.best_reciprocal_blast('Homo sapiens', query, thresh.euk),
            Reciprocal.best_reciprocal_blast('Escherichia coli', query, thresh.bac),
            Reciprocal.best_reciprocal_blast('Haloferax volcanii', query, thresh.arch)]
print(init_acc)
runs = []
count = 0
orgs = ['Homo sapiens', 'Escherichia coli', 'Haloferax volcanii']
for e in init_acc:
    count += 1
    if not e:
        continue
    print("Pass " + repr(count))
    for o in orgs:
        runs.append(Reciprocal.best_reciprocal_blast(o, list(e.values())[0][0], 5))

print(runs)

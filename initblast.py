import sys
from multiprocessing import Process, Queue

import FetchUtil
import Reciprocal

try:
    query = sys.argv[1]
except IndexError:
    query = input('Query: ')
dom = FetchUtil.fetch_organism(query)
if dom == 'Archaea':
    thresh1 = 1e-10
    thresh2 = 1e-5
    thresh3 = 5
elif dom == 'Eukaryota':
    thresh1 = 5
    thresh2 = 5
    thresh3 = 1e-10
else:
    thresh1 = 1e-5
    thresh2 = 1e-10
    thresh3 = 5

init_acc = [Reciprocal.bestrecipblast('Homo sapiens', query, thresh3, None),
            Reciprocal.bestrecipblast('Escherichia coli', query, thresh2, None),
            Reciprocal.bestrecipblast('Haloferax volcanii', query, thresh1, None)]
runs = []
count = 0
orgs = ['Homo sapiens', 'Escherichia coli', 'Haloferax volcanii']
for e in init_acc:
    count += 1
    if e == {}:
        continue
    print("Pass " + repr(count))
    for o in orgs:
        q = Queue()
        p = Process(target=Reciprocal.bestrecipblast, args=(o, e.values()[0][0], 5, q))
        p.start()
        p.join()
        acs = []
        while not q.empty():
            acs.append(q.get())
        runs.append(acs)

print(runs)

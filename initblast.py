#!/usr/bin/env python
import Reciprocal, Fetchutil
import sys, subprocess
from multiprocessing import Process,Queue
if '-h' in sys.argv:
  print '''Syntax: ./initblast.py query [-local|-h]\n
	\t-h\tThis help screen\n
	\t-local\tRun this script on local BLAST db's'''
  sys.exit()
try:
    query=sys.argv[1]
except:
    query=raw_input('Query: ')
local='-local' in sys.argv
orgs=['Homo sapiens','Escherichia coli','Haloferax volcanii']
if local:
  params=['python','Local_Ebot.py']+orgs
  subprocess.call(params)
dom=Fetchutil.orgfetch(query,local)
if dom=='Archaea':
    thresh1=1e-10
    thresh2=1e-5
    thresh3=5
elif dom=='Eukaryota':
    thresh1=5
    thresh2=5
    thresh3=1e-10
else:
    thresh1=1e-5
    thresh2=1e-10
    thresh3=5


init_acc=[]
init_acc.append(Reciprocal.bestrecipblast('Homo sapiens', query, thresh3,None,local))
init_acc.append(Reciprocal.bestrecipblast('Escherichia coli', query, thresh2,None,local))
init_acc.append(Reciprocal.bestrecipblast('Haloferax volcanii', query, thresh1,None,local))
runs=[]
count=0

for e in init_acc:
    count+=1
    if e=={}:
        continue
    print "Pass "+repr(count)
    for o in orgs:
      q=Queue()
      p=Process(target=Reciprocal.bestrecipblast, args=(o, e.values()[0][0], 5, q, local))
      p.start()
      p.join()
      acs=[]
      while not q.empty():
        acs.append(q.get())
      runs.append(acs)

print runs

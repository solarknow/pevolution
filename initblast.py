import Reciprocal, Fetchutil
import sys

try:
    query=sys.argv[1]
except:
    query=raw_input('Query: ')
dom=Fetchutil.orgfetch(query)
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
init_acc.append(Reciprocal.bestrecipblast(['Homo sapiens'], query, thresh3))
init_acc.append(Reciprocal.bestrecipblast(['Escherichia coli'], query, thresh2))
init_acc.append(Reciprocal.bestrecipblast(['Haloferax volcanii'], query, thresh1))
runs=[]
count=0
for e in init_acc:
    count+=1
    if e=={}:
        continue
    print "Pass "+repr(count)
    res=Reciprocal.bestrecipblast(['Homo sapiens','Escherichia coli','Haloferax volcanii'], e.values()[0][0], 5)
    acs=[]
    for l in res:
        acs.append(res[l][0])
    runs.append(acs)

print runs

import recip_edit, Fetchutil
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
init_acc.append(recip_edit.bestrecipblast(['Homo sapiens'], query, thresh3))
init_acc.append(recip_edit.bestrecipblast(['Escherichia coli'], query, thresh2))
init_acc.append(recip_edit.bestrecipblast(['Haloferax volcanii'], query, thresh1))
runs=[]
for e in init_acc:
    if e=={}:
        continue
    print "Pass 2"
    res=recip_edit.bestrecipblast(['Homo sapiens','Escherichia coli','Haloferax volcanii'], e.values()[0][0], 5)
    acs=[]
    for l in res:
        acs.append(res[l][0])
    runs.append(acs)

print runs

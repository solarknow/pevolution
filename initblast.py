"""
Script to test "blastability" of a given sequence
Usage:
    Main.py [options...] --email <email>
    -e, --email <email> Sets your email for remote blasting
    -h, --help This help text
        --local Sets the blast to run locally, currently requires a blastdb called 'nr'
    -q, --query <sequence id> An NCBI protein sequence id
"""
import getopt
import sys

from Bio import Entrez

import FetchUtil
import Reciprocal

try:
    options, remainders = getopt.getopt(sys.argv[1:], 'q:e:h', ['query=', 'email=', 'help', 'local'])
except getopt.GetoptError as err:
    # print help information and exit:
    print(str(err))  # will print something like "option -a not recognized"
    print(__doc__)
    sys.exit(2)
query = ''
out = ''
dom = ''
local = False
for opt, arg in options:
    if opt in ('-q', '--query'):
        query = arg
    elif opt in ('-e', '--email'):
        FetchUtil.set_email(arg)
    elif opt in ('-h', '--help'):
        print(__doc__)
        sys.exit(2)
    elif opt == '--local':
        local = True
# if len(sys.argv) > 1:
#     query = sys.argv[1]
#     local = False if sys.argv[2] == '--remote' else True
#     if Entrez.email is None:
#         FetchUtil.set_email(input("Email: "))
# else:
#     query = input('Query: ')
#     local = True if input('Local? [y/n] ') == 'y' else False
#     if Entrez.email is None:
#         FetchUtil.set_email(input("Email: "))
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

init_acc = [Reciprocal.bestrecipblast([9606, 'Homo sapiens'], query, euk_e_thresh, local),
            Reciprocal.bestrecipblast([358, 'Agrobacterium tumefaciens'], query, bac_e_thresh, local),
            Reciprocal.bestrecipblast([309800, 'Haloferax volcanii'], query, arch_e_thresh, local)]
runs = []
count = 0
orgs = [[9606, 'Homo sapiens'], [358, 'Agrobacterium tumefaciens'], [309800, 'Haloferax volcanii']]
for acc in init_acc:
    if acc == {}:
        continue
    count += 1
    print("Pass " + repr(count))
    acs = [Reciprocal.bestrecipblast(o, list(acc.values())[0][0], 5, True) for o in orgs]
    runs.append(acs)

print(runs)

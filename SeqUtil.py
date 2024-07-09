import ast
import os
import random

import psutil

from helpers.commands import run_prottest, clustal_align
from helpers.file_formats import nexus_fmt

bayesmodels = ['poisson', 'jtt', 'mtrev', 'mtmam', 'wag', 'rtrev', 'cprev', 'vt', 'blosum', 'dayhoff']


def find_longest_key_length(dicto):
    """finds the length of the longest key in dicto"""
    pivot = ''
    for k in dicto:
        if len(k) > len(pivot):
            pivot = k
    return pivot, len(pivot)


def clustal_to_fasta(clust, out):
    """converts a clustal alignment file to a fasta file out"""
    with open(clust) as hand:
        seqs = {}
        hand.readline()
        hand.readline()
        hand.readline()
        while hand:
            line = hand.readline()
            spl=line.split()
            if spl == [] and len(line) > 0:
                continue
            if line == '':
                break
            elif spl[0][0] in [':', '.', '*']:
                hand.readline()
                continue
            seqs[spl[0]] = seqs.get(spl[0], '') + spl[1]
        for j in seqs:
            with open(out, 'a') as outfile:
                outfile.write('>' + j + '\n' + seqs[j] + '\n\n')


def read_dict(fil):
    """Reads in a dict object from a given file"""
    with open(fil) as p:
        return ast.literal_eval(p.read())


def rename_seqs(fas, out=None):
    """Takes in a fasta file fas that has sequences and shortens the names, assigning to a new file"""
    if not out:
        out = fas.split('.')[0]
    orgs = {}
    with open(fas) as infile:
        with open(out, 'w') as outfile:
            while infile:
                lin = infile.readline()
                if lin == '':
                    break
                elif lin[0] == '>':
                    outfile.write('>')
                    linstr = lin[1:].split()
                    org = linstr[0][0] + linstr[1][0:3]
                    if org in orgs:
                        org += repr(orgs[org] + 1)
                    else:
                        orgs[org] = 1
                    outfile.write(org + '\n')
                else:
                    outfile.write(lin)


def remove_gaps_nexus(fil):
    """Modifies the file of nexus aligned sequences to remove gapped positions"""
    arr = {}
    with open(fil) as infile:
        while infile:
            title = infile.readline()
            if title.startswith('matrix'):
                while infile:
                    title = infile.readline()
                    if title[0] == '\n':
                        continue
                    elif title[0] == ';':
                        break
                    spl = title.split()
                    arr[spl[0]] = arr.get(spl[0], []) + spl[1].strip().split('-')
                break
    return {k: ''.join(v) for k, v in arr.items()}


def nexus_to_proml(seqs, inpath):
    """Converts aligned Nexus file seqs to a ProML input file, inpath"""
    with open(seqs) as seqsin:
        count = 0
        length = 0
        while seqsin:
            lin = seqsin.readline()
            if lin.startswith('dimensions'):
                count = int(lin.split()[1][5:])
                length = int(lin.split()[2][6:len(lin.split()[2]) - 1])
                break

        while seqsin:
            lin = seqsin.readline()
            if lin.startswith('matrix'):
                break
        dicto = {}
        while seqsin:
            lin = seqsin.readline()
            if lin == '\n':
                continue
            elif lin == ';\n':
                break
            name = lin.split()[0]
            seq = lin.split()[1].strip()
            dicto[name] = dicto.get(name, '') + seq
    if not length == len(list(dicto.values())[0]):
        length = len(list(dicto.values())[0])

    with open(inpath, 'w') as out:
        out.write(repr(count) + '  ' + repr(length) + '\n')
        for k in dicto:
            out.write(k)
            for i in range(30 - len(k)):
                out.write(' ')
            out.write(dicto[k] + '\n')


def bayes_in_nex(infile):
    """Modifies a PRANK alignment file to be compatible with MrBayes"""
    with open(infile) as inf:
        lines = inf.readlines()

    for i in range(len(lines)):
        if 'trees' in lines[i]:
            lines = lines[:i]
            break
        else:
            lines[i] = lines[i].replace('\'', '')
    with open(infile, 'w') as fil:
        for j in lines:
            fil.write(j)


def bayesfile(infile, model, outfile):
    """Writes a Nexus file for use as a MrBayes batch file"""
    os.makedirs('Bayes', exist_ok=True)
    for k in model:
        extra = k.split('+')
        if extra[0].lower() == 'jtt':
            extra[0] = 'jones'
        elif extra[0].lower() == 'blosum62':
            extra[0] = 'blosum'
        if extra[0].lower() in bayesmodels:
            with open(outfile, 'w') as han:
                han.write('#NEXUS\n' +
                          'begin mrbayes;\n' +
                          f'\texe {infile};\n' +
                          f'\tprset aamodelpr=fixed({extra[0]});\n')
                if len(extra) > 1:
                    if 'I' in extra and 'G' in extra:
                        han.write('\tlset rates=Invgamma;\n')
                        han.write(f'\tprset shapepr=fixed({model[k][0]});\n')
                    elif 'I' in extra:
                        han.write('\tlset rates=Propinv;\n')
                    elif 'G' in extra:
                        han.write('\tlset rates=Gamma;\n')
                        han.write(f'\tprset shapepr=fixed({model[k][0]});\n')

                han.write(f'\tmcmc ngen=50000 samplefreq=50 file={outfile};\n' +
                          '\tsumt burnin=250;\n' +
                          'end;\n\n')


def boot(infile, outfile, norep):
    """takes in an aligned phylip file and creates a file in phylip format with norep bootstrap datasets"""
    with open(infile) as hand:
        lin = hand.readline()
        num_seqs = int(lin.split()[0])
        length = int(lin.strip().split()[1])
        seqs = {}
        for i in range(num_seqs):
            lin = hand.readline().strip().split()
            seqs[lin[0]] = lin[1]
    with open(outfile, 'w') as hand2:
        for j in range(norep):
            hand2.write(str(num_seqs) + '   ' + str(length) + '\n')
            newseq = {}
            for i in range(length):
                rand_int = random.randint(0, length - 1)
                for k in seqs:
                    newseq[k] = newseq.get(k, '') + seqs[k][rand_int]
            for k in newseq:
                hand2.write(k + '    ' + newseq[k] + '\n')
            hand2.write('\n\n')


def count_fasta_seqs(fil):
    """returns number of seqs in given fil fasta file"""
    try:
        p = open(fil)
    except FileNotFoundError:
        return 0
    num_seqs = sum([1 for lin in p if lin.startswith('>')])
    p.close()
    return num_seqs


def splice_align(inalign, outalign):
    """Splices out evolutionarily uncertain areas in inalign, realigns the new sequence"""
    align = open(inalign)
    dicto = {}
    inalign = inalign.split('.')[0]
    noseq = 0
    align_size = 0
    while align:
        lin = align.readline()
        if lin.startswith('end'):
            break
        if lin.startswith('dimensions'):
            li = lin.strip().split()
            for i in li:
                if i.startswith('ntax'):
                    noseq = int(i.split('=')[1])
                elif i.startswith('nchar'):
                    align_size = int(i.split('=')[1].split(';')[0])
        if lin.startswith('matrix'):
            while align:
                lin = align.readline()
                if lin == '\n':
                    continue
                elif lin == ';\n':
                    break
                else:
                    lo = lin.split()
                    try:
                        dicto[lo[0]] += lo[1]
                    except KeyError:
                        dicto.update({lo[0]: lo[1]})
    # sequences have been read in; determining the positions to be removed
    pos = {}

    for i in dicto.values():
        for c in range(len(i)):
            if i[c] == '-':
                try:
                    pos[c] += 1
                except KeyError:
                    pos.update({c: 1})
    # print pos
    # processing pos; eliminating key:value pairs whose values are less than half len(i)
    for k in pos.keys():
        if pos[k] < (0.5 * noseq):
            pos.pop(k)
    post = list(pos.keys())
    post.sort()
    # initial condition to stop recursion
    if len(post) <= .05 * align_size:
        with open(outalign, 'w') as han:
            han.write(nexus_fmt(repr(noseq), repr(align_size), 'protein'))

            for k in dicto.keys():
                han.write(k)
                for i in range((find_longest_key_length(dicto) + 6) - len(k)):
                    han.write(' ')
                han.write(dicto[k] + '\n')
            han.write(';\nend;\n')
        return
    # adds one residue on either side of the residues indicated in post
    for p in range(len(post)):
        if 0 < post[p] < align_size - 1:
            bef = post[p] - 1
            af = post[p] + 1
            if bef in post:
                post.append(bef)
            if af in post:
                post.append(af)
    post.sort()
    post.reverse()

    print(align_size, len(post))

    # processes dicto; removes residues contained in post
    for i in dicto.keys():
        val = dicto[i]
        new_val = ''
        for j in range(len(post) - 1):
            new_val += val[post[j + 1] + 1:post[j]]
        for k in range(align_size - len(new_val)):
            new_val += '-'
        dicto.update({i: new_val})
    # creates the nex file containing the newly spliced data
    align_size = len(list(dicto.values())[0])
    with open(inalign + '.edit', 'w') as han:
        han.write(nexus_fmt(repr(len(dicto.keys())), repr(align_size), 'protein'))

        for k in dicto.keys():
            han.write(k)
            for i in range((find_longest_key_length(dicto) + 6) - len(k)):
                han.write(' ')
            han.write(dicto[k])
            # print dicto[k]
            han.write('\n')
        han.write(';\nend;\n')
    print('.edit written. aligning')
    clustal_align(inalign + '.edit', inalign + '.ed', "nexus")
    print('.ed written splice aligned')


def bestmod(infile):
    """Takes in alignment file runs protTest, and extracts best model(s)"""
    if not os.path.exists('Prot'):
        os.mkdir('Prot')
    procs = int(round(psutil.cpu_count() / 2.0))
    out = infile.split(os.sep)[1].split('.')[0]
    outfile = 'Prot' + os.sep + out + '.pro'
    run_prottest(infile, outfile, repr(procs))
    with open(outfile) as prot_hand:
        models = {}
        ret = {}
        # reading and processing prottest output
        while prot_hand:
            lin = prot_hand.readline()
            if lin.startswith('Model.'):
                lsplit = lin.split()
                mod = lsplit[2]
                modmod = mod.split('+')
                para = ['0', '0']  # index 0 is G index 1 is I
                if len(modmod) > 1:
                    if 'G' in modmod:
                        while 1:
                            lin = prot_hand.readline()
                            # print lin,'g'
                            if lin.split()[0] == 'gamma':
                                para[0] = lin.split()[6]
                                break
                    if 'I' in modmod:
                        while 1:
                            lin = prot_hand.readline()
                            # print lin.split(),'i'
                            if lin.split()[0] == 'proportion':
                                para[1] = lin.split()[5]
                                break
                models.update({mod: para})
                # print models
                continue
            elif lin.startswith('Best model'):
                lsplit = lin.split()
                mod = lsplit[5]
                ret.update({mod: models[mod]})
                # print ret
                if mod.lower().split('+')[0] not in bayesmodels:
                    while 1:
                        lin = prot_hand.readline()
                        # print lin,1
                        if lin.endswith('-\n'):
                            while 2:
                                lin = prot_hand.readline().split()
                                # print lin,2
                                if lin[0].split('+')[0].lower() in bayesmodels:
                                    ret.update({lin[0]: models[lin[0]]})
                                    # print ret
                                    break
                            break

                # elif lin.startswith('*') and :
                #    print "Saving Tree"
                #    tree=''
                #   while 9:
                #       lin=prot_hand.readline().strip()
                #      tree+=lin
                #      if lin.endswith(';'):
                #         break
                # open('Prot/'+out+'.tre','w').write(tree)
                return ret


def addseq(oldseq, newseq):
    """Transfers the sequence from newseq to oldseq"""
    old = open(oldseq, 'a')
    new = open(newseq)
    old.write(new.read())
    old.close()
    new.close()

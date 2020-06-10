import json
import os
import psutil
import random
import subprocess

from constants import BAYES_PATH, ALIGNS_PATH, PROT_PATH, PROTTEST_PATH

bayesmodels = ['poisson', 'jtt', 'mtrev', 'mtmam', 'wag', 'rtrev', 'cprev', 'vt', 'blosum', 'dayhoff']
os.makedirs(ALIGNS_PATH, exist_ok=True)
os.makedirs(BAYES_PATH, exist_ok=True)


def longest_key_length(dicto):
    """finds the length of the longest key in dicto"""
    keys = dicto.keys()
    pivot = ''
    for k in keys:
        if len(k) > len(pivot):
            pivot = k
    return len(pivot)


def clustal_to_fasta(clust, out):
    """converts a clustal alignment file to a fasta file out"""
    with open(clust) as hand:
        seqs = {}
        hand.readline()
        hand.readline()
        hand.readline()
        while hand:
            line = hand.readline()
            if line.split() == [] and len(line) > 0:
                continue
            if line == '':
                break
            elif line.split()[0][0] in [':', '.', '*']:
                hand.readline()
                continue
            lin = line.split()
            try:
                seqs[lin[0]] += lin[1]
            except KeyError:
                seqs.update({lin[0]: lin[1]})
        for j in seqs.keys():
            with open(out, 'a') as outfile:
                outfile.write('>' + j + '\n' + seqs[j] + '\n\n')


def dict_extract(fil):
    """Extracts a dict object from a given file"""
    with open(fil) as p:
        dic = json.load(p)
    return dic


def shorten_sequence_names(fas):
    """Takes in a fasta file fas that has sequences and shortens the names, assigning to a new file"""
    out = fas.split('.')[0]

    orgs = {}
    with open(fas) as infile:
        for line in infile:
            with open(out, 'a') as outfile:
                if line.startswith('>'):
                    linstr = line[1:].split()
                    org = linstr[0][0] + linstr[1][0:3]
                    if org in orgs:
                        org += repr(orgs[org] + 1)
                    else:
                        orgs[org] = 1
                    outfile.write('>' + org + '\n')
                else:
                    outfile.write(line + '\n')


def remove_gaps(fil):
    """Modifies the file of aligned sequences to remove gapped positions"""
    arr = {}
    infile = open(fil)
    while infile:
        title = infile.readline()
        if title.startswith('matrix'):
            while infile:
                title = infile.readline()
                if title.startswith('\n'):
                    continue
                elif title.startswith(';'):
                    break
                spl = title.split()
                if spl[0] in arr:
                    arr[spl[0]] += spl[1].strip().split('-')
                else:
                    arr[spl[0]] = spl[1].strip().split('-')
            break
    return {k: ''.join(v) for k, v in arr.items()}


def fasta_to_proml(seqs, inpath):
    """Converts aligned FASTA file seqs to a ProML input file, inpath, Plus 1 space"""
    seqstr = ''
    with open(seqs) as seqsin:
        seqstr += seqsin.readline()[1:].strip()
        count = 1
        length = 0
        for i in range(11 - len(seqstr)):
            seqstr += ' '
        while seqsin:
            line = seqsin.readline()
            if line.startswith('>'):
                count += 1
                seqstr += '\n' + line[1:].strip()
                for i in range(11 - len(line[1:].strip())):
                    seqstr += ' '
            else:
                seqstr += line.strip()
                if count == 1:
                    length += len(line.strip())
            if line == '':
                break
    with open(inpath, 'w') as seqsout:
        seqsout.write(str(count) + '  ' + str(length) + '\n' + seqstr)


def nexus_to_proml(seqs, inpath):
    """Converts aligned Nexus file seqs to a ProML input file, inpath"""
    count = 0
    length = 0
    with open(seqs) as seqsin:
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
            if name in dicto:
                dicto[name] += seq
            else:
                dicto[name] = seq
    if not length == len(list(dicto.values())[0]):
        length = len(list(dicto.values())[0])
    with open(inpath, 'w') as out:
        out.write(str(count) + '  ' + str(length) + '\n')
        for k in dicto:
            out.write(k)
            out.write(' ' * (30 - len(k)))
            out.write(dicto[k] + '\n')


def prank_to_mrbayes(infile):
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


def write_to_mrbayes(infile, model, outfile):
    """Writes a Nexus file for use as a MrBayes batch file"""
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
                          '\texe ' + infile + ';\n' +
                          '\tprset aamodelpr=fixed(' + extra[0] + ');\n')
                if len(extra) > 1:
                    if 'I' in extra and 'G' in extra:
                        han.write('\tlset rates=Invgamma;\n')
                        han.write('\tprset shapepr=fixed(' + model[k][0] + ');\n')
                    elif 'I' in extra:
                        han.write('\tlset rates=Propinv;\n')
                    elif 'G' in extra:
                        han.write('\tlset rates=Gamma;\n')
                        han.write('\tprset shapepr=fixed(' + model[k][0] + ');\n')

                han.write('\tmcmc ngen=50000 samplefreq=50 file=' + outfile + ';\n' +
                          '\tsumt burnin=250;\n' +
                          'end;\n\n')


def create_bootstrap(infile, outfile, reps):
    """takes in an aligned phylip file and creates a file in phylip format with norep bootstrap datasets"""
    with open(infile) as hand:
        lin = hand.readline()
        noseqs = int(lin.split()[0])
        length = int(lin.strip().split()[1])
        seqs = {}
        for i in range(noseqs):
            lin = hand.readline().strip().split()
            seqs.update({lin[0]: lin[1]})
    with open(outfile, 'w') as hand2:
        for j in range(reps):
            hand2.write(str(noseqs) + '   ' + str(length) + '\n')
            newseqs = {}
            for _ in range(length):
                rand_int = random.randint(0, length - 1)
                for k in seqs:
                    if k in newseqs:
                        newseqs[k] += seqs[k][rand_int]
                    else:
                        newseqs[k] = seqs[k][rand_int]
            for m in newseqs:
                hand2.write(m + '    ' + newseqs[m] + '\n')
            hand2.write('\n\n')


def sequences_count(fil):
    """returns number of seqs in given fil fasta file"""
    if not os.path.exists(fil):
        return 0
    count = 0
    with open(fil) as p:
        while p:
            lin = p.readline()
            if lin == '':
                break
            elif lin.startswith('>'):
                count += 1
    return count


def clean_and_realign(inalign, outalign):
    """Splices out evolutionarily uncertain areas in inalign, realigns the new sequence"""
    with open(inalign) as align:
        dicto = {}
        inalign = inalign.split('.')[0]
        noseq = 0
        align_size = 0
        while align:
            line = align.readline()
            if line.startswith('end'):
                break
            if line.startswith('dimensions'):
                li = line.strip().split()
                for i in li:
                    if i.startswith('ntax'):
                        noseq = int(i.split('=')[1])
                    elif i.startswith('nchar'):
                        align_size = int(i.split('=')[1].split(';')[0])
            if line.startswith('matrix'):
                while align:
                    line = align.readline()
                    if line == '\n':
                        continue
                    elif line == ';\n':
                        break
                    else:
                        lo = line.split()
                        if lo[0] in dicto:
                            dicto[lo[0]] += lo[1]
                        else:
                            dicto[lo[0]] = lo[1]
    # sequences have been read in; determining the positions to be removed
    positions = {}

    for i in dicto.values():
        for c in range(len(i)):
            if i[c] == '-':
                if c in positions:
                    positions[c] += 1
                else:
                    positions[c] = 1
    # processing positions; eliminating key:value pairs whose values are less than half len(i)
    for k in positions:
        if positions[k] < (0.5 * noseq):
            positions.pop(k)
    post = list(positions.keys())
    post.sort()
    # initial condition to stop recursion
    if len(post) <= .05 * align_size:
        with open(outalign, 'w') as han:
            han.write('#NEXUS\n' +
                      'begin data;\n' +
                      'dimensions ntax=' + repr(noseq) + ' nchar=' + repr(align_size) + ';\n' +
                      'format datatype=protein interleave=no gap=-;\n' +
                      'matrix\n\n')

            for k in dicto:
                han.write(k)
                for i in range((longest_key_length(dicto) + 6) - len(k)):
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
    for i in dicto:
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
        han.write('#NEXUS\n' +
                  'begin data;\n' +
                  'dimensions ntax=' + repr(len(dicto)) + ' nchar=' + repr(align_size) + ';\n' +
                  'format datatype=protein interleave=no gap=-;\n' +
                  'matrix\n\n')

        for k in dicto.keys():
            han.write(k)
            for i in range((longest_key_length(dicto) + 6) - len(k)):
                han.write(' ')
            han.write(dicto[k])
            han.write('\n')
            han.write(';\nend;\n')
    print('.edit written. aligning')
    os.system('clustalw -align -infile=' + inalign + '.edit -outfile=' + inalign + '.ed -output=nexus -quiet')
    print('.ed written splice aligning')


def prottest_best_models(infile):
    """Takes in alignment file runs protTest, and extracts best model(s)"""
    os.makedirs(PROT_PATH, exist_ok=True)
    procs = int(round(psutil.cpu_count() / 2.0))
    out = infile.split(os.sep)[1].split('.')[0]
    prot_cmd = ['java', '-jar', PROTTEST_PATH, '-i', infile, '-o', PROT_PATH + out + '.pro', '-all-distributions',
                '-all', '-S', 1, '-threads', procs, '-BIC']
    subprocess.run(prot_cmd)
    models = {}
    ret = {}
    # reading and processing prottest output
    with open(PROT_PATH + out + '.pro') as prot_hand:
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
                            if lin.split()[0] == 'proportion':
                                para[1] = lin.split()[5]
                                break
                models.update({mod: para})
                continue
            elif lin.startswith('Best model'):
                lsplit = lin.split()
                mod = lsplit[5]
                ret.update({mod: models[mod]})
                if mod.lower().split('+')[0] not in bayesmodels:
                    while 1:
                        lin = prot_hand.readline()
                        if lin.endswith('-\n'):
                            while 2:
                                lin = prot_hand.readline().split()
                                if lin[0].split('+')[0].lower() in bayesmodels:
                                    ret.update({lin[0]: models[lin[0]]})
                                    break
                            break
            return ret


def fasta_add_sequence(oldseq, newseq):
    """Transfers the sequence from newseq to oldseq"""
    with open(oldseq, 'a') as old:
        with open(newseq) as new:
            while new:
                add = new.readline()
                old.write(add)
                if add == '':
                    return

import ast
import random

from helpers.constants import nexus_fmt

bayesmodels = ['poisson', 'jtt', 'mtrev', 'mtmam', 'wag', 'rtrev', 'cprev', 'vt', 'blosum', 'dayhoff']


def find_longest_key_length(dicto):
    """finds the length of the longest key in dicto"""
    pivot = ''
    for k in dicto:
        if len(k) > len(pivot):
            pivot = k
    return pivot, len(pivot)


def rename_seqs(fas, out=None):
    """Takes in a fasta file path fas that has sequences and shortens the names, assigning to a new file"""
    if not out:
        out = str(fas).split('.')[0]
    orgs = {}
    with open(str(fas)) as infile:
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


def nexus_to_proml(seqs, inpath):
    """Converts aligned Nexus file seqs to a ProML input file, inpath"""
    with open(str(seqs)) as seqsin:
        count = 0
        length = 0
        while seqsin:
            lin = seqsin.readline()
            if lin.startswith('dimensions'):
                dims = lin.split()
                count = int(dims[1].split('=')[1])
                length = int(dims[2].split('=')[1][:-1])
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
            line_split = lin.split()
            name = line_split[0]
            seq = line_split[1].strip()
            dicto[name] = dicto.get(name, '') + seq
    if not length == len(list(dicto.values())[0]):
        length = len(list(dicto.values())[0])

    with open(str(inpath), 'w') as out:
        out.write(repr(count) + '  ' + repr(length) + '\n')
        for k in dicto:
            out.write(k)
            for i in range(30 - len(k)):
                out.write(' ')
            out.write(dicto[k] + '\n')


def bayes_in_nex(infile):
    """Modifies a PRANK alignment file to be compatible with MrBayes"""
    with open(str(infile)) as inf:
        lines = inf.readlines()

    for i in range(len(lines)):
        if 'trees' in lines[i]:
            lines = lines[:i]
            break
        else:
            lines[i] = lines[i].replace('\'', '')
    with open(str(infile), 'w') as fil:
        for j in lines:
            fil.write(j)


def bayesfile(infile, model, outfile):
    """Writes a Nexus file for use as a MrBayes batch file"""
    for k in model:
        extra = k.split('+')
        if extra[0].lower() == 'jtt':
            extra[0] = 'jones'
        elif extra[0].lower() == 'blosum62':
            extra[0] = 'blosum'
        if extra[0].lower() in bayesmodels:
            with open(str(outfile), 'w') as han:
                han.write('#NEXUS\n' +
                          'begin mrbayes;\n' +
                          f'\texe {str(infile)};\n' +
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

                han.write(f'\tmcmc ngen=50000 samplefreq=50 file={str(outfile)};\n' +
                          '\tsumt burnin=250;\n' +
                          'end;\n\n')


def best_model(outfile):
    """Takes in alignment file, runs protTest, and extracts best model(s)
    @returns {
    model: [gamma, proportion]
    }
    """
    with open(str(outfile)) as prot_hand:
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
                return ret


def addseq(oldseq, newseq):
    """Transfers the sequence from newseq to oldseq"""
    old = open(str(oldseq), 'a')
    new = open(str(newseq))
    old.write(new.read())
    old.close()
    new.close()